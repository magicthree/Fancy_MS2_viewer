#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# fancyMS2viewer - a viewer for MS/MS (MS2) spectral libraries.
# Copyright (C) 2026  magicthree
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
fancyMS2viewer - a viewer for .msp / .mgf / .mzML MS/MS spectral libraries.

Features
--------
- Load one or more .msp files (NIST-style text format).
- Browse all spectra in a list.
- Search / filter spectra by any metadata field (Name, Precursor m/z, etc.).
- Interactive stick-plot of the MS2 peaks:
    * hover to read m/z / intensity of the nearest peak
    * drag left-to-right to zoom the m/z range
    * mouse wheel to zoom in/out, right-click (or 'R') to reset the view
- Metadata panel showing every field of the selected spectrum.

No third-party dependencies: standard-library tkinter only.

Usage
-----
    py msp_viewer.py [file1.msp file2.msp ...]

Files can also be opened from the File menu after launch.
"""

import os
import re
import sys
import gzip
import base64
import zlib
import xml.etree.ElementTree as ET
from array import array
import tkinter as tk
from tkinter import ttk, filedialog

# Optional OS drag-and-drop support (drop files onto the window).
# Falls back gracefully to a plain Tk window if tkinterdnd2 is unavailable.
try:
    from tkinterdnd2 import TkinterDnD, DND_FILES
    _TkBase = TkinterDnD.Tk
    _DND_AVAILABLE = True
except Exception:  # noqa: BLE001
    _TkBase = tk.Tk
    DND_FILES = None
    _DND_AVAILABLE = False


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #

class _PeakView:
    """Sequence view over two parallel arrays, yielding (mz, intensity) tuples.

    Lets peaks be stored compactly as array('d') while all existing code that
    iterates/indexes/len()s a list of (mz, intensity) tuples keeps working.
    """

    __slots__ = ("mz", "it")

    def __init__(self, mz, it):
        self.mz = mz
        self.it = it

    def __len__(self):
        return len(self.mz)

    def __bool__(self):
        return len(self.mz) > 0

    def __iter__(self):
        return zip(self.mz, self.it)

    def __getitem__(self, i):
        if isinstance(i, slice):
            return list(zip(self.mz[i], self.it[i]))
        return (self.mz[i], self.it[i])


class Spectrum:
    """One MS/MS record. Peaks are stored as two compact array('d')."""

    __slots__ = ("meta", "mz", "inten", "source", "num", "_blob")

    def __init__(self, meta, peaks=None, source="", mz=None, inten=None):
        self.meta = meta            # dict: field name -> value (str)
        self.source = source        # originating file
        self.num = 0                # 1-based load order (set by the app)
        self._blob = None           # cached lower-cased search text
        if mz is not None:          # fast path: arrays supplied directly
            n = min(len(mz), len(inten))
            self.mz = mz if n == len(mz) else mz[:n]
            self.inten = inten if n == len(inten) else inten[:n]
        else:                       # from a list of (mz, intensity) tuples
            peaks = peaks or []
            self.mz = array("d", [p[0] for p in peaks])
            self.inten = array("d", [p[1] for p in peaks])

    @property
    def peaks(self):
        return _PeakView(self.mz, self.inten)

    # convenience accessors -------------------------------------------------
    @property
    def name(self):
        for key in ("Name", "NAME", "name", "Compound_name", "Title",
                    "TITLE", "title", "id"):
            if key in self.meta:
                return self.meta[key]
        return "(unnamed)"

    @property
    def precursor(self):
        for key in ("PrecursorMZ", "PRECURSORMZ", "Precursor_mz",
                    "precursor_mz", "PRECURSOR_MZ", "selected ion m/z",
                    "PEPMASS", "pepmass", "ExactMass", "MW"):
            if key in self.meta:
                val = self.meta[key]
                # PEPMASS may be "mz intensity" — keep just the m/z.
                return val.split()[0] if isinstance(val, str) and val else val
        return ""

    def field(self, key):
        return self.meta.get(key, "")

    def search_blob(self):
        """All searchable text, lower-cased, cached lazily."""
        if self._blob is None:
            self._blob = " \n ".join(
                f"{k}: {v}" for k, v in self.meta.items()).lower()
        return self._blob


def _split_peak_line(line):
    """Yield (mz, intensity) pairs from one peak line.

    Handles the common variants:
        "123.456 789"
        "123.456\t789"
        "123.456:789"
        "123.456 789; 200.0 50; ..."   (multiple pairs per line)
    A trailing annotation token (peak label) after the intensity is ignored.
    """
    # Replace common inner separators with spaces, then read numbers greedily.
    tokens = (line.replace(";", " ")
                  .replace(",", " ")
                  .replace(":", " ")
                  .replace("\t", " ")
                  .split())
    nums = []
    for tok in tokens:
        try:
            nums.append(float(tok))
        except ValueError:
            # non-numeric token (annotation) - stop reading this pair's extras
            # but keep going in case it's noise between numbers
            continue
    # pair them up
    for i in range(0, len(nums) - 1, 2):
        yield (nums[i], nums[i + 1])


def parse_msp(path, progress=None, cancel=None, sink=None):
    """Parse an .msp file into a list of Spectrum objects.

    The format is a sequence of records. Each record is a block of
    "Key: Value" metadata lines followed by peak lines. Records are
    delimited by a new "Name:" field and/or blank lines. If `sink` is given
    each Spectrum is delivered to it as parsed (streaming).
    """
    spectra = []
    src = os.path.basename(path)
    count = [0]
    meta = {}
    peaks = []
    expecting_peaks = False
    num_peaks_declared = None

    def flush():
        nonlocal meta, peaks, expecting_peaks, num_peaks_declared
        if meta or peaks:
            sp = Spectrum(dict(meta), list(peaks), source=src)
            if sink is not None:
                sink(sp)
            else:
                spectra.append(sp)
            count[0] += 1
            if progress is not None and count[0] % 2000 == 0:
                progress(count[0])
        meta = {}
        peaks = []
        expecting_peaks = False
        num_peaks_declared = None

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n").rstrip("\r")
            stripped = line.strip()

            if not stripped:
                # blank line: usually ends a record (only if we have content)
                if meta or peaks:
                    flush()
                continue

            # A new "Name:" starts a new record.
            low = stripped.lower()
            if low.startswith("name:") and (meta or peaks):
                flush()

            # metadata line?  "Key: value" — but not while reading peaks,
            # unless it clearly looks like "word: ..." with a non-numeric key.
            colon = stripped.find(":")
            is_meta = False
            if colon > 0 and not expecting_peaks:
                key = stripped[:colon].strip()
                # a metadata key shouldn't start with a digit (that's a peak)
                if key and not key[0].isdigit():
                    is_meta = True

            if is_meta:
                key = stripped[:colon].strip()
                value = stripped[colon + 1:].strip()
                meta[key] = value
                if key.lower().replace(" ", "") in ("numpeaks", "num_peaks"):
                    expecting_peaks = True
                    try:
                        num_peaks_declared = int(value)
                    except ValueError:
                        num_peaks_declared = None
            else:
                # peak line
                expecting_peaks = True
                for mz, inten in _split_peak_line(stripped):
                    peaks.append((mz, inten))

    flush()
    if progress is not None:
        progress(count[0])
    return spectra


def _open_text(path):
    """Open a (optionally gzip-compressed) text file for reading."""
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_mgf(path, progress=None, cancel=None, sink=None):
    """Parse a Mascot Generic Format (.mgf) file into Spectrum objects.

    Format: one or more BEGIN IONS / END IONS blocks. Inside a block,
    "KEY=value" lines are metadata and "mz intensity [charge]" lines are peaks.
    If `sink` is given each Spectrum is delivered to it as parsed (streaming).
    """
    spectra = []
    src = os.path.basename(path)
    count = 0
    meta = {}
    peaks = []
    in_block = False

    with _open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            upper = line.upper()
            if upper == "BEGIN IONS":
                meta, peaks, in_block = {}, [], True
                continue
            if upper == "END IONS":
                if meta or peaks:
                    sp = Spectrum(meta, peaks, source=src)
                    if sink is not None:
                        sink(sp)
                    else:
                        spectra.append(sp)
                    count += 1
                    if progress is not None and count % 2000 == 0:
                        progress(count)
                        if cancel is not None and cancel.is_set():
                            break
                meta, peaks, in_block = {}, [], False
                continue
            if not in_block:
                continue
            # metadata "KEY=value" (key starts with a letter)
            if "=" in line and line[0].isalpha():
                key, val = line.split("=", 1)
                meta[key.strip()] = val.strip()
                continue
            # otherwise a peak line: mz intensity [charge]
            parts = line.replace("\t", " ").split()
            try:
                mz = float(parts[0])
                inten = float(parts[1]) if len(parts) > 1 else 0.0
                peaks.append((mz, inten))
            except (ValueError, IndexError):
                pass

    if progress is not None:
        progress(count)
    return spectra


# cvParam accessions that describe the encoding of a binary array (not data
# we want to surface as metadata).
_BINARY_CV = {
    "MS:1000523", "MS:1000521", "MS:1000522", "MS:1000519",  # precisions
    "MS:1000574", "MS:1000576",                              # compression
    "MS:1000514", "MS:1000515",                              # array types
}


def _local(tag):
    """Strip an XML namespace from a tag name."""
    return tag.rsplit("}", 1)[-1]


def _decode_binary(text, is_64, zlib_comp):
    """Decode an mzML base64 (optionally zlib) binary array to array('d').

    Uses array.frombytes (C speed) instead of struct.unpack + list, which is
    much faster and far lighter for large arrays.
    """
    if not text:
        return array("d")
    data = base64.b64decode(text)
    if zlib_comp:
        data = zlib.decompress(data)
    width = 8 if is_64 else 4
    count = len(data) // width
    if count == 0:
        return array("d")
    arr = array("d" if is_64 else "f")
    arr.frombytes(data[:count * width])
    if sys.byteorder == "big":          # mzML arrays are little-endian
        arr.byteswap()
    return arr if is_64 else array("d", arr)


def _mzml_ms_level(elem):
    """Cheaply read the ms level, stopping before the (huge) binary arrays.

    The 'ms level' cvParam always precedes binaryDataArrayList in mzML, so we
    can decide whether to skip a spectrum without decoding its peak data.
    """
    for node in elem.iter():
        lt = _local(node.tag)
        if lt == "cvParam" and node.get("accession") == "MS:1000511":
            try:
                return int(node.get("value"))
            except (ValueError, TypeError):
                return None
        if lt == "binaryDataArray":
            break
    return None


def _parse_mzml_spectrum(elem, src):
    """Build a Spectrum from one <spectrum> element (MS2+ only).

    MS1 spectra are skipped *before* decoding their binary arrays — on
    profile/large MS1 scans (e.g. timsTOF) that decode is the dominant cost.
    """
    ms_level = _mzml_ms_level(elem)
    if ms_level is not None and ms_level < 2:
        return None

    meta = {}
    mz_arr = None
    int_arr = None

    for node in elem.iter():
        lt = _local(node.tag)
        if lt == "cvParam":
            acc = node.get("accession", "")
            name = node.get("name", "")
            value = node.get("value", "")
            if acc in _BINARY_CV or not name:
                continue
            # keep the first occurrence of a human-readable parameter.
            # Intern the key (few distinct names) to share strings across the
            # (possibly hundreds of thousands of) spectra.
            if name not in meta:
                meta[sys.intern(name)] = value
        elif lt == "binaryDataArray":
            is_64 = False
            zcomp = False
            kind = None
            bin_text = None
            for child in node.iter():
                clt = _local(child.tag)
                if clt == "cvParam":
                    acc = child.get("accession", "")
                    if acc == "MS:1000523":
                        is_64 = True
                    elif acc == "MS:1000521":
                        is_64 = False
                    elif acc == "MS:1000574":
                        zcomp = True
                    elif acc == "MS:1000514":
                        kind = "mz"
                    elif acc == "MS:1000515":
                        kind = "intensity"
                elif clt == "binary":
                    bin_text = child.text
            if kind == "mz":
                mz_arr = _decode_binary(bin_text, is_64, zcomp)
            elif kind == "intensity":
                int_arr = _decode_binary(bin_text, is_64, zcomp)

    meta[sys.intern("id")] = elem.get("id", "")
    if ms_level is not None:
        meta[sys.intern("ms level")] = str(ms_level)

    if mz_arr is None:
        mz_arr = array("d")
    if int_arr is None:
        int_arr = array("d")
    return Spectrum(meta, source=src, mz=mz_arr, inten=int_arr)


def parse_mzml(path, progress=None, cancel=None, sink=None):
    """Parse an mzML (.mzML / .mzML.gz) file into Spectrum objects (MS2+).

    MS1 spectra are detected and skipped without decoding their (large) peak
    arrays. If `sink` is given, each MS2 Spectrum is delivered to it as soon as
    it is parsed (streaming) and the return value is empty; otherwise all are
    collected and returned.
    """
    spectra = []
    src = os.path.basename(path)
    if path.lower().endswith(".gz"):
        fh = gzip.open(path, "rb")
    else:
        fh = open(path, "rb")
    scanned = 0
    try:
        for _event, elem in ET.iterparse(fh, events=("end",)):
            if _local(elem.tag) == "spectrum":
                spec = _parse_mzml_spectrum(elem, src)
                if spec is not None:
                    if sink is not None:
                        sink(spec)
                    else:
                        spectra.append(spec)
                elem.clear()
                scanned += 1
                if scanned % 2000 == 0:
                    if progress is not None:
                        progress(len(spectra))
                    if cancel is not None and cancel.is_set():
                        break
    finally:
        fh.close()
    if progress is not None:
        progress(len(spectra))
    return spectra


def parse_file(path, progress=None, cancel=None, sink=None):
    """Dispatch to the right parser based on file extension."""
    kw = dict(progress=progress, cancel=cancel, sink=sink)
    name = path.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    if name.endswith(".mgf"):
        return parse_mgf(path, **kw)
    if name.endswith(".mzml"):
        return parse_mzml(path, **kw)
    if name.endswith(".msp"):
        return parse_msp(path, **kw)
    # Unknown extension: sniff the first non-blank line.
    try:
        with _open_text(path) as fh:
            head = ""
            for line in fh:
                if line.strip():
                    head = line.strip()
                    break
    except Exception:  # noqa: BLE001
        head = ""
    if head.upper().startswith("BEGIN IONS"):
        return parse_mgf(path, **kw)
    if head.startswith("<?xml") or "<mzML" in head or "<indexedmzML" in head:
        return parse_mzml(path, **kw)
    return parse_msp(path, **kw)


# --------------------------------------------------------------------------- #
# Spectrum plot canvas
# --------------------------------------------------------------------------- #

class SpectrumPlot(tk.Canvas):
    """A lightweight stick-plot for a mass spectrum, drawn on a Canvas."""

    PAD_L = 70   # left margin (y axis labels)
    PAD_R = 20
    PAD_T = 30   # top margin (title)
    PAD_B = 50   # bottom margin (x axis labels)

    def __init__(self, master, **kw):
        super().__init__(master, background="white",
                         highlightthickness=0, **kw)
        self.peaks = []          # [(mz, inten)]
        self.title = ""
        self.relative = True     # show intensity as % of base peak
        self.scheme = DEFAULT_SCHEME   # peak colour scheme
        self.bar_width = 2       # stick thickness (px)
        self.transparency = 0    # 0 (opaque) .. 100 (invisible)
        self.theme = BG_THEMES[DEFAULT_BG]   # background colour theme
        self.precursor = None    # precursor m/z of current spectrum
        self.limit_prec = False  # show only peaks <= precursor + 5
        self.on_context = None   # callback(event) for right-click menu

        # current spectrum's metadata (for annotation value lookup)
        self.meta = {}
        self._name = ""
        self._source = ""
        self._num = 0

        # metadata annotations on the plot (persist across spectra; positions
        # are fractions of the canvas so they ignore zoom and window resize).
        # Each: {field, fx, fy, size, show_field}.  "Name" is a special field.
        self.annotations = [dict(field="Name", fx=0.5, fy=0.06, size=12,
                                 show_field=False, bold=True, italic=False,
                                 color="")]
        self._annot_items = {}   # canvas item id -> annotation index
        self._drag_annot = None
        self._annot_off = (0, 0)
        self._pan = None         # 'x' or 'y' while panning an axis
        self._pan_start = None
        self._base = 1.0         # fixed intensity reference (no y rescale)
        self._xtitle_item = None  # canvas ids of axis-title text (for dblclick)
        self._ytitle_item = None

        # axis / label styling (configured from a separate menu)
        self.show_xticks = True
        self.show_yticks = True
        self.show_titles = True        # axis titles (m/z, intensity)
        self.show_peaklabels = True    # per-peak m/z labels
        self.peaklabel_auto = True     # auto count by spacing (else fixed n)
        self.peaklabel_n = 10
        self.axis_fontsize = 8         # tick-label font
        self.peaklabel_fontsize = 8
        self.xtitle = ""               # custom x-axis title ("" = auto m/z)
        self.ytitle = ""               # custom y-axis title ("" = auto)
        # X and Y axis titles are styled independently
        self.xtitle_fontsize = 9
        self.ytitle_fontsize = 9
        self.xtitle_bold = False
        self.xtitle_italic = False
        self.ytitle_bold = False
        self.ytitle_italic = False
        self.tick_bold = False
        self.tick_italic = False
        self.peaklabel_bold = False
        self.peaklabel_italic = False
        # per-element text colour overrides ("" = use the background theme)
        self.tick_color = ""
        self.xtitle_color = ""
        self.ytitle_color = ""
        self.peaklabel_color = ""
        self._dlg = None               # single open plot dialog
        self.view_min = None     # current x (m/z) view window
        self.view_max = None
        self._full_min = None
        self._full_max = None
        self.y_lo = 0.0          # current y view window, as a fraction of base
        self.y_hi = 1.0
        self._drag_start = None
        self._drag_rect = None
        self._tip = None

        self.bind("<Configure>", lambda e: self.redraw())
        self.bind("<Motion>", self._on_motion)
        self.bind("<Leave>", lambda e: self._hide_tip())
        self.bind("<ButtonPress-1>", self._on_press)
        self.bind("<B1-Motion>", self._on_drag)
        self.bind("<ButtonRelease-1>", self._on_release)
        self.bind("<Double-Button-1>", self._on_double)
        self.bind("<Button-3>", self._on_right_click)
        self.bind("<MouseWheel>", self._on_wheel)          # Windows / macOS
        self.bind("<Button-4>",
                  lambda e: self._wheel_zoom(0.8, e.x, e.y, True))  # Linux
        self.bind("<Button-5>",
                  lambda e: self._wheel_zoom(1.25, e.x, e.y, False))

    # -- public API ---------------------------------------------------------
    def set_spectrum(self, spec):
        self.peaks = list(spec.peaks) if spec else []
        self.title = spec.name if spec else ""
        self.precursor = _prec_float(spec) if spec else None
        self.meta = dict(spec.meta) if spec else {}
        self._name = spec.name if spec else ""
        self._source = spec.source if spec else ""
        self._num = spec.num if spec else 0
        self._recompute_range()
        self.y_lo, self.y_hi = 0.0, 1.0
        self.redraw()

    def annotation_value(self, ann):
        """Text value of an annotation for the current spectrum ('' if none)."""
        field = ann["field"]
        if field == "Name":
            v = self._name
            return "" if v in ("", "(unnamed)") else v
        return str(self.meta.get(field, ""))

    def add_annotation(self, field, show_field=False):
        n = len(self.annotations)
        self.annotations.append(dict(
            field=field, fx=0.5, fy=min(0.06 + 0.06 * n, 0.9),
            size=12, show_field=show_field, bold=True, italic=False,
            color=""))
        self.redraw()

    def clear_annotations(self):
        self.annotations = []
        self.redraw()

    def _recompute_range(self):
        """Set the full x range + fixed intensity base from the peaks.

        The base (100 %% reference) is taken over the whole visible range so
        that zooming in x does NOT rescale the y axis.
        """
        if self.peaks:
            mzs = [p[0] for p in self.peaks]
            lo, hi = min(mzs), max(mzs)
            if self.limit_prec and self.precursor is not None:
                hi = min(hi, self.precursor + 5.0)
            span = max(hi - lo, 1.0)
            self._full_min = lo - span * 0.03
            self._full_max = hi + span * 0.03
            in_range = [it for mz, it in self.peaks
                        if self._full_min <= mz <= self._full_max]
            self._base = max(in_range, default=1.0) or 1.0
        else:
            self._full_min, self._full_max = 0.0, 100.0
            self._base = 1.0
        self.view_min, self.view_max = self._full_min, self._full_max

    def set_relative(self, flag):
        self.relative = bool(flag)
        self.redraw()

    def set_scheme(self, scheme):
        self.scheme = scheme
        self.redraw()

    def set_bg_theme(self, name):
        self.theme = BG_THEMES.get(name, BG_THEMES[DEFAULT_BG])
        self.configure(background=_hex(self.theme["bg"]))
        self.redraw()

    def set_limit_prec(self, flag):
        self.limit_prec = bool(flag)
        self._recompute_range()
        self.redraw()

    def set_axis_style(self, **kw):
        for k, v in kw.items():
            if hasattr(self, k):
                setattr(self, k, v)
        self.redraw()

    def set_bar_style(self, width, transparency):
        self.bar_width = max(1, int(width))
        self.transparency = max(0, min(100, int(transparency)))
        self.redraw()

    def reset_x(self):
        self.view_min, self.view_max = self._full_min, self._full_max
        self.redraw()

    def reset_y(self):
        self.y_lo, self.y_hi = 0.0, 1.0
        self.redraw()

    def reset_view(self):
        self.view_min, self.view_max = self._full_min, self._full_max
        self.y_lo, self.y_hi = 0.0, 1.0
        self.redraw()

    def _on_right_click(self, event):
        idx = self._annot_at(event.x, event.y)
        if idx is not None:
            self._annotation_menu(event, idx)
            return
        if self.on_context is not None:
            self.on_context(event)
        else:
            self.reset_view()

    def _annotation_menu(self, event, idx):
        a = self.annotations[idx]
        m = tk.Menu(self, tearoff=0)
        label = ("Show value only" if a.get("show_field")
                 else "Show \"Field: value\"")

        def toggle():
            a["show_field"] = not a.get("show_field")
            self.redraw()

        def bigger():
            a["size"] = min(60, int(a["size"]) + 2)
            self.redraw()

        def smaller():
            a["size"] = max(6, int(a["size"]) - 2)
            self.redraw()

        def delete():
            del self.annotations[idx]
            self.redraw()

        m.add_command(label=f"Label: {a['field']}", state="disabled")
        m.add_separator()
        m.add_command(label=label, command=toggle)
        m.add_command(label="Set size…  (or double-click)",
                      command=lambda: self._annotation_size_dialog(idx))
        m.add_command(label="Bigger", command=bigger)
        m.add_command(label="Smaller", command=smaller)
        m.add_separator()
        m.add_command(label="Delete label", command=delete)
        try:
            m.tk_popup(event.x_root, event.y_root)
        finally:
            m.grab_release()

    # -- coordinate helpers -------------------------------------------------
    def _plot_area(self):
        w = self.winfo_width()
        h = self.winfo_height()
        x0 = self.PAD_L
        x1 = w - self.PAD_R
        y0 = self.PAD_T
        y1 = h - self.PAD_B
        return x0, y0, x1, y1

    def _mz_to_x(self, mz, x0, x1):
        lo, hi = self.view_min, self.view_max
        if hi <= lo:
            return x0
        return x0 + (mz - lo) / (hi - lo) * (x1 - x0)

    def _x_to_mz(self, x, x0, x1):
        lo, hi = self.view_min, self.view_max
        if x1 <= x0:
            return lo
        return lo + (x - x0) / (x1 - x0) * (hi - lo)

    def _frac_to_y(self, frac, y0, y1):
        span = self.y_hi - self.y_lo
        if span <= 0:
            return y1
        return y1 - (frac - self.y_lo) / span * (y1 - y0)

    def _y_to_frac(self, py, y0, y1):
        if y1 <= y0:
            return self.y_lo
        py = min(max(py, y0), y1)
        return self.y_lo + (y1 - py) / (y1 - y0) * (self.y_hi - self.y_lo)

    # -- drawing ------------------------------------------------------------
    def redraw(self):
        self.delete("all")
        self._xtitle_item = self._ytitle_item = None
        w = self.winfo_width()
        h = self.winfo_height()
        if w < 120 or h < 90:
            return
        x0, y0, x1, y1 = self._plot_area()
        axis_c = _hex(self.theme["axis"])

        # axes
        self.create_line(x0, y1, x1, y1, fill=axis_c)   # x axis
        self.create_line(x0, y0, x0, y1, fill=axis_c)   # y axis

        if not self.peaks:
            self.create_text((x0 + x1) / 2, (y0 + y1) / 2,
                             text="No peaks to display",
                             fill=_hex(self.theme["text"]), font=("Arial", 11))
            return

        # peaks within current view; y uses a FIXED base (no rescale on zoom)
        vis = [(mz, it) for (mz, it) in self.peaks
               if self.view_min <= mz <= self.view_max]
        base = self._base if self._base > 0 else 1.0

        # y grid + labels
        self._draw_y_axis(x0, y0, x1, y1, base)
        # x ticks
        self._draw_x_axis(x0, y0, x1, y1)

        # sticks (coloured by scheme; thickness + transparency configurable;
        # clipped to the current y view window)
        for mz, it in vis:
            frac = it / base
            if frac <= self.y_lo:          # top is below the visible band
                continue
            x = self._mz_to_x(mz, x0, x1)
            y = self._frac_to_y(min(frac, self.y_hi), y0, y1)
            color = _blend_white(scheme_color(self.scheme, frac),
                                 self.transparency)
            self.create_line(x, y1, x, y, fill=_hex(color),
                             width=self.bar_width)

        # peak m/z labels: auto (by spacing) or a fixed count
        if self.show_peaklabels:
            cand = sorted((p for p in vis if p[1] / base > self.y_lo),
                          key=lambda p: p[1], reverse=True)
            placed = []
            for mz, it in cand:
                x = self._mz_to_x(mz, x0, x1)
                if self.peaklabel_auto:
                    if any(abs(x - px) < 42 for px in placed):
                        continue
                elif len(placed) >= max(0, self.peaklabel_n):
                    break
                placed.append(x)
                frac = it / base
                y = self._frac_to_y(min(frac, self.y_hi), y0, y1)
                self.create_text(
                    x, y - 6, text=f"{mz:.4f}".rstrip("0").rstrip("."),
                    anchor="s",
                    fill=self.peaklabel_color or _hex(self.theme["label"]),
                    font=("Arial", self.peaklabel_fontsize,
                          _tk_fontstyle(self.peaklabel_bold,
                                        self.peaklabel_italic)))

        # metadata annotations (draggable, fixed under zoom)
        self._draw_annotations(w, h)

    def _draw_annotations(self, w, h):
        self._annot_items = {}
        default = _hex(self.theme["title"])
        for idx, a in enumerate(self.annotations):
            val = self.annotation_value(a)
            text = f"{a['field']}: {val}" if a.get("show_field") else val
            if not text.strip():
                continue
            fam = UI_FONT if _has_cjk(text) else "Arial"
            style = _tk_fontstyle(a.get("bold", True), a.get("italic", False))
            item = self.create_text(
                a["fx"] * w, a["fy"] * h, text=text, anchor="center",
                fill=a.get("color") or default,
                font=(fam, int(a["size"]), style), tags=("annot",))
            self._annot_items[item] = idx

    def _annot_at(self, x, y):
        for item in reversed(self.find_overlapping(x - 2, y - 2, x + 2, y + 2)):
            if item in self._annot_items:
                return self._annot_items[item]
        return None

    def _draw_y_axis(self, x0, y0, x1, y1, base):
        axis_c = _hex(self.theme["axis"])
        text_c = _hex(self.theme["text"])
        tick_c = self.tick_color or text_c
        title_c = self.ytitle_color or text_c
        tick_style = _tk_fontstyle(self.tick_bold, self.tick_italic)
        title_style = _tk_fontstyle(self.ytitle_bold, self.ytitle_italic)
        if self.show_yticks:
            for i in range(0, 6):
                frac = self.y_lo + (i / 5.0) * (self.y_hi - self.y_lo)
                y = self._frac_to_y(frac, y0, y1)
                self.create_line(x0 - 4, y, x0, y, fill=axis_c)
                if self.relative:
                    label = f"{frac * 100:.0f}%"
                else:
                    label = _fmt_si(frac * base)
                self.create_text(x0 - 8, y, anchor="e", text=label,
                                 fill=tick_c,
                                 font=("Arial", self.axis_fontsize,
                                       tick_style))
        if self.show_titles:
            ytitle = self.ytitle or ("Relative intensity" if self.relative
                                     else "Intensity")
            self._ytitle_item = self.create_text(
                16, (y0 + y1) / 2, angle=90, text=ytitle, fill=title_c,
                font=("Arial", self.ytitle_fontsize, title_style),
                tags=("ytitle",))

    def _draw_x_axis(self, x0, y0, x1, y1):
        axis_c = _hex(self.theme["axis"])
        text_c = _hex(self.theme["text"])
        tick_c = self.tick_color or text_c
        title_c = self.xtitle_color or text_c
        tick_style = _tk_fontstyle(self.tick_bold, self.tick_italic)
        title_style = _tk_fontstyle(self.xtitle_bold, self.xtitle_italic)
        lo, hi = self.view_min, self.view_max
        span = hi - lo
        if span <= 0:
            return
        if self.show_xticks:
            step = _nice_step(span / 8.0)
            v = (int(lo / step) + 1) * step
            while v < hi:
                x = self._mz_to_x(v, x0, x1)
                self.create_line(x, y1, x, y1 + 4, fill=axis_c)
                self.create_text(x, y1 + 7, anchor="n", text=f"{v:g}",
                                 fill=tick_c,
                                 font=("Arial", self.axis_fontsize,
                                       tick_style))
                v += step
        if self.show_titles:
            self._xtitle_item = self.create_text(
                (x0 + x1) / 2, y1 + 30, text=self.xtitle or "m/z",
                fill=title_c, font=("Arial", self.xtitle_fontsize, title_style),
                tags=("xtitle",))

    # -- interaction --------------------------------------------------------
    def _nearest_peak(self, px, py):
        x0, y0, x1, y1 = self._plot_area()
        if not (x0 <= px <= x1):
            return None
        best = None
        best_dx = 1e9
        for mz, it in self.peaks:
            if not (self.view_min <= mz <= self.view_max):
                continue
            x = self._mz_to_x(mz, x0, x1)
            dx = abs(x - px)
            if dx < best_dx:
                best_dx = dx
                best = (mz, it)
        if best and best_dx <= 12:
            return best
        return None

    def _on_motion(self, e):
        peak = self._nearest_peak(e.x, e.y)
        if peak:
            mz, it = peak
            self._show_tip(e.x, e.y, f"m/z {mz:.4f}\nint {it:g}")
        else:
            self._hide_tip()

    def _show_tip(self, x, y, text):
        self._hide_tip()
        tx, ty = x + 12, y + 8
        self._tip = self.create_text(
            tx, ty, anchor="nw", text=text, fill="#000",
            font=("Consolas", 9), tags="tip")
        bbox = self.bbox(self._tip)
        if bbox:
            rect = self.create_rectangle(
                bbox[0] - 3, bbox[1] - 2, bbox[2] + 3, bbox[3] + 2,
                fill="#ffffe0", outline="#aaa", tags="tip")
            self.tag_lower(rect, self._tip)

    def _hide_tip(self):
        self.delete("tip")
        self._tip = None

    def _on_press(self, e):
        idx = self._annot_at(e.x, e.y)
        if idx is not None:                    # start moving an annotation
            self._drag_annot = idx
            w, h = self.winfo_width(), self.winfo_height()
            a = self.annotations[idx]
            self._annot_off = (e.x - a["fx"] * w, e.y - a["fy"] * h)
            return
        x0, y0, x1, y1 = self._plot_area()
        # dragging on an axis pans it (useful after zoom)
        if e.x <= x0 and y0 <= e.y <= y1:
            self._pan = "y"
            self._pan_start = (e.x, e.y, self.y_lo, self.y_hi)
            return
        if e.y >= y1 and x0 <= e.x <= x1:
            self._pan = "x"
            self._pan_start = (e.x, e.y, self.view_min, self.view_max)
            return
        self._drag_start = (e.x, e.y)

    def _on_drag(self, e):
        if self._drag_annot is not None:
            w = max(self.winfo_width(), 1)
            h = max(self.winfo_height(), 1)
            ox, oy = self._annot_off
            a = self.annotations[self._drag_annot]
            a["fx"] = min(max((e.x - ox) / w, 0.0), 1.0)
            a["fy"] = min(max((e.y - oy) / h, 0.0), 1.0)
            self.redraw()
            return
        if self._pan is not None:
            self._do_pan(e)
            return
        if self._drag_start is None:
            return
        if self._drag_rect:
            self.delete(self._drag_rect)
        x0, y0, x1, y1 = self._plot_area()
        ax, ay = self._drag_start
        self._drag_rect = self.create_rectangle(
            ax, ay, e.x, e.y,
            outline="#3498db", fill="#3498db", stipple="gray25")

    def _do_pan(self, e):
        x0, y0, x1, y1 = self._plot_area()
        sx, sy, a, b = self._pan_start
        if self._pan == "x":
            span = b - a
            dmz = -(e.x - sx) / max(x1 - x0, 1) * span
            lo, hi = a + dmz, b + dmz
            if lo < self._full_min:
                lo, hi = self._full_min, self._full_min + span
            if hi > self._full_max:
                lo, hi = self._full_max - span, self._full_max
            self.view_min, self.view_max = lo, hi
        else:
            span = b - a
            dfr = (e.y - sy) / max(y1 - y0, 1) * span
            lo, hi = a + dfr, b + dfr
            if lo < 0.0:
                lo, hi = 0.0, span
            if hi > 1.0:
                lo, hi = 1.0 - span, 1.0
            self.y_lo, self.y_hi = lo, hi
        self.redraw()

    def _on_release(self, e):
        if self._drag_annot is not None:
            self._drag_annot = None
            return
        if self._pan is not None:
            self._pan = None
            return
        if self._drag_rect:
            self.delete(self._drag_rect)
            self._drag_rect = None
        if self._drag_start is None:
            return
        x0, y0, x1, y1 = self._plot_area()
        ax, ay = self._drag_start
        bx, by = e.x, e.y
        self._drag_start = None
        # A drag wide enough zooms x; tall enough zooms y; a box does both.
        # (horizontal drag → x only, vertical drag → y only, box → both)
        if abs(bx - ax) < 6 and abs(by - ay) < 6:
            return                       # treat as a click
        if abs(bx - ax) >= 6:
            lo = self._x_to_mz(min(ax, bx), x0, x1)
            hi = self._x_to_mz(max(ax, bx), x0, x1)
            if hi - lo > 1e-6:
                self.view_min, self.view_max = lo, hi
        if abs(by - ay) >= 6:
            f1 = self._y_to_frac(ay, y0, y1)
            f2 = self._y_to_frac(by, y0, y1)
            lo, hi = max(0.0, min(f1, f2)), min(1.0, max(f1, f2))
            if hi - lo > 0.005:
                self.y_lo, self.y_hi = lo, hi
        self.redraw()

    def _title_at(self, x, y):
        items = self.find_overlapping(x - 2, y - 2, x + 2, y + 2)
        if self._xtitle_item in items:
            return "x"
        if self._ytitle_item in items:
            return "y"
        return None

    def _on_double(self, e):
        idx = self._annot_at(e.x, e.y)
        if idx is not None:                    # double-click a label → size box
            self._annotation_size_dialog(idx)
            return
        which = self._title_at(e.x, e.y)
        if which:                              # double-click axis title → edit
            self._axis_title_dialog(which)
            return
        x0, y0, x1, y1 = self._plot_area()
        px, py = e.x, e.y
        if x0 < px < x1 and y0 < py < y1:      # inside plot → reset both
            self.reset_view()
        elif px <= x0 and y0 - 6 <= py <= y1 + 6:   # on the y axis
            self.reset_y()
        elif py >= y1 and x0 - 6 <= px <= x1 + 6:   # on the x axis
            self.reset_x()
        else:
            self.reset_view()

    def _open_dialog(self, title):
        """Create a single, centered, modal dialog (freezes other windows)."""
        if self._dlg is not None and self._dlg.winfo_exists():
            self._dlg.lift()
            return None
        dlg = tk.Toplevel(self)
        dlg.title(title)
        dlg.transient(self.winfo_toplevel())
        dlg.resizable(False, False)
        self._dlg = dlg
        return dlg

    def _finish_dialog(self, dlg):
        """Centre the dialog on the main window and grab focus (modal)."""
        dlg.update_idletasks()
        root = self.winfo_toplevel()
        x = root.winfo_rootx() + (root.winfo_width() - dlg.winfo_width()) // 2
        y = root.winfo_rooty() + (root.winfo_height() - dlg.winfo_height()) // 2
        dlg.geometry(f"+{max(x, 0)}+{max(y, 0)}")
        try:
            dlg.grab_set()
        except Exception:  # noqa: BLE001
            pass

    def _annotation_size_dialog(self, idx):
        a = self.annotations[idx]
        dlg = self._open_dialog("Label size")
        if dlg is None:
            return
        size = tk.IntVar(value=int(a["size"]))
        bold = tk.BooleanVar(value=a.get("bold", True))
        italic = tk.BooleanVar(value=a.get("italic", False))
        color = tk.StringVar(value=a.get("color", ""))

        def apply(*_a):
            try:
                v = int(float(size.get()))
            except (ValueError, tk.TclError):
                return
            a["size"] = max(6, min(72, v))
            a["bold"] = bold.get()
            a["italic"] = italic.get()
            a["color"] = color.get()
            self.redraw()

        top = ttk.Frame(dlg)
        top.pack(padx=14, pady=(14, 4))
        ttk.Label(top, text=f"“{a['field']}”  font size (px):").pack(
            side="left")
        sp = tk.Spinbox(top, from_=6, to=72, width=5, textvariable=size,
                        command=apply)
        sp.pack(side="left", padx=(6, 0))
        sp.bind("<KeyRelease>", apply)
        tk.Scale(dlg, from_=6, to=72, orient="horizontal", length=240,
                 variable=size, command=apply).pack(padx=14)
        sty = ttk.Frame(dlg)
        sty.pack(padx=14, pady=(2, 0), anchor="w")
        ttk.Checkbutton(sty, text="Bold", variable=bold,
                        command=apply).pack(side="left")
        ttk.Checkbutton(sty, text="Italic", variable=italic,
                        command=apply).pack(side="left", padx=(10, 0))
        _color_row(dlg, "Colour", color, self.theme["title"],
                   self.winfo_toplevel()._pick_color, apply).pack(
            anchor="w", padx=14, pady=(4, 0))
        ttk.Button(dlg, text="Close", command=dlg.destroy).pack(
            anchor="e", padx=14, pady=(6, 14))
        self._finish_dialog(dlg)

    def _axis_title_dialog(self, which):
        default = "m/z" if which == "x" else (
            "Relative intensity" if self.relative else "Intensity")
        cur = self.xtitle if which == "x" else self.ytitle
        dlg = self._open_dialog(f"{which.upper()}-axis title")
        if dlg is None:
            return
        is_x = which == "x"
        var = tk.StringVar(value=cur)
        fsize = tk.IntVar(value=int(self.xtitle_fontsize if is_x
                                    else self.ytitle_fontsize))
        bold = tk.BooleanVar(value=self.xtitle_bold if is_x
                             else self.ytitle_bold)
        italic = tk.BooleanVar(value=self.xtitle_italic if is_x
                               else self.ytitle_italic)
        color = tk.StringVar(value=self.xtitle_color if is_x
                             else self.ytitle_color)

        def apply(*_a):
            try:
                fs = max(5, min(40, int(fsize.get())))
            except (ValueError, tk.TclError):
                fs = None
            if is_x:
                self.xtitle = var.get()
                if fs:
                    self.xtitle_fontsize = fs
                self.xtitle_bold = bold.get()
                self.xtitle_italic = italic.get()
                self.xtitle_color = color.get()
            else:
                self.ytitle = var.get()
                if fs:
                    self.ytitle_fontsize = fs
                self.ytitle_bold = bold.get()
                self.ytitle_italic = italic.get()
                self.ytitle_color = color.get()
            self.redraw()

        ttk.Label(dlg, text=f"{which.upper()}-axis title "
                            f"(blank = “{default}”):").pack(
            anchor="w", padx=14, pady=(14, 4))
        ent = ttk.Entry(dlg, textvariable=var, width=34)
        ent.pack(padx=14)
        ent.bind("<KeyRelease>", apply)
        ent.focus_set()
        row = ttk.Frame(dlg)
        row.pack(anchor="w", padx=14, pady=(8, 2))
        ttk.Label(row, text="Font size (px):").pack(side="left")
        sp = tk.Spinbox(row, from_=5, to=40, width=5, textvariable=fsize,
                        command=apply)
        sp.pack(side="left", padx=(6, 0))
        sp.bind("<KeyRelease>", apply)
        tk.Scale(dlg, from_=5, to=40, orient="horizontal", length=240,
                 variable=fsize, command=apply).pack(padx=14)
        sty = ttk.Frame(dlg)
        sty.pack(anchor="w", padx=14, pady=(2, 0))
        ttk.Checkbutton(sty, text="Bold", variable=bold,
                        command=apply).pack(side="left")
        ttk.Checkbutton(sty, text="Italic", variable=italic,
                        command=apply).pack(side="left", padx=(10, 0))
        _color_row(dlg, "Colour", color, self.theme["text"],
                   self.winfo_toplevel()._pick_color, apply).pack(
            anchor="w", padx=14, pady=(4, 0))
        ttk.Button(dlg, text="Close", command=dlg.destroy).pack(
            anchor="e", padx=14, pady=(6, 14))
        self._finish_dialog(dlg)

    def _on_wheel(self, e):
        factor = 0.8 if e.delta > 0 else 1.25
        self._wheel_zoom(factor, e.x, e.y, up=(e.delta > 0))

    def _wheel_zoom(self, factor, px, py, up=True):
        idx = self._annot_at(px, py)
        if idx is not None:                    # resize the label under cursor
            a = self.annotations[idx]
            a["size"] = min(60, max(6, int(a["size"]) + (1 if up else -1)))
            self.redraw()
            return
        x0, y0, x1, y1 = self._plot_area()
        if px < x0:                # over the y-axis area → zoom y
            self._zoom_y(factor, py)
        else:                      # otherwise zoom x
            self._zoom_x(factor, px)

    # If the cursor is within this fraction of an axis end, pin that end.
    EDGE_ANCHOR = 0.15

    def _zoom_x(self, factor, px):
        x0, y0, x1, y1 = self._plot_area()
        t = (px - x0) / max(x1 - x0, 1)
        if t <= self.EDGE_ANCHOR:          # near left end → pin view_min
            center = self.view_min
        elif t >= 1 - self.EDGE_ANCHOR:    # near right end → pin view_max
            center = self.view_max
        else:
            center = self._x_to_mz(px, x0, x1)
        lo = center - (center - self.view_min) * factor
        hi = center + (self.view_max - center) * factor
        lo = max(lo, self._full_min)
        hi = min(hi, self._full_max)
        if hi - lo > 1e-6:
            self.view_min, self.view_max = lo, hi
            self.redraw()

    def _zoom_y(self, factor, py):
        x0, y0, x1, y1 = self._plot_area()
        t = (py - y0) / max(y1 - y0, 1)    # 0 at top (y_hi), 1 at bottom (y_lo)
        if t <= self.EDGE_ANCHOR:          # near top → pin y_hi
            center = self.y_hi
        elif t >= 1 - self.EDGE_ANCHOR:    # near bottom → pin y_lo (baseline)
            center = self.y_lo
        else:
            center = self._y_to_frac(py, y0, y1)
        lo = center - (center - self.y_lo) * factor
        hi = center + (self.y_hi - center) * factor
        lo = max(0.0, lo)
        hi = min(1.0, hi)
        if hi - lo > 0.005:
            self.y_lo, self.y_hi = lo, hi
            self.redraw()


def _nice_step(raw):
    """Return a 'nice' axis step (1/2/5 * 10^n) close to raw."""
    import math
    if raw <= 0:
        return 1.0
    exp = math.floor(math.log10(raw))
    base = 10 ** exp
    for m in (1, 2, 5, 10):
        if raw <= m * base:
            return m * base
    return 10 * base


def _fmt_si(v):
    """Compact number formatting for the intensity axis."""
    a = abs(v)
    if a >= 1e9:
        return f"{v/1e9:.1f}G"
    if a >= 1e6:
        return f"{v/1e6:.1f}M"
    if a >= 1e3:
        return f"{v/1e3:.1f}k"
    return f"{v:.0f}"


# --------------------------------------------------------------------------- #
# Chemistry: monoisotopic mass from a formula + adduct m/z
# --------------------------------------------------------------------------- #

# Monoisotopic masses (most-abundant isotope), u.
_ELEMENT_MASS = {
    "H": 1.0078250319, "D": 2.0141017779, "He": 4.002603254,
    "Li": 7.0160034366, "B": 11.0093054, "C": 12.0, "N": 14.0030740052,
    "O": 15.9949146221, "F": 18.9984031627, "Na": 22.989769282,
    "Mg": 23.985041697, "Al": 26.98153853, "Si": 27.9769265327,
    "P": 30.97376151, "S": 31.97207069, "Cl": 34.968852682,
    "K": 38.9637064864, "Ca": 39.962590863, "Fe": 55.9349363,
    "Co": 58.9331943, "Ni": 57.9353424, "Cu": 62.9295977,
    "Zn": 63.9291422, "As": 74.9215946, "Se": 79.9165218,
    "Br": 78.9183376, "Ag": 106.9050916, "I": 126.9044719,
    "Pt": 194.9647917, "Au": 196.9665688, "Hg": 201.9706434,
}

_ELECTRON_MASS = 0.00054857990907

# adduct: label -> (M multiplier, add formula, remove formula, signed charge)
# Preset adduct labels for the dropdown. Any label in the general form
# [<n>M ±A ±B …]<z><sign> also works when typed in by hand (custom adducts).
ADDUCT_PRESETS = [
    "[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M+H-H2O]+",
    "[M+2H]2+", "[M+3H]3+", "[M+2NH4]2+", "[M+H+NH4]2+", "[2M+H]+",
    "[2M+Na]+", "[2M+NH4]+", "[M]+",
    "[M-H]-", "[M+Cl]-", "[M+HCOO]-", "[M+CH3COO]-",
    "[M-2H]2-", "[2M-H]-", "[M]-",
]

_FORMULA_TOKEN = re.compile(r"[A-Z][a-z]?|\(|\)|\d+")


def parse_formula(formula):
    """Parse a chemical formula (with optional parentheses) to {element: count}."""
    text = formula.replace(" ", "")
    tokens = _FORMULA_TOKEN.findall(text)
    if "".join(tokens) != text:
        raise ValueError(f"Invalid characters in formula: {formula!r}")
    stack = [{}]
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok == "(":
            stack.append({})
            i += 1
        elif tok == ")":
            grp = stack.pop()
            i += 1
            mult = 1
            if i < len(tokens) and tokens[i].isdigit():
                mult = int(tokens[i])
                i += 1
            top = stack[-1]
            for el, c in grp.items():
                top[el] = top.get(el, 0) + c * mult
        elif tok.isdigit():
            i += 1  # stray number — ignore
        else:
            i += 1
            cnt = 1
            if i < len(tokens) and tokens[i].isdigit():
                cnt = int(tokens[i])
                i += 1
            stack[-1][tok] = stack[-1].get(tok, 0) + cnt
    if len(stack) != 1:
        raise ValueError(f"Unbalanced parentheses in formula: {formula!r}")
    return stack[0]


def formula_mass(formula):
    """Monoisotopic neutral mass of a formula (raises on unknown element)."""
    total = 0.0
    for el, cnt in parse_formula(formula).items():
        if el not in _ELEMENT_MASS:
            raise ValueError(f"Unknown element: {el}")
        total += _ELEMENT_MASS[el] * cnt
    return total


_ADDUCT_TERM = re.compile(r"([+-])(\d*)([A-Za-z][A-Za-z0-9]*)")


def parse_adduct(adduct):
    """Parse an adduct label like [2M+H]+, [M+2NH4]2+, [M+H-H2O]+.

    Returns (M-multiplier, added mass, removed mass, signed charge).
    Accepts ']' or '}' as the closing bracket and the charge either before or
    after the sign (2+ or +2).
    """
    s = adduct.strip().replace(" ", "")
    if not s.startswith("["):
        raise ValueError(f"Adduct must start with '[': {adduct!r}")
    close = max(s.find("]"), s.find("}"))
    if close < 0:
        raise ValueError(f"Adduct missing ']': {adduct!r}")
    core, chg = s[1:close], s[close + 1:]
    if "+" in chg:
        sign = 1
    elif "-" in chg:
        sign = -1
    else:
        raise ValueError(f"Adduct needs a charge sign: {adduct!r}")
    digits = "".join(c for c in chg if c.isdigit())
    z = sign * (int(digits) if digits else 1)

    cm = re.match(r"^(\d*)M(.*)$", core)
    if not cm:
        raise ValueError(f"Adduct must contain 'M': {adduct!r}")
    mult = int(cm.group(1)) if cm.group(1) else 1
    rest = cm.group(2)
    if rest and not re.fullmatch(r"([+-]\d*[A-Za-z][A-Za-z0-9]*)+", rest):
        raise ValueError(f"Cannot parse adduct: {adduct!r}")
    add = rem = 0.0
    for op, cnt, frm in _ADDUCT_TERM.findall(rest):
        mass = (int(cnt) if cnt else 1) * formula_mass(frm)
        if op == "+":
            add += mass
        else:
            rem += mass
    return mult, add, rem, z


def adduct_mz(neutral_mass, adduct):
    """m/z for a neutral monoisotopic mass under the given adduct label."""
    mult, add, rem, z = parse_adduct(adduct)
    m = neutral_mass * mult + add - rem
    # positive charge removes electrons, negative charge adds them
    m -= z * _ELECTRON_MASS
    return m / abs(z)


def formula_adduct_mz(formula, adduct):
    """Convenience: compute m/z straight from a formula + adduct label."""
    return adduct_mz(formula_mass(formula), adduct)


# --------------------------------------------------------------------------- #
# Export: write spectra to .msp / .mgf, render to a PIL image
# --------------------------------------------------------------------------- #

def _fmt_mz(v):
    return f"{v:.6f}".rstrip("0").rstrip(".")


def _to_float(s):
    """Parse a float or return None (empty/invalid)."""
    try:
        return float(str(s).strip())
    except (ValueError, TypeError):
        return None


def spectrum_to_msp(spec):
    """Serialize one Spectrum to NIST .msp text."""
    lines = []
    if "Name" not in spec.meta:
        nm = spec.name
        lines.append(f"Name: {nm if nm != '(unnamed)' else spec.source}")
    for k, v in spec.meta.items():
        if k.lower().replace(" ", "") in ("numpeaks", "num_peaks"):
            continue
        lines.append(f"{k}: {v}")
    lines.append(f"Num Peaks: {len(spec.peaks)}")
    for mz, inten in spec.peaks:
        lines.append(f"{_fmt_mz(mz)} {inten:g}")
    return "\n".join(lines) + "\n"


_MGF_SAFE_KEY = re.compile(r"^[A-Za-z0-9_]+$")


def spectrum_to_mgf(spec):
    """Serialize one Spectrum to Mascot Generic Format text."""
    lines = ["BEGIN IONS"]
    title = spec.name if spec.name != "(unnamed)" else \
        f"{spec.source} #{spec.num}"
    lines.append(f"TITLE={title}")
    prec = spec.precursor
    if prec:
        try:
            lines.append(f"PEPMASS={float(prec):.6f}")
        except (ValueError, TypeError):
            pass
    for ck in ("CHARGE", "charge", "charge state"):
        if ck in spec.meta:
            lines.append(f"CHARGE={spec.meta[ck]}")
            break
    skip = {"title", "pepmass", "charge", "charge state"}
    for k, v in spec.meta.items():
        if k.lower() in skip:
            continue
        if _MGF_SAFE_KEY.match(k) and "\n" not in str(v):
            lines.append(f"{k.upper()}={v}")
    for mz, inten in spec.peaks:
        lines.append(f"{_fmt_mz(mz)} {inten:g}")
    lines.append("END IONS")
    return "\n".join(lines) + "\n"


# -- peak colour schemes --------------------------------------------------- #
# Each scheme maps a relative intensity (0..1) to an (r, g, b) colour. Single-
# stop schemes are flat; multi-stop schemes give a "fancy" intensity gradient.
COLOR_SCHEMES = {
    "Classic blue": [(31, 119, 180)],
    "Ocean": [(140, 200, 255), (8, 60, 150)],
    "Sunset (fancy)": [(255, 214, 102), (255, 122, 0),
                       (214, 40, 57), (120, 20, 110)],
    "Viridis": [(68, 1, 84), (59, 82, 139), (33, 145, 140),
                (94, 201, 98), (253, 231, 37)],
    "Magma": [(0, 0, 4), (81, 18, 124), (183, 55, 121),
              (252, 137, 97), (252, 253, 191)],
    "Forest": [(180, 225, 150), (16, 110, 50)],
    "Mono grey": [(170, 170, 170), (40, 40, 40)],
}
DEFAULT_SCHEME = "Sunset (fancy)"

# User-defined scheme (edited via the colour dialog); starts as a copy of the
# default so selecting it before editing still produces something sensible.
CUSTOM_SCHEME = "Custom…"
COLOR_SCHEMES[CUSTOM_SCHEME] = list(COLOR_SCHEMES[DEFAULT_SCHEME])


def _positioned(stops):
    """Normalise stops to a sorted [(pos, (r,g,b)), …] list.

    Accepts either plain [(r,g,b), …] (spread evenly over 0..1) or already
    positioned [(pos, (r,g,b)), …].
    """
    if not stops:
        return [(0.0, (31, 119, 180))]
    if isinstance(stops[0][1], (tuple, list)):          # positioned
        return sorted(((float(p), tuple(c)) for p, c in stops),
                      key=lambda x: x[0])
    n = len(stops)                                       # plain rgb list
    if n == 1:
        return [(0.0, tuple(stops[0]))]
    return [(i / (n - 1), tuple(stops[i])) for i in range(n)]


def _interp_stops(stops, frac):
    """Interpolate colour stops (plain or positioned) at frac -> (r,g,b)."""
    pts = _positioned(stops)
    if len(pts) == 1:
        return pts[0][1]
    frac = 0.0 if frac < 0 else 1.0 if frac > 1 else frac
    if frac <= pts[0][0]:
        return pts[0][1]
    if frac >= pts[-1][0]:
        return pts[-1][1]
    for i in range(len(pts) - 1):
        p0, c0 = pts[i]
        p1, c1 = pts[i + 1]
        if p0 <= frac <= p1:
            t = (frac - p0) / (p1 - p0) if p1 > p0 else 0.0
            return tuple(int(round(c0[k] + (c1[k] - c0[k]) * t))
                         for k in range(3))
    return pts[-1][1]


def scheme_color(scheme, frac):
    """Interpolate a colour scheme at fractional intensity frac -> (r,g,b)."""
    stops = COLOR_SCHEMES.get(scheme) or COLOR_SCHEMES["Classic blue"]
    return _interp_stops(stops, frac)


def _hex(rgb):
    return "#%02x%02x%02x" % rgb


def _hex_to_rgb(s):
    s = s.lstrip("#")
    return (int(s[0:2], 16), int(s[2:4], 16), int(s[4:6], 16))


def _color_row(parent, label, var, default_rgb, pick, on_change):
    """A label + swatch + 'Colour…' / 'Default' row bound to a hex StringVar.

    An empty var means 'use the theme default' (default_rgb).
    """
    row = ttk.Frame(parent)
    swatch = tk.Label(row, width=3, relief="solid", borderwidth=1)

    def refresh():
        swatch.configure(background=var.get() or _hex(default_rgb))

    def choose():
        cur = var.get() or _hex(default_rgb)
        rgb = pick(_hex_to_rgb(cur))
        if rgb:
            var.set(_hex(rgb))
            refresh()
            on_change()

    def clear():
        var.set("")
        refresh()
        on_change()

    ttk.Label(row, text=label, width=12).pack(side="left")
    swatch.pack(side="left", padx=(0, 6))
    ttk.Button(row, text="Colour…", width=8, command=choose).pack(side="left")
    ttk.Button(row, text="Default", width=8, command=clear).pack(
        side="left", padx=(4, 0))
    refresh()
    return row


# -- plot background themes (bg + matching axis/text/title/label colours) --- #
BG_THEMES = {
    "White": dict(bg=(255, 255, 255), axis=(68, 68, 68), text=(85, 85, 85),
                  title=(17, 17, 17), label=(192, 60, 43)),
    "Light grey": dict(bg=(244, 245, 247), axis=(90, 90, 90),
                       text=(80, 80, 80), title=(20, 20, 20),
                       label=(192, 60, 43)),
    "Sepia": dict(bg=(247, 241, 227), axis=(120, 100, 70), text=(105, 88, 58),
                  title=(60, 45, 20), label=(170, 60, 40)),
    "Dark": dict(bg=(34, 38, 43), axis=(150, 150, 150), text=(200, 200, 200),
                 title=(236, 236, 236), label=(255, 138, 128)),
    "Black": dict(bg=(0, 0, 0), axis=(130, 130, 130), text=(190, 190, 190),
                  title=(255, 255, 255), label=(255, 107, 107)),
}
DEFAULT_BG = "White"


def _prec_float(spec):
    """Precursor m/z of a spectrum as a float, or None."""
    try:
        return float(str(spec.precursor).split()[0])
    except (ValueError, TypeError, IndexError):
        return None


# -- fonts: Arial for Latin, Microsoft YaHei for CJK ------------------------ #

def _resource_path(name):
    """Locate a bundled resource, whether running from source or a frozen exe."""
    base = getattr(sys, "_MEIPASS",
                   os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base, name)


def _find_font(*names):
    dirs = [os.path.join(os.environ.get("WINDIR", r"C:\Windows"), "Fonts"),
            r"C:\Windows\Fonts",
            "/usr/share/fonts", os.path.expanduser("~/.fonts"),
            os.path.expanduser("~/Library/Fonts")]
    for d in dirs:
        for n in names:
            p = os.path.join(d, n)
            if os.path.isfile(p):
                return p
    return None


# Arial for Latin text; Microsoft YaHei (msyh.ttc) for Chinese.
ARIAL_PATH = _find_font("arial.ttf", "Arial.ttf", "DejaVuSans.ttf")
YAHEI_PATH = _find_font("msyh.ttc", "msyh.ttf", "msyhl.ttc",
                        "SourceHanSansSC-Regular.otf",
                        "NotoSansCJK-Regular.ttc")
# Latin bold / italic variants (fall back to the regular face if absent)
ARIAL_BOLD = _find_font("arialbd.ttf", "DejaVuSans-Bold.ttf") or ARIAL_PATH
ARIAL_ITALIC = _find_font("ariali.ttf", "DejaVuSans-Oblique.ttf") or ARIAL_PATH
ARIAL_BOLDIT = _find_font("arialbi.ttf",
                          "DejaVuSans-BoldOblique.ttf") or ARIAL_PATH
YAHEI_BOLD = _find_font("msyhbd.ttc", "msyhbd.ttf") or YAHEI_PATH
UI_FONT = "Microsoft YaHei UI" if _find_font("msyh.ttc") else "Arial"


def _latin_font_path(bold, italic):
    if bold and italic:
        return ARIAL_BOLDIT
    if bold:
        return ARIAL_BOLD
    if italic:
        return ARIAL_ITALIC
    return ARIAL_PATH


def _tk_fontstyle(bold, italic):
    s = ("bold " if bold else "") + ("italic" if italic else "")
    return s.strip() or "normal"


def _has_cjk(s):
    """True if the string contains CJK / fullwidth characters."""
    for ch in s:
        o = ord(ch)
        if (0x2E80 <= o <= 0x9FFF or 0x3000 <= o <= 0x30FF
                or 0xAC00 <= o <= 0xD7AF or 0xF900 <= o <= 0xFAFF
                or 0xFF00 <= o <= 0xFFEF):
            return True
    return False


# -- drawing backends (shared layout, raster + vector output) --------------- #

class _PILBackend:
    """Raster drawing backend over PIL, top-left origin, y downwards.

    Latin text is drawn in Arial, CJK text in Microsoft YaHei.
    """

    def __init__(self, width, height):
        from PIL import Image, ImageDraw, ImageFont
        self._Image = Image
        self._ImageDraw = ImageDraw
        self._ImageFont = ImageFont
        self.w, self.h = width, height
        self.img = Image.new("RGB", (width, height), "white")
        self.dr = ImageDraw.Draw(self.img)
        self._cache = {}

    def fill_bg(self, rgb):
        self.dr.rectangle([0, 0, self.w, self.h], fill=tuple(rgb))

    def _f(self, size, cjk, bold=False, italic=False):
        key = (size, cjk, bold, italic)
        if key in self._cache:
            return self._cache[key]
        ImageFont = self._ImageFont
        fnt = None
        if cjk:
            path = YAHEI_BOLD if bold else YAHEI_PATH
        else:
            path = _latin_font_path(bold, italic)
        if path:
            try:
                if path.lower().endswith(".ttc"):
                    fnt = ImageFont.truetype(path, size, index=0)
                else:
                    fnt = ImageFont.truetype(path, size)
            except Exception:  # noqa: BLE001
                fnt = None
        if fnt is None:
            for name in ("arial.ttf", "DejaVuSans.ttf", "segoeui.ttf"):
                try:
                    fnt = ImageFont.truetype(name, size)
                    break
                except Exception:  # noqa: BLE001
                    continue
        if fnt is None:
            fnt = ImageFont.load_default()
        self._cache[key] = fnt
        return fnt

    def line(self, x0, y0, x1, y1, rgb, width=1):
        self.dr.line([(x0, y0), (x1, y1)], fill=rgb, width=int(width))

    def text(self, x, y, s, rgb, size=15, anchor="lt", rotate=0,
             bold=False, italic=False):
        fnt = self._f(size, _has_cjk(s), bold, italic)
        if rotate:
            tw = int(self.dr.textlength(s, font=fnt))
            tmp = self._Image.new("RGBA", (tw + 4, size + 8), (0, 0, 0, 0))
            self._ImageDraw.Draw(tmp).text((2, 2), s, fill=rgb, font=fnt)
            tmp = tmp.rotate(rotate, expand=True)
            # anchor 'cm' style: centre on (x, y)
            self.img.paste(tmp, (int(x - tmp.width / 2),
                                 int(y - tmp.height / 2)), tmp)
            return
        pil_anchor = {"l": "l", "c": "m", "r": "r"}[anchor[0]] + \
                     {"t": "a", "m": "m", "b": "s"}[anchor[1]]
        self.dr.text((x, y), s, fill=rgb, font=fnt, anchor=pil_anchor)

    def save(self, path):
        self.img.save(path)


# Helvetica (AFM) glyph widths /1000em for the printable ASCII range.
_HELV_W = {
    ' ': 278, '!': 278, '"': 355, '#': 556, '$': 556, '%': 889, '&': 667,
    "'": 191, '(': 333, ')': 333, '*': 389, '+': 584, ',': 278, '-': 333,
    '.': 278, '/': 278, ':': 278, ';': 278, '<': 584, '=': 584, '>': 584,
    '?': 556, '@': 1015, '[': 278, '\\': 278, ']': 278, '^': 469, '_': 556,
    '`': 333, '{': 334, '|': 260, '}': 334, '~': 584,
    'A': 667, 'B': 667, 'C': 722, 'D': 722, 'E': 667, 'F': 611, 'G': 778,
    'H': 722, 'I': 278, 'J': 500, 'K': 667, 'L': 556, 'M': 833, 'N': 722,
    'O': 778, 'P': 667, 'Q': 778, 'R': 722, 'S': 667, 'T': 611, 'U': 722,
    'V': 667, 'W': 944, 'X': 667, 'Y': 667, 'Z': 611,
    'a': 556, 'b': 556, 'c': 500, 'd': 556, 'e': 556, 'f': 278, 'g': 556,
    'h': 556, 'i': 222, 'j': 222, 'k': 500, 'l': 222, 'm': 833, 'n': 556,
    'o': 556, 'p': 556, 'q': 556, 'r': 333, 's': 500, 't': 278, 'u': 556,
    'v': 500, 'w': 722, 'x': 500, 'y': 500, 'z': 500,
}
for _d in "0123456789":
    _HELV_W[_d] = 556


class _PDFBackend:
    """Minimal single-font (Helvetica) vector PDF backend, one op list/page.

    Fallback used only when reportlab is unavailable (Latin text only).
    """

    def __init__(self, path, width, height):
        self.path = path
        self.w = width
        self.h = height
        self.pages = []
        self.ops = []
        self.pages.append(self.ops)

    def new_page(self):
        self.ops = []
        self.pages.append(self.ops)

    def fill_bg(self, rgb):
        r, g, b = (c / 255.0 for c in rgb)
        self.ops.append(
            f"{r:.3f} {g:.3f} {b:.3f} rg 0 0 {self.w} {self.h} re f")

    @staticmethod
    def _text_width(s, size):
        return sum(_HELV_W.get(ch, 556) for ch in s) / 1000.0 * size

    def line(self, x0, y0, x1, y1, rgb, width=1):
        r, g, b = (c / 255.0 for c in rgb)
        Y = self.h
        self.ops.append(
            f"{r:.3f} {g:.3f} {b:.3f} RG {width:.2f} w "
            f"{x0:.2f} {Y - y0:.2f} m {x1:.2f} {Y - y1:.2f} l S")

    @staticmethod
    def _esc(s):
        return s.replace("\\", r"\\").replace("(", r"\(").replace(")", r"\)")

    def text(self, x, y, s, rgb, size=15, anchor="lt", rotate=0,
             bold=False, italic=False):
        r, g, b = (c / 255.0 for c in rgb)
        Y = self.h
        tw = self._text_width(s, size)
        asc, cap = 0.718 * size, 0.70 * size
        # baseline offset from the anchor point (screen y, downwards +)
        vy = {"t": asc, "m": cap * 0.5, "b": 0.0}[anchor[1]]
        by_top = y + vy               # baseline, screen coords
        if rotate:
            # rotate 90° CCW, centred on (x, y)
            px, py = x, Y - y
            self.ops.append(
                f"{r:.3f} {g:.3f} {b:.3f} rg BT "
                f"0 1 -1 0 {px:.2f} {py:.2f} Tm /F1 {size:.2f} Tf "
                f"-{tw/2:.2f} -{size*0.35:.2f} Td ({self._esc(s)}) Tj ET")
            return
        hx = {"l": 0.0, "c": tw / 2.0, "r": tw}[anchor[0]]
        bx = x - hx
        self.ops.append(
            f"{r:.3f} {g:.3f} {b:.3f} rg BT /F1 {size:.2f} Tf "
            f"{bx:.2f} {Y - by_top:.2f} Td ({self._esc(s)}) Tj ET")

    def save(self):
        path = self.path
        objs = []          # object bodies (1-indexed via append order)

        def add(body):
            objs.append(body)
            return len(objs)

        font_id = add("<< /Type /Font /Subtype /Type1 "
                      "/BaseFont /Helvetica >>")
        kids_ph = None
        pages_id = add("")  # placeholder, patched after we know kids
        page_ids = []
        for ops in self.pages:
            stream = "\n".join(ops)
            content_id = add(
                f"<< /Length {len(stream.encode('latin-1', 'replace'))} >>\n"
                f"stream\n{stream}\nendstream")
            page_id = add(
                f"<< /Type /Page /Parent {pages_id} 0 R "
                f"/MediaBox [0 0 {self.w} {self.h}] "
                f"/Resources << /Font << /F1 {font_id} 0 R >> >> "
                f"/Contents {content_id} 0 R >>")
            page_ids.append(page_id)
        kids = " ".join(f"{pid} 0 R" for pid in page_ids)
        objs[pages_id - 1] = (f"<< /Type /Pages /Count {len(page_ids)} "
                              f"/Kids [{kids}] >>")
        catalog_id = add(f"<< /Type /Catalog /Pages {pages_id} 0 R >>")

        out = bytearray(b"%PDF-1.4\n")
        offsets = [0]
        for i, body in enumerate(objs, start=1):
            offsets.append(len(out))
            out += f"{i} 0 obj\n{body}\nendobj\n".encode("latin-1", "replace")
        xref_pos = len(out)
        n = len(objs) + 1
        out += f"xref\n0 {n}\n".encode()
        out += b"0000000000 65535 f \n"
        for off in offsets[1:]:
            out += f"{off:010d} 00000 n \n".encode()
        out += (f"trailer\n<< /Size {n} /Root {catalog_id} 0 R >>\n"
                f"startxref\n{xref_pos}\n%%EOF").encode()
        with open(path, "wb") as fh:
            fh.write(out)


_RL_STATE = {"ready": False, "arial": False, "yahei": False,
             "arial_bd": False, "arial_it": False, "arial_bi": False,
             "yahei_bd": False}


def _ensure_rl_fonts():
    """Register Arial + Microsoft YaHei (with bold/italic variants) once."""
    if _RL_STATE["ready"]:
        return
    _RL_STATE["ready"] = True
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont

    def reg(name, path, key, idx=None):
        if not path:
            return
        try:
            if path.lower().endswith(".ttc"):
                pdfmetrics.registerFont(TTFont(name, path,
                                               subfontIndex=idx or 0))
            else:
                pdfmetrics.registerFont(TTFont(name, path))
            _RL_STATE[key] = True
        except Exception:  # noqa: BLE001
            pass

    reg("Arial", ARIAL_PATH, "arial")
    reg("Arial-Bold", ARIAL_BOLD, "arial_bd")
    reg("Arial-Italic", ARIAL_ITALIC, "arial_it")
    reg("Arial-BoldItalic", ARIAL_BOLDIT, "arial_bi")
    reg("YaHei", YAHEI_PATH, "yahei")
    reg("YaHei-Bold", YAHEI_BOLD, "yahei_bd")


class _RLBackend:
    """Vector PDF backend via reportlab, embedding Arial + YaHei subsets."""

    def __init__(self, path, width, height):
        from reportlab.pdfgen import canvas
        _ensure_rl_fonts()
        self.c = canvas.Canvas(path, pagesize=(width, height))
        self.w, self.h = width, height

    def new_page(self):
        self.c.showPage()

    def fill_bg(self, rgb):
        self.c.setFillColorRGB(*(v / 255.0 for v in rgb))
        self.c.rect(0, 0, self.w, self.h, fill=1, stroke=0)

    def _font_for(self, s, bold=False, italic=False):
        if _has_cjk(s) and _RL_STATE["yahei"]:
            return "YaHei-Bold" if bold and _RL_STATE["yahei_bd"] else "YaHei"
        if not _RL_STATE["arial"]:
            # base-14 Helvetica variants
            v = ("Helvetica-BoldOblique" if bold and italic else
                 "Helvetica-Bold" if bold else
                 "Helvetica-Oblique" if italic else "Helvetica")
            return v
        if bold and italic and _RL_STATE["arial_bi"]:
            return "Arial-BoldItalic"
        if bold and _RL_STATE["arial_bd"]:
            return "Arial-Bold"
        if italic and _RL_STATE["arial_it"]:
            return "Arial-Italic"
        return "Arial"

    def line(self, x0, y0, x1, y1, rgb, width=1):
        self.c.setStrokeColorRGB(*(v / 255.0 for v in rgb))
        self.c.setLineWidth(width)
        Y = self.h
        self.c.line(x0, Y - y0, x1, Y - y1)

    def text(self, x, y, s, rgb, size=15, anchor="lt", rotate=0,
             bold=False, italic=False):
        self.c.setFillColorRGB(*(v / 255.0 for v in rgb))
        font = self._font_for(s, bold, italic)
        Y = self.h
        asc, cap = 0.718 * size, 0.70 * size
        vy = {"t": asc, "m": cap * 0.5, "b": 0.0}[anchor[1]]
        if rotate:
            self.c.saveState()
            self.c.translate(x, Y - y)
            self.c.rotate(90)
            self.c.setFont(font, size)
            self.c.drawCentredString(0, -size * 0.35, s)
            self.c.restoreState()
            return
        self.c.setFont(font, size)
        by = Y - (y + vy)
        if anchor[0] == "l":
            self.c.drawString(x, by, s)
        elif anchor[0] == "c":
            self.c.drawCentredString(x, by, s)
        else:
            self.c.drawRightString(x, by, s)

    def save(self):
        self.c.showPage()
        self.c.save()


def _make_pdf_backend(path, width, height):
    """Prefer the reportlab vector backend (embeds Arial/YaHei); fall back."""
    try:
        import reportlab  # noqa: F401
        return _RLBackend(path, width, height)
    except Exception:  # noqa: BLE001
        return _PDFBackend(path, width, height)


def _blend_white(rgb, transparency):
    """Blend a colour toward white to fake transparency over a white ground.

    transparency: 0 (opaque) .. 100 (invisible).
    """
    if not transparency:
        return tuple(rgb)
    f = max(0, min(100, transparency)) / 100.0
    return tuple(int(round(c + (255 - c) * f)) for c in rgb)


DEFAULT_STYLE = dict(show_xticks=True, show_yticks=True, show_titles=True,
                     show_peaklabels=True, peaklabel_auto=True,
                     peaklabel_n=10, xtitle="", ytitle="",
                     axis_fontsize=8, xtitle_fontsize=9, ytitle_fontsize=9,
                     peaklabel_fontsize=8,
                     tick_bold=False, tick_italic=False,
                     xtitle_bold=False, xtitle_italic=False,
                     ytitle_bold=False, ytitle_italic=False,
                     peaklabel_bold=False, peaklabel_italic=False,
                     tick_color="", xtitle_color="", ytitle_color="",
                     peaklabel_color="")


def _annotation_value(spec, field):
    if field == "Name":
        v = spec.name
        return "" if v in ("", "(unnamed)") else v
    return str(spec.meta.get(field, ""))


def _draw_spectrum(be, width, height, spec, relative=True,
                   scheme=DEFAULT_SCHEME, bar_width=2, transparency=0,
                   theme=None, limit_prec=False, annotations=None, style=None):
    """Draw a spectrum stick-plot onto a backend (shared by PNG and PDF)."""
    theme = theme or BG_THEMES[DEFAULT_BG]
    st = {**DEFAULT_STYLE, **(style or {})}
    axis, grey, label = theme["axis"], theme["text"], theme["label"]
    be.fill_bg(theme["bg"])
    pad_l, pad_r, pad_t, pad_b = 90, 30, 55, 70
    x0, y0 = pad_l, pad_t
    x1, y1 = width - pad_r, height - pad_b
    # scale on-screen font-size settings up to the (larger) export canvas
    sf = width / 900.0

    def fs(v, floor=6):
        return max(floor, int(round(v * sf)))

    def col(custom, default_rgb):
        return _hex_to_rgb(custom) if custom else default_rgb

    be.line(x0, y1, x1, y1, axis, 1)
    be.line(x0, y0, x0, y1, axis, 1)

    def draw_annotations():
        for a in (annotations or []):
            val = _annotation_value(spec, a["field"])
            text = f"{a['field']}: {val}" if a.get("show_field") else val
            if not text.strip():
                continue
            size = max(8, int(a["size"] * width / 900.0))
            be.text(a["fx"] * width, a["fy"] * height, text,
                    col(a.get("color", ""), theme["title"]),
                    size=size, anchor="cm",
                    bold=a.get("bold", False), italic=a.get("italic", False))

    peaks = spec.peaks or []
    if not peaks:
        be.text((x0 + x1) / 2, (y0 + y1) / 2, "No peaks", grey,
                size=15, anchor="cm")
        draw_annotations()
        return

    mzs = [p[0] for p in peaks]
    lo, hi = min(mzs), max(mzs)
    pf = _prec_float(spec)
    if limit_prec and pf is not None:      # cap at precursor + 5 m/z
        hi = min(hi, pf + 5.0)
    span = max(hi - lo, 1.0)
    vmin, vmax = lo - span * 0.03, hi + span * 0.03
    draw = [(mz, it) for (mz, it) in peaks if vmin <= mz <= vmax]
    if not draw:
        draw = peaks
    base = max((it for _, it in draw), default=1.0) or 1.0

    def mz2x(mz):
        return x0 + (mz - vmin) / (vmax - vmin) * (x1 - x0)

    if st["show_yticks"]:
        for i in range(6):
            frac = i / 5.0
            y = y1 - frac * (y1 - y0)
            be.line(x0 - 5, y, x0, y, axis, 1)
            lab = f"{frac * 100:.0f}%" if relative else _fmt_si(frac * base)
            be.text(x0 - 10, y, lab, col(st["tick_color"], grey),
                    size=fs(st["axis_fontsize"]), anchor="rm",
                    bold=st["tick_bold"], italic=st["tick_italic"])

    if st["show_xticks"]:
        step = _nice_step((vmax - vmin) / 8.0)
        v = (int(vmin / step) + 1) * step
        while v < vmax:
            x = mz2x(v)
            be.line(x, y1, x, y1 + 5, axis, 1)
            be.text(x, y1 + 8, f"{v:g}", col(st["tick_color"], grey),
                    size=fs(st["axis_fontsize"]), anchor="ct",
                    bold=st["tick_bold"], italic=st["tick_italic"])
            v += step
    if st["show_titles"]:
        xt = st["xtitle"] or "m/z"
        yt = st["ytitle"] or ("Relative intensity" if relative
                              else "Intensity")
        be.text((x0 + x1) / 2, y1 + 42, xt, col(st["xtitle_color"], grey),
                size=fs(st["xtitle_fontsize"]), anchor="ct",
                bold=st["xtitle_bold"], italic=st["xtitle_italic"])
        be.text(24, (y0 + y1) / 2, yt, col(st["ytitle_color"], grey),
                size=fs(st["ytitle_fontsize"]), anchor="cm", rotate=90,
                bold=st["ytitle_bold"], italic=st["ytitle_italic"])

    for mz, it in draw:
        x = mz2x(mz)
        frac = it / base
        y = y1 - frac * (y1 - y0)
        be.line(x, y1, x, y,
                _blend_white(scheme_color(scheme, frac), transparency),
                bar_width)

    if st["show_peaklabels"]:
        gap = max(40.0, (x1 - x0) / 22.0)
        placed = []
        for mz, it in sorted(draw, key=lambda p: p[1], reverse=True):
            x = mz2x(mz)
            if st["peaklabel_auto"]:
                if any(abs(x - px) < gap for px in placed):
                    continue
            elif len(placed) >= max(0, st["peaklabel_n"]):
                break
            placed.append(x)
            y = y1 - (it / base) * (y1 - y0)
            be.text(x, y - 6, f"{mz:.4f}".rstrip("0").rstrip("."),
                    col(st["peaklabel_color"], label),
                    size=fs(st["peaklabel_fontsize"]), anchor="cb",
                    bold=st["peaklabel_bold"], italic=st["peaklabel_italic"])

    draw_annotations()


def render_spectrum_image(spec, relative=True, width=1200, height=700,
                          scheme=DEFAULT_SCHEME, bar_width=2, transparency=0,
                          bg=DEFAULT_BG, limit_prec=False, annotations=None,
                          style=None):
    """Render a spectrum stick-plot to a PIL RGB image (raster)."""
    theme = BG_THEMES.get(bg, BG_THEMES[DEFAULT_BG])
    be = _PILBackend(width, height)
    _draw_spectrum(be, width, height, spec, relative, scheme,
                   bar_width, transparency, theme, limit_prec,
                   annotations, style)
    return be.img


def save_spectra_pdf(specs, path, relative=True, scheme=DEFAULT_SCHEME,
                     width=1200, height=700, bar_width=2, transparency=0,
                     bg=DEFAULT_BG, limit_prec=False, annotations=None,
                     style=None):
    """Write one or more spectra to a multi-page vector PDF.

    Text is embedded as Arial (Latin) / Microsoft YaHei (CJK) subsets.
    """
    theme = BG_THEMES.get(bg, BG_THEMES[DEFAULT_BG])
    be = _make_pdf_backend(path, width, height)
    for i, s in enumerate(specs):
        if i > 0:
            be.new_page()
        _draw_spectrum(be, width, height, s, relative, scheme,
                       bar_width, transparency, theme, limit_prec,
                       annotations, style)
    be.save()


# --------------------------------------------------------------------------- #
# Main application
# --------------------------------------------------------------------------- #

class App(_TkBase):
    # special sort options (not real metadata fields)
    SORT_FILE_ORDER = "(file order)"
    SORT_NAME = "Name"
    SORT_NUM_PEAKS = "# peaks"

    # list-label options
    LABEL_AUTO = "(auto)"

    def __init__(self, initial_files=None):
        super().__init__()
        self.title("fancyMS2viewer")
        self.geometry("1100x680")
        self.minsize(820, 500)

        # application icon (falls back to a blank icon if the file is missing)
        try:
            icon_path = _resource_path("icon.png")
            if os.path.isfile(icon_path):
                self._app_icon = tk.PhotoImage(file=icon_path)
            else:
                self._app_icon = tk.PhotoImage(width=1, height=1)
            self.iconphoto(True, self._app_icon)
        except Exception:  # noqa: BLE001
            pass

        # UI font: Arial for Latin (Windows font-linking covers CJK fallback)
        import tkinter.font as tkfont
        for _fn in ("TkDefaultFont", "TkTextFont", "TkMenuFont",
                    "TkHeadingFont", "TkIconFont"):
            try:
                tkfont.nametofont(_fn).configure(family="Arial")
            except Exception:  # noqa: BLE001
                pass

        self.spectra = []          # all loaded spectra
        self.filtered = []         # currently shown (after search)
        self.hidden_fields = set()  # metadata field names hidden from the panel

        # shared display settings (used by menu + plot + exporters)
        self.relative_var = tk.BooleanVar(value=True)
        self.scheme_name = tk.StringVar(value=DEFAULT_SCHEME)
        self.bar_width = tk.IntVar(value=2)          # stick thickness (px)
        self.bar_trans = tk.IntVar(value=0)          # transparency % (0=opaque)
        self.bg_scheme = tk.StringVar(value=DEFAULT_BG)   # plot background
        self.limit_prec_var = tk.BooleanVar(value=True)   # peaks <= prec + 5

        self._build_ui()
        self._build_menu()
        self._setup_dnd()

        self.plot.set_scheme(self.scheme_name.get())
        self.plot.set_bar_style(self.bar_width.get(), self.bar_trans.get())
        self.plot.set_bg_theme(self.bg_scheme.get())
        self.plot.limit_prec = self.limit_prec_var.get()
        self.plot.on_context = self._plot_context_menu

        if initial_files:
            self.load_files(initial_files)

    # -- drag & drop --------------------------------------------------------
    def _setup_dnd(self):
        """Register the window (and main widgets) as file drop targets."""
        if not _DND_AVAILABLE:
            return
        # Register on the areas a user is likely to drop onto. In tkdnd a
        # drop only fires on a registered widget, so cover the big ones.
        for widget in (self, self.listbox, self.plot, self.meta_tree):
            try:
                widget.drop_target_register(DND_FILES)
                widget.dnd_bind("<<Drop>>", self._on_drop)
            except Exception:  # noqa: BLE001
                pass

    def _on_drop(self, event):
        # event.data is a Tcl list of paths; brace-wrapped if they contain
        # spaces. splitlist() parses it correctly.
        try:
            paths = self.tk.splitlist(event.data)
        except Exception:  # noqa: BLE001
            paths = event.data.split()
        files = [p for p in paths if os.path.isfile(p)]
        if files:
            self.load_files(files)
        else:
            self.status.config(text="Dropped item is not a readable file.")
        return event.action

    # -- UI construction ----------------------------------------------------
    def _build_menu(self):
        menubar = tk.Menu(self)
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label="Open… (.msp/.mgf/.mzML)",
                             command=self.open_dialog, accelerator="Ctrl+O")
        filemenu.add_command(label="Clear all", command=self.clear_all)
        filemenu.add_separator()
        filemenu.add_command(label="Save spectra as .msp…",
                             command=lambda: self.save_spectra("msp"))
        filemenu.add_command(label="Save spectra as .mgf…",
                             command=lambda: self.save_spectra("mgf"))
        filemenu.add_command(label="Save plot as PNG…",
                             command=lambda: self.save_image("png"))
        filemenu.add_command(label="Save plot as PDF… (vector)",
                             command=lambda: self.save_image("pdf"))
        filemenu.add_separator()
        filemenu.add_command(
            label="Copy peaks to clipboard (m/z + intensity)",
            command=self.copy_peaks, accelerator="Ctrl+C")
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=self.destroy)
        menubar.add_cascade(label="File", menu=filemenu)

        # -- Settings menu --
        setmenu = tk.Menu(menubar, tearoff=0)
        colormenu = tk.Menu(setmenu, tearoff=0)
        for name in COLOR_SCHEMES:
            if name == CUSTOM_SCHEME:
                continue
            colormenu.add_radiobutton(label=name, value=name,
                                      variable=self.scheme_name,
                                      command=self._on_scheme_change)
        colormenu.add_radiobutton(label=CUSTOM_SCHEME, value=CUSTOM_SCHEME,
                                  variable=self.scheme_name,
                                  command=self._on_scheme_change)
        colormenu.add_separator()
        colormenu.add_command(label="Custom colour / gradient…",
                              command=self.edit_custom_color)
        setmenu.add_cascade(label="Peak colour", menu=colormenu)
        setmenu.add_command(label="Bar style (thickness / transparency)…",
                            command=self.edit_bar_style)
        bgmenu = tk.Menu(setmenu, tearoff=0)
        for name in BG_THEMES:
            bgmenu.add_radiobutton(label=name, value=name,
                                   variable=self.bg_scheme,
                                   command=self._on_bg_change)
        setmenu.add_cascade(label="Plot background", menu=bgmenu)
        setmenu.add_separator()
        setmenu.add_checkbutton(
            label="Show relative intensity (%)  (else absolute)",
            variable=self.relative_var, command=self._toggle_relative)
        setmenu.add_checkbutton(
            label="Show only peaks ≤ precursor + 5 m/z",
            variable=self.limit_prec_var, command=self._on_limit_change)
        setmenu.add_separator()
        setmenu.add_command(label="Metadata fields to show…",
                            command=self.choose_meta_fields)
        menubar.add_cascade(label="Settings", menu=setmenu)

        # -- Plot menu: metadata labels on the figure + axis styling --
        plotmenu = tk.Menu(menubar, tearoff=0)
        lblmenu = tk.Menu(plotmenu, tearoff=0)
        lblmenu.add_command(label="Add label…", command=self.add_label_dialog)
        lblmenu.add_command(label="Clear all labels",
                            command=lambda: self.plot.clear_annotations())
        lblmenu.add_separator()
        lblmenu.add_command(
            label="On the plot: drag = move · wheel = resize · "
                  "right-click = edit/delete", state="disabled")
        plotmenu.add_cascade(label="Metadata labels", menu=lblmenu)
        plotmenu.add_command(label="Axes & tick labels…",
                             command=self.edit_axis_style)
        menubar.add_cascade(label="Plot", menu=plotmenu)

        self.config(menu=menubar)
        self.bind("<Control-o>", lambda e: self.open_dialog())
        self.bind("<Control-s>", lambda e: self.save_spectra("msp"))
        self.bind("<Control-c>", lambda e: self.copy_peaks())

    def _on_scheme_change(self):
        self.plot.set_scheme(self.scheme_name.get())

    def _on_bg_change(self):
        self.plot.set_bg_theme(self.bg_scheme.get())

    def _on_limit_change(self):
        self.plot.set_limit_prec(self.limit_prec_var.get())

    def add_label_dialog(self):
        """Add a draggable metadata label to the plot."""
        fields = ["Name"] + [f for f in self._all_meta_fields()
                             if f not in ("# peaks", "source file", "Name")]
        dlg = tk.Toplevel(self)
        dlg.title("Add label")
        dlg.transient(self)
        dlg.resizable(False, False)
        ttk.Label(dlg, text="Metadata field:").grid(row=0, column=0,
                                                    sticky="w", padx=10,
                                                    pady=(12, 4))
        fvar = tk.StringVar(value=fields[0])
        ttk.Combobox(dlg, textvariable=fvar, values=fields, width=26,
                     state="readonly").grid(row=0, column=1, padx=(0, 10),
                                            pady=(12, 4))
        mode = tk.StringVar(value="value")
        ttk.Radiobutton(dlg, text="Value only", value="value",
                        variable=mode).grid(row=1, column=1, sticky="w")
        ttk.Radiobutton(dlg, text="Field: value", value="field",
                        variable=mode).grid(row=2, column=1, sticky="w")
        btns = ttk.Frame(dlg)
        btns.grid(row=3, column=0, columnspan=2, sticky="e", padx=10,
                  pady=(8, 12))

        def add():
            self.plot.add_annotation(fvar.get(),
                                     show_field=(mode.get() == "field"))
            dlg.destroy()
        ttk.Button(btns, text="Add", command=add).pack(side="right")
        ttk.Button(btns, text="Cancel", command=dlg.destroy).pack(
            side="right", padx=(0, 6))

    def edit_axis_style(self):
        """Edit axis ticks / titles / peak labels: show, content & size."""
        p = self.plot
        dlg = tk.Toplevel(self)
        dlg.title("Axes & tick labels")
        dlg.transient(self)
        dlg.resizable(False, False)
        v_xt = tk.BooleanVar(value=p.show_xticks)
        v_yt = tk.BooleanVar(value=p.show_yticks)
        v_ti = tk.BooleanVar(value=p.show_titles)
        v_pl = tk.BooleanVar(value=p.show_peaklabels)
        v_auto = tk.BooleanVar(value=p.peaklabel_auto)
        n_pl = tk.IntVar(value=p.peaklabel_n)
        f_ax = tk.IntVar(value=p.axis_fontsize)
        f_xti = tk.IntVar(value=p.xtitle_fontsize)
        f_yti = tk.IntVar(value=p.ytitle_fontsize)
        f_pl = tk.IntVar(value=p.peaklabel_fontsize)
        xt = tk.StringVar(value=p.xtitle)
        yt = tk.StringVar(value=p.ytitle)
        b_tk = tk.BooleanVar(value=p.tick_bold)
        i_tk = tk.BooleanVar(value=p.tick_italic)
        b_xti = tk.BooleanVar(value=p.xtitle_bold)
        i_xti = tk.BooleanVar(value=p.xtitle_italic)
        b_yti = tk.BooleanVar(value=p.ytitle_bold)
        i_yti = tk.BooleanVar(value=p.ytitle_italic)
        b_pl = tk.BooleanVar(value=p.peaklabel_bold)
        i_pl = tk.BooleanVar(value=p.peaklabel_italic)
        c_tk = tk.StringVar(value=p.tick_color)
        c_xti = tk.StringVar(value=p.xtitle_color)
        c_yti = tk.StringVar(value=p.ytitle_color)
        c_pl = tk.StringVar(value=p.peaklabel_color)

        def apply(*_a):
            p.set_axis_style(
                show_xticks=v_xt.get(), show_yticks=v_yt.get(),
                show_titles=v_ti.get(), show_peaklabels=v_pl.get(),
                peaklabel_auto=v_auto.get(), peaklabel_n=n_pl.get(),
                axis_fontsize=max(5, f_ax.get()),
                xtitle_fontsize=max(5, f_xti.get()),
                ytitle_fontsize=max(5, f_yti.get()),
                peaklabel_fontsize=max(5, f_pl.get()),
                xtitle=xt.get(), ytitle=yt.get(),
                tick_bold=b_tk.get(), tick_italic=i_tk.get(),
                xtitle_bold=b_xti.get(), xtitle_italic=i_xti.get(),
                ytitle_bold=b_yti.get(), ytitle_italic=i_yti.get(),
                peaklabel_bold=b_pl.get(), peaklabel_italic=i_pl.get(),
                tick_color=c_tk.get(), xtitle_color=c_xti.get(),
                ytitle_color=c_yti.get(), peaklabel_color=c_pl.get())

        pad = dict(padx=12, pady=2)
        ttk.Label(dlg, text="Show:", font=("Arial", 9, "bold")).pack(
            anchor="w", padx=12, pady=(12, 2))
        ttk.Checkbutton(dlg, text="X-axis tick labels", variable=v_xt,
                        command=apply).pack(anchor="w", **pad)
        ttk.Checkbutton(dlg, text="Y-axis tick labels", variable=v_yt,
                        command=apply).pack(anchor="w", **pad)
        ttk.Checkbutton(dlg, text="Axis titles", variable=v_ti,
                        command=apply).pack(anchor="w", **pad)
        ttk.Checkbutton(dlg, text="Peak m/z labels", variable=v_pl,
                        command=apply).pack(anchor="w", **pad)

        ttk.Separator(dlg).pack(fill="x", padx=10, pady=6)
        ttk.Label(dlg, text="Axis titles (blank = default):",
                  font=("Arial", 9, "bold")).pack(anchor="w", padx=12)
        for lab, var in (("X title", xt), ("Y title", yt)):
            row = ttk.Frame(dlg)
            row.pack(anchor="w", **pad)
            ttk.Label(row, text=lab, width=7).pack(side="left")
            e = ttk.Entry(row, textvariable=var, width=26)
            e.pack(side="left")
            e.bind("<KeyRelease>", apply)

        ttk.Separator(dlg).pack(fill="x", padx=10, pady=6)
        ttk.Label(dlg, text="Sizes & counts:",
                  font=("Arial", 9, "bold")).pack(anchor="w", padx=12)
        for lab, var, lo, hi in (("Tick font", f_ax, 5, 20),
                                 ("X-title font", f_xti, 5, 24),
                                 ("Y-title font", f_yti, 5, 24),
                                 ("Peak-label font", f_pl, 5, 20)):
            row = ttk.Frame(dlg)
            row.pack(anchor="w", **pad)
            ttk.Label(row, text=lab, width=15).pack(side="left")
            tk.Scale(row, from_=lo, to=hi, orient="horizontal", length=140,
                     variable=var, command=apply).pack(side="left")
        ttk.Checkbutton(dlg, text="Peak labels: auto count (by spacing)",
                        variable=v_auto, command=apply).pack(anchor="w", **pad)
        row = ttk.Frame(dlg)
        row.pack(anchor="w", **pad)
        ttk.Label(row, text="…or fixed count", width=15).pack(side="left")
        tk.Scale(row, from_=0, to=40, orient="horizontal", length=140,
                 variable=n_pl, command=apply).pack(side="left")

        ttk.Separator(dlg).pack(fill="x", padx=10, pady=6)
        ttk.Label(dlg, text="Bold / Italic:",
                  font=("Arial", 9, "bold")).pack(anchor="w", padx=12)
        for lab, bv, iv in (("Tick labels", b_tk, i_tk),
                            ("X title", b_xti, i_xti),
                            ("Y title", b_yti, i_yti),
                            ("Peak labels", b_pl, i_pl)):
            row = ttk.Frame(dlg)
            row.pack(anchor="w", **pad)
            ttk.Label(row, text=lab, width=12).pack(side="left")
            ttk.Checkbutton(row, text="Bold", variable=bv,
                            command=apply).pack(side="left")
            ttk.Checkbutton(row, text="Italic", variable=iv,
                            command=apply).pack(side="left", padx=(8, 0))

        ttk.Separator(dlg).pack(fill="x", padx=10, pady=6)
        ttk.Label(dlg, text="Colours (Default = follow background):",
                  font=("Arial", 9, "bold")).pack(anchor="w", padx=12)
        th = p.theme
        for lab, var, dflt in (("Tick labels", c_tk, th["text"]),
                               ("X title", c_xti, th["text"]),
                               ("Y title", c_yti, th["text"]),
                               ("Peak labels", c_pl, th["label"])):
            _color_row(dlg, lab, var, dflt, self._pick_color, apply).pack(
                anchor="w", **pad)
        ttk.Button(dlg, text="Close", command=dlg.destroy).pack(
            anchor="e", padx=12, pady=(6, 12))

    def copy_peaks(self):
        """Copy the current spectrum's peaks to the clipboard as two columns.

        Tab-separated m/z and intensity, one peak per line — pastes straight
        into Excel or a text editor.
        """
        spec = getattr(self, "_current_spec", None)
        if spec is None:
            sel = self._selected_spectra()
            spec = sel[0] if sel else None
        if spec is None or not len(spec.mz):
            self._alert("Copy peaks", "No spectrum with peaks is selected.")
            return
        text = "\n".join(f"{_fmt_mz(mz)}\t{it:g}" for mz, it in spec.peaks)
        self.clipboard_clear()
        self.clipboard_append(text)
        self.status.config(
            text=f"Copied {len(spec.mz)} peaks (m/z + intensity) "
                 f"to clipboard.")

    def edit_bar_style(self):
        """Adjust stick thickness and transparency (applied live)."""
        dlg = tk.Toplevel(self)
        dlg.title("Bar style")
        dlg.transient(self)
        dlg.resizable(False, False)

        def apply_live(*_a):
            self.plot.set_bar_style(self.bar_width.get(),
                                    self.bar_trans.get())

        ttk.Label(dlg, text="Thickness (px)").grid(row=0, column=0, sticky="w",
                                                   padx=10, pady=(12, 2))
        tk.Scale(dlg, from_=1, to=10, orient="horizontal", length=220,
                 variable=self.bar_width, command=apply_live).grid(
            row=0, column=1, padx=(4, 10), pady=(12, 2))
        ttk.Label(dlg, text="Transparency (%)").grid(row=1, column=0,
                                                     sticky="w", padx=10)
        tk.Scale(dlg, from_=0, to=95, orient="horizontal", length=220,
                 variable=self.bar_trans, command=apply_live).grid(
            row=1, column=1, padx=(4, 10))
        ttk.Label(dlg, text="0 = opaque (default)",
                  foreground="#888").grid(row=2, column=1, sticky="w",
                                          padx=(4, 10))
        ttk.Button(dlg, text="Close", command=dlg.destroy).grid(
            row=3, column=1, sticky="e", padx=(4, 10), pady=(6, 12))

    def _ask_image_size(self, specs=None):
        """Ask for export image size (current or custom W×H) with a live preview.

        Returns (width, height) or None if cancelled.
        """
        from PIL import ImageTk
        pw = max(int(self.plot.winfo_width()), 100)
        ph = max(int(self.plot.winfo_height()), 100)
        spec = (specs or [getattr(self, "_current_spec", None)])[0]
        params = self._export_params()

        dlg = tk.Toplevel(self)
        dlg.title("Image size")
        dlg.transient(self)
        dlg.resizable(False, False)
        result = {"size": None}
        mode = tk.StringVar(value="preview")
        wv = tk.StringVar(value="1200")
        hv = tk.StringVar(value="700")

        left = ttk.Frame(dlg)
        left.pack(side="left", padx=12, pady=12, anchor="n")
        ttk.Radiobutton(left, text=f"Current preview size ({pw} × {ph})",
                        value="preview", variable=mode).grid(
            row=0, column=0, columnspan=3, sticky="w", pady=(0, 4))
        ttk.Radiobutton(left, text="Custom size:", value="custom",
                        variable=mode).grid(row=1, column=0, sticky="w")
        ttk.Entry(left, textvariable=wv, width=7).grid(row=1, column=1)
        ttk.Label(left, text="×").grid(row=1, column=2)
        ttk.Entry(left, textvariable=hv, width=7).grid(row=2, column=1,
                                                       sticky="w")

        preview = ttk.Frame(dlg)
        preview.pack(side="left", padx=(0, 12), pady=12)
        ttk.Label(preview, text="Preview", foreground="#888").pack()
        prev_cv = tk.Canvas(preview, width=400, height=250,
                            background="#fafafa", highlightthickness=1,
                            highlightbackground="#ccc")
        prev_cv.pack()

        def chosen():
            if mode.get() == "preview":
                return (pw, ph)
            w, h = _to_float(wv.get()), _to_float(hv.get())
            if w and h and w >= 50 and h >= 50:
                return (int(w), int(h))
            return None

        def show_msg(msg):
            prev_cv.delete("all")
            prev_cv.create_text(200, 125, text=msg, fill="#999",
                                font=("Arial", 10))

        def update_preview():
            sz = chosen()
            if sz is None or spec is None:
                show_msg("(enter width & height)")
                return
            w, h = sz
            try:
                img = render_spectrum_image(spec, width=w, height=h, **params)
            except Exception:  # noqa: BLE001
                show_msg("(cannot render)")
                return
            img.thumbnail((392, 242))
            self._size_preview_photo = ImageTk.PhotoImage(img)
            prev_cv.delete("all")
            prev_cv.create_image(200, 125, image=self._size_preview_photo,
                                 anchor="center")
            prev_cv.create_text(200, 244, text=f"{w} × {h} px", fill="#888",
                                font=("Arial", 8), anchor="s")

        after = {"id": None}

        def schedule(*_a):
            if after["id"]:
                try:
                    dlg.after_cancel(after["id"])
                except Exception:  # noqa: BLE001
                    pass
            after["id"] = dlg.after(250, update_preview)

        for var in (wv, hv, mode):
            var.trace_add("write", schedule)

        btns = ttk.Frame(dlg)
        btns.pack(side="bottom", anchor="e", padx=12, pady=(0, 12))

        def ok():
            sz = chosen()
            if sz is None:
                self._alert("Image size", "Enter valid width and height "
                                          "(≥ 50 px).")
                return
            result["size"] = sz
            dlg.destroy()
        ttk.Button(btns, text="OK", command=ok).pack(side="right")
        ttk.Button(btns, text="Cancel", command=dlg.destroy).pack(
            side="right", padx=(0, 6))

        dlg.update_idletasks()
        root = self.winfo_toplevel()
        x = root.winfo_rootx() + (root.winfo_width() - dlg.winfo_width()) // 2
        y = root.winfo_rooty() + (root.winfo_height() - dlg.winfo_height()) // 2
        dlg.geometry(f"+{max(x, 0)}+{max(y, 0)}")
        update_preview()
        dlg.grab_set()
        dlg.wait_window()
        return result["size"]

    # -- right-click context menus -----------------------------------------
    def _plot_context_menu(self, event):
        """Right-click on the plot: export the current spectrum (single)."""
        spec = getattr(self, "_current_spec", None)
        m = tk.Menu(self, tearoff=0)
        m.add_command(label="Reset zoom", command=self.plot.reset_view)
        m.add_separator()
        if spec is not None:
            one = [spec]
            m.add_command(label="Export plot as PNG…",
                          command=lambda: self.save_image("png", one))
            m.add_command(label="Export plot as PDF (vector)…",
                          command=lambda: self.save_image("pdf", one))
            m.add_separator()
            m.add_command(label="Export spectrum as .msp…",
                          command=lambda: self.save_spectra("msp", one))
            m.add_command(label="Export spectrum as .mgf…",
                          command=lambda: self.save_spectra("mgf", one))
            m.add_command(label="Copy peaks (m/z + intensity)",
                          command=self.copy_peaks)
        try:
            m.tk_popup(event.x_root, event.y_root)
        finally:
            m.grab_release()

    def _list_context_menu(self, event):
        """Right-click on the list: export the selection.

        To avoid accidental bulk work, images can only be exported for a
        single selected spectrum; multiple selections may only be exported to
        one .msp / .mgf file.
        """
        idx = self.listbox.nearest(event.y)
        if idx not in self.listbox.curselection():
            self.listbox.selection_clear(0, "end")
            self.listbox.selection_set(idx)
            self.listbox.activate(idx)
            if 0 <= idx < len(self.filtered):
                self._show_spectrum(self.filtered[idx])
        specs = self._selected_spectra()
        if not specs:
            return
        m = tk.Menu(self, tearoff=0)
        n = len(specs)
        if n == 1:
            m.add_command(label="Export plot as PNG…",
                          command=lambda: self.save_image("png", specs))
            m.add_command(label="Export plot as PDF (vector)…",
                          command=lambda: self.save_image("pdf", specs))
            m.add_separator()
            m.add_command(label="Export spectrum as .msp…",
                          command=lambda: self.save_spectra("msp", specs))
            m.add_command(label="Export spectrum as .mgf…",
                          command=lambda: self.save_spectra("mgf", specs))
            m.add_command(label="Copy peaks (m/z + intensity)",
                          command=self.copy_peaks)
        else:
            m.add_command(label=f"Export {n} spectra to one .msp…",
                          command=lambda: self.save_spectra("msp", specs))
            m.add_command(label=f"Export {n} spectra to one .mgf…",
                          command=lambda: self.save_spectra("mgf", specs))
            m.add_separator()
            m.add_command(label="(image export disabled for multiple)",
                          state="disabled")
        try:
            m.tk_popup(event.x_root, event.y_root)
        finally:
            m.grab_release()

    # -- small English replacements for the OS-localised native dialogs -----
    def _alert(self, title, text):
        """Simple English message dialog (avoids the OS-localised messagebox)."""
        dlg = tk.Toplevel(self)
        dlg.title(title)
        dlg.transient(self)
        dlg.resizable(False, False)
        ttk.Label(dlg, text=text, wraplength=380,
                  justify="left").pack(padx=20, pady=(18, 10))
        ttk.Button(dlg, text="OK", command=dlg.destroy).pack(pady=(0, 14))
        dlg.bind("<Return>", lambda e: dlg.destroy())
        dlg.update_idletasks()
        x = self.winfo_rootx() + (self.winfo_width() - dlg.winfo_width()) // 2
        y = self.winfo_rooty() + (self.winfo_height() - dlg.winfo_height()) // 2
        dlg.geometry(f"+{max(x, 0)}+{max(y, 0)}")
        try:
            dlg.grab_set()
        except Exception:  # noqa: BLE001
            pass

    def _pick_color(self, initial=(255, 0, 0)):
        """English RGB colour picker (replaces the OS-localised colour dialog).

        Returns an (r, g, b) tuple, or None if cancelled.
        """
        dlg = tk.Toplevel(self)
        dlg.title("Choose colour")
        dlg.transient(self)
        dlg.resizable(False, False)
        result = {"rgb": None}
        rv = tk.IntVar(value=initial[0])
        gv = tk.IntVar(value=initial[1])
        bv = tk.IntVar(value=initial[2])
        hexv = tk.StringVar()

        def cur():
            return (rv.get(), gv.get(), bv.get())

        def refresh(*_a):
            c = cur()
            swatch.configure(background=_hex(c))
            hexv.set("%02X%02X%02X" % c)

        # -- colour palette (色盘): hue across x, white→vivid→black down y --
        import colorsys
        from PIL import Image, ImageTk
        PW, PH = 224, 120
        pimg = Image.new("RGB", (PW, PH))
        ppx = pimg.load()
        for x in range(PW):
            h = x / (PW - 1)
            for y in range(PH):
                t = y / (PH - 1)
                if t < 0.5:
                    s, v = t * 2.0, 1.0            # white → vivid
                else:
                    s, v = 1.0, 1.0 - (t - 0.5) * 2.0   # vivid → black
                r, g, b = colorsys.hsv_to_rgb(h, s, v)
                ppx[x, y] = (int(r * 255), int(g * 255), int(b * 255))
        self._pal_photo = ImageTk.PhotoImage(pimg)   # keep a reference
        pal = tk.Canvas(dlg, width=PW, height=PH, highlightthickness=1,
                        highlightbackground="#999", cursor="crosshair")
        pal.create_image(0, 0, anchor="nw", image=self._pal_photo)
        pal.grid(row=0, column=0, columnspan=2, padx=(10, 6), pady=(12, 6))

        def on_pal(ev):
            x = min(max(int(ev.x), 0), PW - 1)
            y = min(max(int(ev.y), 0), PH - 1)
            r, g, b = pimg.getpixel((x, y))
            rv.set(r)
            gv.set(g)
            bv.set(b)
            refresh()
        pal.bind("<Button-1>", on_pal)
        pal.bind("<B1-Motion>", on_pal)

        swatch = tk.Label(dlg, width=8, height=6, relief="solid",
                          borderwidth=1)
        swatch.grid(row=0, column=2, padx=(0, 10), pady=(12, 6))

        for i, (lab, var) in enumerate((("R", rv), ("G", gv), ("B", bv))):
            ttk.Label(dlg, text=lab).grid(row=1 + i, column=0, padx=(10, 2))
            tk.Scale(dlg, from_=0, to=255, orient="horizontal", length=220,
                     variable=var, command=refresh).grid(row=1 + i, column=1,
                                                         sticky="we")
            ttk.Label(dlg, textvariable=var, width=4).grid(row=1 + i, column=2,
                                                          padx=(2, 10))

        ttk.Label(dlg, text="Hex").grid(row=4, column=0, padx=(10, 2),
                                        pady=(6, 0))
        hexbox = ttk.Entry(dlg, textvariable=hexv, width=10)
        hexbox.grid(row=4, column=1, sticky="w", pady=(6, 0))

        def from_hex(*_a):
            s = hexv.get().strip().lstrip("#")
            if len(s) == 6:
                try:
                    rv.set(int(s[0:2], 16))
                    gv.set(int(s[2:4], 16))
                    bv.set(int(s[4:6], 16))
                    refresh()
                except ValueError:
                    pass
        hexbox.bind("<Return>", from_hex)

        # quick palette
        palette = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
                   "#17becf", "#e377c2", "#000000", "#7f7f7f", "#ffffff"]
        palrow = ttk.Frame(dlg)
        palrow.grid(row=5, column=0, columnspan=3, padx=10, pady=8)
        for hx in palette:
            r = int(hx[1:3], 16)
            g = int(hx[3:5], 16)
            b = int(hx[5:7], 16)

            def setc(rr=r, gg=g, bb=b):
                rv.set(rr)
                gv.set(gg)
                bv.set(bb)
                refresh()
            tk.Button(palrow, background=hx, width=2, relief="solid",
                      borderwidth=1, command=setc).pack(side="left", padx=1)

        btns = ttk.Frame(dlg)
        btns.grid(row=6, column=0, columnspan=3, sticky="e", padx=10,
                  pady=(4, 12))

        def ok():
            result["rgb"] = cur()
            dlg.destroy()
        ttk.Button(btns, text="OK", command=ok).pack(side="right")
        ttk.Button(btns, text="Cancel", command=dlg.destroy).pack(
            side="right", padx=(0, 6))

        refresh()
        dlg.grab_set()
        dlg.wait_window()
        return result["rgb"]

    def edit_custom_color(self):
        """Build a custom single colour or gradient with adjustable stops."""
        dlg = tk.Toplevel(self)
        dlg.title("Custom peak colour")
        dlg.transient(self)
        dlg.resizable(False, False)

        # working copy: list of [pos(0..1), [r,g,b]]
        stops = [[p, list(c)] for p, c in _positioned(
            COLOR_SCHEMES[CUSTOM_SCHEME])]
        mode = tk.StringVar(value="single" if len(stops) == 1 else "gradient")

        ttk.Label(dlg, text="Mode:").grid(row=0, column=0, sticky="w",
                                          padx=10, pady=(12, 4))
        modebar = ttk.Frame(dlg)
        modebar.grid(row=0, column=1, sticky="w", pady=(12, 4))
        ttk.Radiobutton(modebar, text="Single colour", value="single",
                        variable=mode).pack(side="left")
        ttk.Radiobutton(modebar, text="Gradient (low → high)",
                        value="gradient", variable=mode).pack(side="left",
                                                             padx=(8, 0))

        rows_frame = ttk.Frame(dlg)
        rows_frame.grid(row=1, column=0, columnspan=2, sticky="we",
                        padx=10, pady=4)
        preview = tk.Canvas(dlg, height=28, width=360, highlightthickness=1,
                            highlightbackground="#999")
        preview.grid(row=2, column=0, columnspan=2, sticky="we",
                     padx=10, pady=(4, 8))

        def redraw_preview():
            preview.delete("all")
            w = preview.winfo_width() or 360
            data = ([[0.0, stops[0][1]]] if mode.get() == "single"
                    else [[s[0], s[1]] for s in stops])
            for x in range(w):
                frac = x / max(w - 1, 1)
                preview.create_line(x, 0, x, 28,
                                    fill=_hex(_interp_stops(data, frac)))

        def rebuild():
            for c in rows_frame.winfo_children():
                c.destroy()
            gradient = mode.get() == "gradient"
            shown = stops if gradient else stops[:1]
            for idx, (pos, rgb) in enumerate(shown):
                row = ttk.Frame(rows_frame)
                row.pack(fill="x", pady=2)
                sw = tk.Label(row, width=4, background=_hex(tuple(rgb)),
                              relief="solid", borderwidth=1)
                sw.pack(side="left")

                def pick(i=idx):
                    c = self._pick_color(tuple(stops[i][1]))
                    if c:
                        stops[i][1] = list(c)
                        rebuild()
                        redraw_preview()
                ttk.Button(row, text="Pick…", command=pick).pack(
                    side="left", padx=6)

                if gradient:
                    ttk.Label(row, text="pos").pack(side="left")
                    pvar = tk.IntVar(value=int(round(pos * 100)))

                    def moved(val, i=idx, v=None):
                        stops[i][0] = int(float(val)) / 100.0
                        redraw_preview()
                    tk.Scale(row, from_=0, to=100, orient="horizontal",
                             length=140, variable=pvar,
                             command=moved).pack(side="left", padx=(2, 6))
                    if len(stops) > 2:
                        def rem(i=idx):
                            del stops[i]
                            rebuild()
                            redraw_preview()
                        ttk.Button(row, text="✕", width=3,
                                   command=rem).pack(side="left")

            if gradient and len(stops) < 8:
                def add():
                    stops.append([1.0, list(stops[-1][1])])
                    rebuild()
                    redraw_preview()
                ttk.Button(rows_frame, text="+ add stop",
                           command=add).pack(anchor="w", pady=(4, 0))

        mode.trace_add("write", lambda *a: (rebuild(), redraw_preview()))
        preview.bind("<Configure>", lambda e: redraw_preview())

        btns = ttk.Frame(dlg)
        btns.grid(row=3, column=0, columnspan=2, sticky="e", padx=10,
                  pady=(4, 12))

        def apply():
            if mode.get() == "single":
                COLOR_SCHEMES[CUSTOM_SCHEME] = [tuple(stops[0][1])]
            else:
                COLOR_SCHEMES[CUSTOM_SCHEME] = [
                    (s[0], tuple(s[1])) for s in stops]
            self.scheme_name.set(CUSTOM_SCHEME)
            self.plot.set_scheme(CUSTOM_SCHEME)
            dlg.destroy()

        ttk.Button(btns, text="Apply", command=apply).pack(side="right")
        ttk.Button(btns, text="Cancel", command=dlg.destroy).pack(
            side="right", padx=(0, 6))

        dlg.columnconfigure(1, weight=1)
        rebuild()
        dlg.update_idletasks()
        redraw_preview()

    def _build_ui(self):
        # top bar (two rows: search, then sort)
        top = ttk.Frame(self, padding=(8, 6))
        top.pack(side="top", fill="x")

        # -- row 1: search controls (match mode + add button) --
        row1 = ttk.Frame(top)
        row1.pack(side="top", fill="x")

        ttk.Label(row1, text="Search — match").pack(side="left")
        self.match_mode = tk.StringVar(value="all")
        ttk.Radiobutton(row1, text="All (AND)", value="all",
                        variable=self.match_mode,
                        command=self.apply_filter).pack(side="left",
                                                        padx=(6, 0))
        ttk.Radiobutton(row1, text="Any (OR)", value="any",
                        variable=self.match_mode,
                        command=self.apply_filter).pack(side="left",
                                                        padx=(2, 8))
        ttk.Button(row1, text="+ Add search box",
                   command=self.add_search_row).pack(side="left")

        ttk.Label(row1, text="(display options in the Settings menu)",
                  foreground="#999").pack(side="right")

        # container that holds the dynamic list of search boxes
        self.search_rows = []          # list of dicts: field_var, text_var, frame
        self.search_container = ttk.Frame(top)
        self.search_container.pack(side="top", fill="x", pady=(4, 0))
        self.add_search_row()          # start with one box

        # -- row 2: sort --
        row2 = ttk.Frame(top)
        row2.pack(side="top", fill="x", pady=(6, 0))

        ttk.Label(row2, text="Sort by:").pack(side="left")
        self.sort_var = tk.StringVar(value=self.SORT_FILE_ORDER)
        self.sort_combo = ttk.Combobox(row2, textvariable=self.sort_var,
                                        width=22, state="readonly",
                                        values=[self.SORT_FILE_ORDER])
        self.sort_combo.pack(side="left", padx=(4, 8))
        self.sort_combo.bind("<<ComboboxSelected>>",
                             lambda e: self.apply_filter())

        self.sort_desc_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(row2, text="Descending",
                        variable=self.sort_desc_var,
                        command=self.apply_filter).pack(side="left")

        ttk.Label(row2, text="Label by:").pack(side="left", padx=(16, 0))
        self.label_var = tk.StringVar(value=self.LABEL_AUTO)
        self.label_combo = ttk.Combobox(row2, textvariable=self.label_var,
                                         width=22, state="readonly",
                                         values=[self.LABEL_AUTO])
        self.label_combo.pack(side="left", padx=(4, 8))
        self.label_combo.bind("<<ComboboxSelected>>",
                              lambda e: self.apply_filter())

        # -- row 3: m/z filter --
        row3 = ttk.Frame(top)
        row3.pack(side="top", fill="x", pady=(6, 0))

        self.mz_enable = tk.BooleanVar(value=False)
        ttk.Checkbutton(row3, text="m/z filter",
                        variable=self.mz_enable).pack(side="left")
        self.mz_enable.trace_add("write", lambda *a: self._on_mz_enable())

        # all options live in a frame that is only shown while enabled
        self.mz_options = ttk.Frame(row3)
        opt = self.mz_options

        # target m/z (typed directly, or computed from formula+adduct)
        self.mz_target = tk.StringVar()
        self.mz_target.trace_add("write", lambda *a: self._on_mz_filter_change())
        ttk.Label(opt, text="m/z").pack(side="left", padx=(8, 2))
        ttk.Entry(opt, textvariable=self.mz_target, width=12).pack(side="left")

        ttk.Label(opt, text="±").pack(side="left", padx=(6, 2))
        self.mz_tol = tk.StringVar(value="10")
        self.mz_tol.trace_add("write", lambda *a: self._on_mz_filter_change())
        ttk.Entry(opt, textvariable=self.mz_tol, width=7).pack(side="left")
        self.mz_tol_unit = tk.StringVar(value="ppm")
        ttk.Combobox(opt, textvariable=self.mz_tol_unit, width=5,
                     state="readonly", values=["ppm", "Da"]).pack(
            side="left", padx=(2, 8))
        self.mz_tol_unit.trace_add("write",
                                   lambda *a: self._on_mz_filter_change())

        ttk.Label(opt, text="match").pack(side="left")
        self.mz_target_kind = tk.StringVar(value="precursor")
        ttk.Combobox(opt, textvariable=self.mz_target_kind, width=10,
                     state="readonly",
                     values=["precursor", "any peak"]).pack(side="left",
                                                            padx=(4, 12))
        self.mz_target_kind.trace_add("write",
                                      lambda *a: self._on_mz_filter_change())

        # formula (+ extra mass) + adduct -> computes the m/z above
        ttk.Label(opt, text="Formula").pack(side="left")
        self.mz_formula = tk.StringVar()
        self.mz_formula.trace_add("write", lambda *a: self._compute_mz())
        ttk.Entry(opt, textvariable=self.mz_formula, width=14).pack(
            side="left", padx=(4, 4))
        ttk.Label(opt, text="+Δmass").pack(side="left")
        self.mz_extra = tk.StringVar()
        self.mz_extra.trace_add("write", lambda *a: self._compute_mz())
        ttk.Entry(opt, textvariable=self.mz_extra, width=9).pack(
            side="left", padx=(4, 8))
        ttk.Label(opt, text="Adduct").pack(side="left")
        self.mz_adduct = tk.StringVar(value="[M+H]+")
        # editable → custom adducts like [2M+H]+, [M+2NH4]2+ can be typed in
        ttk.Combobox(opt, textvariable=self.mz_adduct, width=14,
                     values=ADDUCT_PRESETS).pack(side="left", padx=(4, 0))
        self.mz_adduct.trace_add("write", lambda *a: self._compute_mz())

        # main paned area
        paned = ttk.Panedwindow(self, orient="horizontal")
        paned.pack(side="top", fill="both", expand=True)

        # left: list of spectra
        left = ttk.Frame(paned)
        self.count_label = ttk.Label(left, text="0 spectra")
        self.count_label.pack(side="top", anchor="w", padx=6, pady=(6, 2))
        list_frame = ttk.Frame(left)
        list_frame.pack(side="top", fill="both", expand=True,
                        padx=6, pady=(0, 6))
        self.listbox = tk.Listbox(list_frame, activestyle="dotbox",
                                  exportselection=False, selectmode="extended")
        sb = ttk.Scrollbar(list_frame, orient="vertical",
                           command=self.listbox.yview)
        self.listbox.configure(yscrollcommand=sb.set)
        self.listbox.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")
        self.listbox.bind("<<ListboxSelect>>", self._on_select)
        self.listbox.bind("<Button-3>", self._list_context_menu)
        paned.add(left, weight=1)

        # right: plot (top) + metadata (bottom)
        right = ttk.Panedwindow(paned, orient="vertical")

        plot_frame = ttk.Frame(right)
        self.plot = SpectrumPlot(plot_frame, width=600, height=360)
        self.plot.pack(fill="both", expand=True)
        hint = ttk.Label(plot_frame,
                         text="Drag = zoom box (x/y)   •   Wheel = zoom "
                              "(y over axis)   •   Double-click = reset "
                              "(axis or all)   •   Right-click = menu",
                         foreground="#888")
        hint.pack(side="bottom", anchor="w", padx=4, pady=2)
        right.add(plot_frame, weight=3)

        meta_frame = ttk.Frame(right)
        meta_head = ttk.Frame(meta_frame)
        meta_head.pack(fill="x", padx=6, pady=(4, 0))
        ttk.Label(meta_head, text="Metadata",
                  font=("Arial", 9, "bold")).pack(side="left")
        ttk.Button(meta_head, text="Fields…", width=8,
                   command=self.choose_meta_fields).pack(side="right")
        cols = ("field", "value")
        self.meta_tree = ttk.Treeview(meta_frame, columns=cols,
                                      show="headings", height=8)
        self.meta_tree.heading("field", text="Field")
        self.meta_tree.heading("value", text="Value")
        self.meta_tree.column("field", width=180, anchor="w", stretch=False)
        self.meta_tree.column("value", width=400, anchor="w")
        msb = ttk.Scrollbar(meta_frame, orient="vertical",
                            command=self.meta_tree.yview)
        self.meta_tree.configure(yscrollcommand=msb.set)
        self.meta_tree.pack(side="left", fill="both", expand=True,
                           padx=(6, 0), pady=(0, 6))
        msb.pack(side="right", fill="y", pady=(0, 6))
        right.add(meta_frame, weight=2)

        paned.add(right, weight=3)

        # status bar with an embedded (initially hidden) progress bar
        statusbar = ttk.Frame(self)
        statusbar.pack(side="bottom", fill="x")
        self.progress = ttk.Progressbar(statusbar, mode="indeterminate",
                                        length=160)
        # packed only while loading (see _start_progress / _stop_progress)
        ready = ("Ready — open an .msp / .mgf / .mzML file"
                 + (" or drag one onto the window." if _DND_AVAILABLE
                    else "."))
        self.status = ttk.Label(statusbar, text=ready,
                                relief="sunken", anchor="w")
        self.status.pack(side="left", fill="x", expand=True)

    def _start_progress(self):
        self.progress.pack(side="right", padx=6, pady=2)
        self.progress.start(12)

    def _stop_progress(self):
        try:
            self.progress.stop()
            self.progress.pack_forget()
        except Exception:  # noqa: BLE001
            pass

    # -- dynamic search boxes ----------------------------------------------
    def _current_field_values(self):
        fields = set()
        for s in self.spectra:
            fields.update(s.meta.keys())
        return ["(any field)"] + sorted(fields, key=str.lower)

    # search operators
    OP_CONTAINS = "contains"
    OP_RANGE = "range"
    OP_TOL = "value ± tol"

    def add_search_row(self):
        row = ttk.Frame(self.search_container)
        row.pack(side="top", fill="x", pady=1)

        field_var = tk.StringVar(value="(any field)")
        field_combo = ttk.Combobox(row, textvariable=field_var, width=20,
                                    state="readonly",
                                    values=self._current_field_values())
        field_combo.pack(side="left")
        field_combo.bind("<<ComboboxSelected>>",
                         lambda e: self.apply_filter())

        op_var = tk.StringVar(value=self.OP_CONTAINS)
        op_combo = ttk.Combobox(row, textvariable=op_var, width=12,
                                state="readonly",
                                values=[self.OP_CONTAINS, self.OP_RANGE,
                                        self.OP_TOL])
        op_combo.pack(side="left", padx=(4, 4))

        # value entries: text_var (contains / range-min / value),
        #                text_var2 (range-max / tol)
        text_var = tk.StringVar()
        text_var2 = tk.StringVar()
        text_var.trace_add("write", lambda *a: self.schedule_filter())
        text_var2.trace_add("write", lambda *a: self.schedule_filter())

        entry1 = ttk.Entry(row, textvariable=text_var, width=16)
        entry1.pack(side="left")
        mid_lbl = ttk.Label(row, text="")   # "to" or "±" between the entries
        entry2 = ttk.Entry(row, textvariable=text_var2, width=10)

        entry1.bind("<Return>", lambda e: self.add_search_row())
        entry2.bind("<Return>", lambda e: self.add_search_row())

        record = {"frame": row, "field_var": field_var, "combo": field_combo,
                  "op_var": op_var, "text_var": text_var,
                  "text_var2": text_var2, "mid_lbl": mid_lbl, "entry2": entry2}

        def on_op(*_a):
            self._layout_search_row(record)
            self.apply_filter()
        op_combo.bind("<<ComboboxSelected>>", on_op)

        remove_btn = ttk.Button(row, text="✕", width=3,
                                command=lambda: self._remove_search_row(record))
        remove_btn.pack(side="right")
        record["remove_btn"] = remove_btn

        self._layout_search_row(record)
        self.search_rows.append(record)
        self._update_remove_buttons()
        entry1.focus_set()

    def _layout_search_row(self, record):
        """Show/hide the second entry + mid label based on the operator."""
        op = record["op_var"].get()
        record["mid_lbl"].pack_forget()
        record["entry2"].pack_forget()
        if op == self.OP_CONTAINS:
            return
        record["mid_lbl"].configure(text="to" if op == self.OP_RANGE else "±")
        record["mid_lbl"].pack(side="left", padx=4)
        record["entry2"].pack(side="left", padx=(0, 6))

    def _remove_search_row(self, record):
        if len(self.search_rows) <= 1:
            # never remove the last row — just clear it
            record["text_var"].set("")
            record["text_var2"].set("")
            record["op_var"].set(self.OP_CONTAINS)
            record["field_var"].set("(any field)")
            self._layout_search_row(record)
            return
        record["frame"].destroy()
        self.search_rows.remove(record)
        self._update_remove_buttons()
        self.apply_filter()

    def _update_remove_buttons(self):
        # disable the remove button when only one row remains
        only_one = len(self.search_rows) <= 1
        for r in self.search_rows:
            r["remove_btn"].configure(
                state="disabled" if only_one else "normal")

    # -- file handling ------------------------------------------------------
    def open_dialog(self):
        paths = filedialog.askopenfilenames(
            title="Open spectrum file(s)",
            filetypes=[("Spectrum files", "*.msp *.mgf *.mzML *.mzml *.gz"),
                       ("MSP library", "*.msp"),
                       ("Mascot Generic Format", "*.mgf"),
                       ("mzML", "*.mzML *.mzml *.mzML.gz"),
                       ("All files", "*.*")])
        if paths:
            self.load_files(paths)

    def load_files(self, paths):
        """Parse files on a background thread, streaming spectra into the list.

        Spectra appear (and are browsable/searchable) as they are read, so the
        user can start working immediately; a progress bar sits in the status
        bar. MS1 spectra are skipped in the mzML reader (never stored).
        """
        import threading
        if getattr(self, "_loading", False):
            return
        self._loading = True
        self._cancel = threading.Event()
        self._consumed = 0
        self._start_count = len(self.spectra)
        self._load_state = {"buf": [], "file": "", "done": False,
                            "error": None, "skipped": 0}

        t = threading.Thread(target=self._load_worker, args=(list(paths),),
                             daemon=True)
        t.start()
        self._start_progress()
        self.after(80, self._poll_load)

    def _load_worker(self, paths):
        st = self._load_state
        buf = st["buf"]

        def sink(sp):
            if len(sp.mz):                    # ignore peak-less spectra
                buf.append(sp)
            else:
                st["skipped"] += 1
        try:
            for p in paths:
                if self._cancel.is_set():
                    break
                st["file"] = os.path.basename(p)
                parse_file(p, sink=sink, cancel=self._cancel)
        except Exception as exc:  # noqa: BLE001
            st["error"] = str(exc)
        st["done"] = True

    def _poll_load(self):
        st = self._load_state
        buf = st["buf"]
        n = len(buf)
        if n > self._consumed:
            self._ingest(buf[self._consumed:n])
            self._consumed = n
        self.status.config(
            text=f"Loading {st['file']} … {len(self.spectra):,} spectra "
                 f"(you can browse/search now)")
        if not st["done"]:
            self.after(120, self._poll_load)
            return
        # final drain, then finalise
        n = len(buf)
        if n > self._consumed:
            self._ingest(buf[self._consumed:n])
            self._consumed = n
        self._finish_load()

    def _ingest(self, batch):
        """Append a batch of freshly-parsed spectra to the model + list view."""
        conds = self._row_conditions()
        match_all = self.match_mode.get() == "all"
        bounds = self._mz_filter_bounds()
        first_empty = not self.filtered
        labels = []
        for sp in batch:
            self.spectra.append(sp)
            sp.num = len(self.spectra)
            if self._passes(sp, conds, match_all, bounds):
                self.filtered.append(sp)
                labels.append(self._label_for(sp))
        if labels:
            self.listbox.insert("end", *labels)
        self.count_label.config(
            text=f"{len(self.filtered)} / {len(self.spectra)} spectra")
        if first_empty and self.filtered:
            self.listbox.selection_set(0)
            self.listbox.see(0)
            self._show_spectrum(self.filtered[0])

    def _finish_load(self):
        st = self._load_state
        self._loading = False
        self._stop_progress()
        if st["error"]:
            self._alert("Parse error", st["error"])
        self._rebuild_field_list()
        self.apply_filter()          # one clean pass: applies sort + relabels
        added = len(self.spectra) - self._start_count
        if self._cancel.is_set():
            msg = f"Load cancelled — {len(self.spectra):,} spectra loaded."
        else:
            msg = f"Loaded {added:,} spectra ({len(self.spectra):,} total)."
            if st["skipped"]:
                msg += f" Ignored {st['skipped']:,} without peaks."
        self.status.config(text=msg)

    def clear_all(self):
        self.spectra.clear()
        self.filtered.clear()
        self.listbox.delete(0, "end")
        self.meta_tree.delete(*self.meta_tree.get_children())
        self.plot.set_spectrum(None)
        self._rebuild_field_list()
        self.count_label.config(text="0 spectra")
        self.status.config(text="Cleared.")

    def _rebuild_field_list(self):
        fields = set()
        for s in self.spectra:
            fields.update(s.meta.keys())
        sorted_fields = sorted(fields, key=str.lower)
        # search field selectors (one per dynamic search box)
        values = ["(any field)"] + sorted_fields
        for r in self.search_rows:
            r["combo"]["values"] = values
            if r["field_var"].get() not in values:
                r["field_var"].set("(any field)")
        # sort field selector
        sort_values = [self.SORT_FILE_ORDER, self.SORT_NAME,
                       self.SORT_NUM_PEAKS] + sorted_fields
        self.sort_combo["values"] = sort_values
        if self.sort_var.get() not in sort_values:
            self.sort_var.set(self.SORT_FILE_ORDER)
        # list-label field selector
        label_values = [self.LABEL_AUTO] + sorted_fields
        self.label_combo["values"] = label_values
        if self.label_var.get() not in label_values:
            self.label_var.set(self.LABEL_AUTO)

    # -- list labelling -----------------------------------------------------
    def _label_for(self, spec):
        """Build the listbox label for a spectrum, honouring 'Label by'.

        In auto mode use the spectrum's name (with a fallback for files that
        have no name field); otherwise use the chosen field's value, falling
        back to a compact "source #n" tag when that field is empty. The
        precursor m/z is appended unless it is already part of the label.
        """
        field = self.label_var.get()
        if field == self.LABEL_AUTO:
            text = spec.name
            if text == "(unnamed)":
                text = f"{spec.source or 'spectrum'} #{spec.num}"
        else:
            text = str(spec.meta.get(field, "")).strip()
            if not text:
                text = f"{spec.source or 'spectrum'} #{spec.num}"
        prec = spec.precursor
        if prec and str(prec) not in text:
            text = f"{text}   [{prec}]"
        return text

    # -- search / filter ----------------------------------------------------
    def _sort_key(self, spec, field):
        """Sort key that orders numbers numerically and text alphabetically.

        Returns a (kind, number, text) tuple so mixed values stay comparable:
        numeric values sort first (kind 0) by magnitude, non-numeric next
        (kind 1) alphabetically. Missing values sort last.
        """
        if field == self.SORT_NUM_PEAKS:
            return (0, float(len(spec.peaks)), "")
        if field == self.SORT_NAME:
            return (1, 0.0, spec.name.lower())
        raw = spec.meta.get(field, "")
        if raw == "" or raw is None:
            return (2, 0.0, "")          # empties last
        try:
            return (0, float(raw), "")
        except (ValueError, TypeError):
            return (1, 0.0, str(raw).lower())

    def _sort_filtered(self):
        field = self.sort_var.get()
        if field == self.SORT_FILE_ORDER:
            return  # keep original file order
        self.filtered.sort(key=lambda s: self._sort_key(s, field),
                           reverse=self.sort_desc_var.get())

    def _numeric_field_value(self, spec, field):
        """Float value of a field for numeric search (any-field -> precursor)."""
        raw = spec.precursor if field == "(any field)" \
            else spec.meta.get(field, "")
        try:
            return float(str(raw).split()[0])
        except (ValueError, TypeError, IndexError):
            return None

    def _row_conditions(self):
        """Collect active search conditions from the rows."""
        conds = []
        for r in self.search_rows:
            op = r["op_var"].get()
            field = r["field_var"].get()
            t1 = r["text_var"].get().strip()
            t2 = r["text_var2"].get().strip()
            if op == self.OP_CONTAINS:
                if t1:
                    conds.append((op, field, t1.lower(), None))
            elif op == self.OP_RANGE:
                lo = _to_float(t1)
                hi = _to_float(t2)
                if lo is not None or hi is not None:
                    conds.append((op, field, lo, hi))
            elif op == self.OP_TOL:
                val = _to_float(t1)
                tol = _to_float(t2) or 0.0
                if val is not None:
                    conds.append((op, field, val, tol))
        return conds

    def _match_condition(self, spec, cond):
        op, field, a, b = cond
        if op == self.OP_CONTAINS:
            hay = spec.search_blob() if field == "(any field)" \
                else str(spec.meta.get(field, "")).lower()
            return a in hay
        v = self._numeric_field_value(spec, field)
        if v is None:
            return False
        if op == self.OP_RANGE:
            if a is not None and v < a:
                return False
            if b is not None and v > b:
                return False
            return True
        # value ± tol
        return abs(v - a) <= b

    def schedule_filter(self):
        """Debounce rapid typing so filtering only runs after a short pause."""
        after = getattr(self, "_filter_after", None)
        if after:
            try:
                self.after_cancel(after)
            except Exception:  # noqa: BLE001
                pass
        self._filter_after = self.after(250, self.apply_filter)

    def _passes(self, spec, conditions, match_all, bounds):
        """True if a spectrum satisfies the text conditions + m/z filter."""
        if conditions:
            if match_all:
                if not all(self._match_condition(spec, c)
                           for c in conditions):
                    return False
            else:
                if not any(self._match_condition(spec, c)
                           for c in conditions):
                    return False
        if bounds:
            low, high = bounds
            if self.mz_target_kind.get() == "precursor":
                try:
                    p = float(spec.precursor)
                except (ValueError, TypeError):
                    return False
                if not (low <= p <= high):
                    return False
            else:  # any peak (iterate the compact m/z array directly)
                if not any(low <= mz <= high for mz in spec.mz):
                    return False
        return True

    def apply_filter(self):
        self._filter_after = None
        conditions = self._row_conditions()
        match_all = self.match_mode.get() == "all"
        bounds = self._mz_filter_bounds()
        if not conditions and not bounds:
            self.filtered = list(self.spectra)
        else:
            self.filtered = [s for s in self.spectra
                             if self._passes(s, conditions, match_all, bounds)]

        self._sort_filtered()

        # bulk-insert labels in one Tcl call (fast even at 100k+ rows)
        self.listbox.delete(0, "end")
        if self.filtered:
            self.listbox.insert("end", *[self._label_for(s)
                                         for s in self.filtered])

        self.count_label.config(
            text=f"{len(self.filtered)} / {len(self.spectra)} spectra")
        if self.filtered:
            self.listbox.selection_clear(0, "end")
            self.listbox.selection_set(0)
            self.listbox.see(0)
            self._show_spectrum(self.filtered[0])
        else:
            self.meta_tree.delete(*self.meta_tree.get_children())
            self.plot.set_spectrum(None)

    # -- m/z filter ---------------------------------------------------------
    def _compute_mz(self):
        formula = self.mz_formula.get().strip()
        if not formula:
            return
        extra = _to_float(self.mz_extra.get()) or 0.0
        try:
            neutral = formula_mass(formula) + extra
            val = adduct_mz(neutral, self.mz_adduct.get())
        except (ValueError, KeyError) as exc:
            self.status.config(text=f"Formula error: {exc}")
            return
        self.mz_enable.set(True)
        self.mz_target.set(f"{val:.4f}")   # triggers apply_filter via trace
        extra_txt = f" +Δ{extra:g}" if extra else ""
        self.status.config(
            text=f"{formula}{extra_txt} {self.mz_adduct.get()} "
                 f"→ m/z {val:.4f} (monoisotopic)")

    def _on_mz_filter_change(self):
        self.schedule_filter()

    def _on_mz_enable(self):
        # show the options only while the filter is enabled
        if self.mz_enable.get():
            self.mz_options.pack(side="left")
        else:
            self.mz_options.pack_forget()
        self.apply_filter()

    def _mz_filter_bounds(self):
        """Return (low, high) m/z window if the filter is active, else None."""
        if not self.mz_enable.get():
            return None
        try:
            target = float(self.mz_target.get())
        except (ValueError, TypeError):
            return None
        try:
            tol = float(self.mz_tol.get())
        except (ValueError, TypeError):
            tol = 0.0
        if self.mz_tol_unit.get() == "ppm":
            window = target * tol / 1e6
        else:
            window = tol
        return (target - window, target + window)

    # -- selection ----------------------------------------------------------
    def _selected_spectra(self):
        """Spectra highlighted in the list (falls back to the displayed one)."""
        sel = self.listbox.curselection()
        if sel:
            return [self.filtered[i] for i in sel
                    if 0 <= i < len(self.filtered)]
        return self.filtered[:1]

    def _on_select(self, _event):
        sel = self.listbox.curselection()
        if not sel:
            return
        idx = sel[0]
        if 0 <= idx < len(self.filtered):
            self._show_spectrum(self.filtered[idx])
        if len(sel) > 1:
            self.status.config(
                text=f"{len(sel)} spectra selected — "
                     f"File ▸ Save to export.")

    def _show_spectrum(self, spec):
        self.plot.set_relative(self.relative_var.get())
        self.plot.set_spectrum(spec)
        self._current_spec = spec
        self._refresh_meta(spec)
        self.status.config(
            text=f"{spec.name} — {len(spec.peaks)} peaks — {spec.source}")

    def _refresh_meta(self, spec):
        self.meta_tree.delete(*self.meta_tree.get_children())
        rows = list(spec.meta.items())
        rows.append(("# peaks", len(spec.peaks)))
        rows.append(("source file", spec.source))
        for k, v in rows:
            if k not in self.hidden_fields:
                self.meta_tree.insert("", "end", values=(k, v))

    def _all_meta_fields(self):
        """Every metadata field name across loaded spectra (+ synthetic ones)."""
        fields = []
        seen = set()
        for s in self.spectra:
            for k in s.meta:
                if k not in seen:
                    seen.add(k)
                    fields.append(k)
        for k in ("# peaks", "source file"):
            if k not in seen:
                fields.append(k)
        return fields

    def choose_meta_fields(self):
        """Dialog with a checkbox per metadata field to show/hide it."""
        fields = self._all_meta_fields()
        if not fields:
            self._alert("Metadata fields", "Load a spectrum first.")
            return
        dlg = tk.Toplevel(self)
        dlg.title("Metadata fields to show")
        dlg.transient(self)
        dlg.geometry("320x420")
        ttk.Label(dlg, text="Check the fields to display:").pack(
            anchor="w", padx=10, pady=(10, 4))

        canvas = tk.Canvas(dlg, highlightthickness=0)
        sb = ttk.Scrollbar(dlg, orient="vertical", command=canvas.yview)
        inner = ttk.Frame(canvas)
        inner.bind("<Configure>",
                   lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=inner, anchor="nw")
        canvas.configure(yscrollcommand=sb.set)
        canvas.pack(side="top", fill="both", expand=True, padx=(10, 0))
        sb.pack(side="right", fill="y")

        vars_ = {}
        for f in fields:
            v = tk.BooleanVar(value=f not in self.hidden_fields)
            vars_[f] = v
            ttk.Checkbutton(inner, text=f, variable=v).pack(anchor="w")

        btns = ttk.Frame(dlg)
        btns.pack(side="bottom", fill="x", padx=10, pady=8)

        def _set_all(state):
            for v in vars_.values():
                v.set(state)

        def _apply():
            self.hidden_fields = {f for f, v in vars_.items() if not v.get()}
            if getattr(self, "_current_spec", None):
                self._refresh_meta(self._current_spec)
            dlg.destroy()

        ttk.Button(btns, text="All", command=lambda: _set_all(True)).pack(
            side="left")
        ttk.Button(btns, text="None", command=lambda: _set_all(False)).pack(
            side="left", padx=(4, 0))
        ttk.Button(btns, text="Apply", command=_apply).pack(side="right")
        ttk.Button(btns, text="Cancel", command=dlg.destroy).pack(
            side="right", padx=(0, 4))

    def _toggle_relative(self):
        self.plot.set_relative(self.relative_var.get())

    # -- saving / export ----------------------------------------------------
    def save_spectra(self, fmt, specs=None):
        if specs is None:
            specs = self._selected_spectra()
        if not specs:
            self._alert("Nothing to save",
                        "Load and select at least one spectrum first.")
            return
        ext = "." + fmt
        path = filedialog.asksaveasfilename(
            defaultextension=ext, title=f"Save {len(specs)} spectra as {ext}",
            filetypes=[(fmt.upper(), "*" + ext), ("All files", "*.*")])
        if not path:
            return
        writer = spectrum_to_msp if fmt == "msp" else spectrum_to_mgf
        try:
            with open(path, "w", encoding="utf-8") as fh:
                fh.write("\n".join(writer(s) for s in specs))
        except Exception as exc:  # noqa: BLE001
            self._alert("Save error", str(exc))
            return
        self.status.config(
            text=f"Saved {len(specs)} spectra to {os.path.basename(path)}")

    def _export_params(self):
        """Render kwargs shared by the export functions and the size preview."""
        return dict(
            relative=self.relative_var.get(), scheme=self.scheme_name.get(),
            bar_width=self.bar_width.get(), transparency=self.bar_trans.get(),
            bg=self.bg_scheme.get(), limit_prec=self.limit_prec_var.get(),
            annotations=[dict(a) for a in self.plot.annotations],
            style=dict(show_xticks=self.plot.show_xticks,
                       show_yticks=self.plot.show_yticks,
                       show_titles=self.plot.show_titles,
                       show_peaklabels=self.plot.show_peaklabels,
                       peaklabel_auto=self.plot.peaklabel_auto,
                       peaklabel_n=self.plot.peaklabel_n,
                       xtitle=self.plot.xtitle, ytitle=self.plot.ytitle,
                       axis_fontsize=self.plot.axis_fontsize,
                       xtitle_fontsize=self.plot.xtitle_fontsize,
                       ytitle_fontsize=self.plot.ytitle_fontsize,
                       peaklabel_fontsize=self.plot.peaklabel_fontsize,
                       tick_bold=self.plot.tick_bold,
                       tick_italic=self.plot.tick_italic,
                       xtitle_bold=self.plot.xtitle_bold,
                       xtitle_italic=self.plot.xtitle_italic,
                       ytitle_bold=self.plot.ytitle_bold,
                       ytitle_italic=self.plot.ytitle_italic,
                       peaklabel_bold=self.plot.peaklabel_bold,
                       peaklabel_italic=self.plot.peaklabel_italic,
                       tick_color=self.plot.tick_color,
                       xtitle_color=self.plot.xtitle_color,
                       ytitle_color=self.plot.ytitle_color,
                       peaklabel_color=self.plot.peaklabel_color))

    def save_image(self, fmt, specs=None):
        if specs is None:
            specs = self._selected_spectra()
        if not specs:
            self._alert("Nothing to save",
                        "Load and select at least one spectrum first.")
            return
        size = self._ask_image_size(specs)
        if size is None:
            return
        w, h = size
        ext = "." + fmt
        path = filedialog.asksaveasfilename(
            defaultextension=ext, title=f"Save plot as {ext}",
            filetypes=[(fmt.upper(), "*" + ext), ("All files", "*.*")])
        if not path:
            return
        pr = self._export_params()
        try:
            if fmt == "png":
                if len(specs) == 1:
                    render_spectrum_image(specs[0], width=w, height=h,
                                          **pr).save(path)
                    n = 1
                else:
                    stem, _ = os.path.splitext(path)
                    for i, s in enumerate(specs, 1):
                        render_spectrum_image(s, width=w, height=h,
                                              **pr).save(f"{stem}_{i}.png")
                    n = len(specs)
            else:  # pdf — one vector page per spectrum
                save_spectra_pdf(specs, path, width=w, height=h, **pr)
                n = len(specs)
        except Exception as exc:  # noqa: BLE001
            self._alert("Export error", str(exc))
            return
        files = f"{n} files" if (fmt == "png" and n > 1) else "1 file"
        self.status.config(
            text=f"Exported {n} {fmt.upper()} ({files}).")


def main():
    files = [a for a in sys.argv[1:] if os.path.isfile(a)]
    app = App(initial_files=files or None)
    app.mainloop()


if __name__ == "__main__":
    main()
