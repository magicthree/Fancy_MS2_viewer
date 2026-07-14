#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MSP Viewer - a simple viewer for .msp MS/MS spectral library files.

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
        self.on_context = None   # callback(event) for right-click menu
        self.view_min = None     # current x (m/z) view window
        self.view_max = None
        self._full_min = None
        self._full_max = None
        self._drag_start = None
        self._drag_rect = None
        self._tip = None

        self.bind("<Configure>", lambda e: self.redraw())
        self.bind("<Motion>", self._on_motion)
        self.bind("<Leave>", lambda e: self._hide_tip())
        self.bind("<ButtonPress-1>", self._on_press)
        self.bind("<B1-Motion>", self._on_drag)
        self.bind("<ButtonRelease-1>", self._on_release)
        self.bind("<Button-3>", self._on_right_click)
        self.bind("<MouseWheel>", self._on_wheel)          # Windows / macOS
        self.bind("<Button-4>", lambda e: self._zoom(0.8, e.x))  # Linux
        self.bind("<Button-5>", lambda e: self._zoom(1.25, e.x))

    # -- public API ---------------------------------------------------------
    def set_spectrum(self, spec):
        self.peaks = list(spec.peaks) if spec else []
        self.title = spec.name if spec else ""
        if self.peaks:
            mzs = [p[0] for p in self.peaks]
            lo, hi = min(mzs), max(mzs)
            span = max(hi - lo, 1.0)
            self._full_min = lo - span * 0.03
            self._full_max = hi + span * 0.03
        else:
            self._full_min, self._full_max = 0.0, 100.0
        self.view_min, self.view_max = self._full_min, self._full_max
        self.redraw()

    def set_relative(self, flag):
        self.relative = bool(flag)
        self.redraw()

    def set_scheme(self, scheme):
        self.scheme = scheme
        self.redraw()

    def set_bar_style(self, width, transparency):
        self.bar_width = max(1, int(width))
        self.transparency = max(0, min(100, int(transparency)))
        self.redraw()

    def reset_view(self):
        self.view_min, self.view_max = self._full_min, self._full_max
        self.redraw()

    def _on_right_click(self, event):
        if self.on_context is not None:
            self.on_context(event)
        else:
            self.reset_view()

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

    # -- drawing ------------------------------------------------------------
    def redraw(self):
        self.delete("all")
        w = self.winfo_width()
        h = self.winfo_height()
        if w < 120 or h < 90:
            return
        x0, y0, x1, y1 = self._plot_area()

        # axes
        self.create_line(x0, y1, x1, y1, fill="#444")   # x axis
        self.create_line(x0, y0, x0, y1, fill="#444")   # y axis

        if not self.peaks:
            self.create_text((x0 + x1) / 2, (y0 + y1) / 2,
                             text="No peaks to display",
                             fill="#999", font=("Arial", 11))
            return

        # peaks within current view
        vis = [(mz, it) for (mz, it) in self.peaks
               if self.view_min <= mz <= self.view_max]
        if vis:
            base = max(it for _, it in vis)
        else:
            base = max((it for _, it in self.peaks), default=1.0)
        if base <= 0:
            base = 1.0

        # title (Microsoft YaHei UI when it contains Chinese, else Arial)
        title_family = UI_FONT if _has_cjk(self.title) else "Arial"
        self.create_text(x0, 14, anchor="w", text=self.title,
                         fill="#222", font=(title_family, 10, "bold"))

        # y grid + labels
        self._draw_y_axis(x0, y0, x1, y1, base)
        # x ticks
        self._draw_x_axis(x0, y0, x1, y1)

        # sticks (coloured by scheme; thickness + transparency configurable)
        for mz, it in vis:
            x = self._mz_to_x(mz, x0, x1)
            frac = it / base
            y = y1 - frac * (y1 - y0)
            color = _blend_white(scheme_color(self.scheme, frac),
                                 self.transparency)
            self.create_line(x, y1, x, y, fill=_hex(color),
                             width=self.bar_width)

        # label the tallest peaks
        vis_sorted = sorted(vis, key=lambda p: p[1], reverse=True)[:8]
        for mz, it in vis_sorted:
            x = self._mz_to_x(mz, x0, x1)
            frac = it / base
            y = y1 - frac * (y1 - y0)
            self.create_text(x, y - 6, text=f"{mz:.4f}".rstrip("0").rstrip("."),
                             anchor="s", fill="#c0392b",
                             font=("Arial", 8))

    def _draw_y_axis(self, x0, y0, x1, y1, base):
        for i in range(0, 6):
            frac = i / 5.0
            y = y1 - frac * (y1 - y0)
            self.create_line(x0 - 4, y, x0, y, fill="#444")
            if self.relative:
                label = f"{frac * 100:.0f}%"
            else:
                label = _fmt_si(frac * base)
            self.create_text(x0 - 8, y, anchor="e", text=label,
                             fill="#555", font=("Arial", 8))
        # axis title
        self.create_text(16, (y0 + y1) / 2, angle=90,
                         text="Relative intensity" if self.relative
                         else "Intensity",
                         fill="#555", font=("Arial", 8))

    def _draw_x_axis(self, x0, y0, x1, y1):
        lo, hi = self.view_min, self.view_max
        span = hi - lo
        if span <= 0:
            return
        step = _nice_step(span / 8.0)
        start = (int(lo / step) + 1) * step
        v = start
        while v < hi:
            x = self._mz_to_x(v, x0, x1)
            self.create_line(x, y1, x, y1 + 4, fill="#444")
            self.create_text(x, y1 + 7, anchor="n",
                             text=f"{v:g}", fill="#555",
                             font=("Arial", 8))
            v += step
        self.create_text((x0 + x1) / 2, y1 + 30, text="m/z",
                         fill="#555", font=("Arial", 9))

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
        self._drag_start = e.x

    def _on_drag(self, e):
        if self._drag_start is None:
            return
        if self._drag_rect:
            self.delete(self._drag_rect)
        x0, y0, x1, y1 = self._plot_area()
        self._drag_rect = self.create_rectangle(
            self._drag_start, y0, e.x, y1,
            outline="#3498db", fill="#3498db", stipple="gray25")

    def _on_release(self, e):
        if self._drag_rect:
            self.delete(self._drag_rect)
            self._drag_rect = None
        if self._drag_start is None:
            return
        x0, y0, x1, y1 = self._plot_area()
        a, b = self._drag_start, e.x
        self._drag_start = None
        if abs(b - a) < 6:      # treat as click, not a zoom
            return
        lo = self._x_to_mz(min(a, b), x0, x1)
        hi = self._x_to_mz(max(a, b), x0, x1)
        if hi - lo > 1e-6:
            self.view_min, self.view_max = lo, hi
            self.redraw()

    def _on_wheel(self, e):
        factor = 0.8 if e.delta > 0 else 1.25
        self._zoom(factor, e.x)

    def _zoom(self, factor, px):
        x0, y0, x1, y1 = self._plot_area()
        center = self._x_to_mz(px, x0, x1)
        lo = center - (center - self.view_min) * factor
        hi = center + (self.view_max - center) * factor
        # clamp to full range
        lo = max(lo, self._full_min)
        hi = min(hi, self._full_max)
        if hi - lo > 1e-6:
            self.view_min, self.view_max = lo, hi
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
ADDUCTS = {
    "[M+H]+": (1, "H", "", 1),
    "[M+Na]+": (1, "Na", "", 1),
    "[M+K]+": (1, "K", "", 1),
    "[M+NH4]+": (1, "NH4", "", 1),
    "[M+H-H2O]+": (1, "H", "H2O", 1),
    "[M+2H]2+": (1, "H2", "", 2),
    "[M+3H]3+": (1, "H3", "", 3),
    "[2M+H]+": (2, "H", "", 1),
    "[M]+": (1, "", "", 1),
    "[M-H]-": (1, "", "H", -1),
    "[M+Cl]-": (1, "Cl", "", -1),
    "[M+HCOO]-": (1, "CHO2", "", -1),
    "[M+CH3COO]-": (1, "C2H3O2", "", -1),
    "[M-2H]2-": (1, "", "H2", -2),
    "[2M-H]-": (2, "", "H", -1),
    "[M]-": (1, "", "", -1),
}

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


def adduct_mz(neutral_mass, adduct):
    """m/z for a neutral monoisotopic mass under the given adduct."""
    mult, add_f, rem_f, z = ADDUCTS[adduct]
    m = neutral_mass * mult
    if add_f:
        m += formula_mass(add_f)
    if rem_f:
        m -= formula_mass(rem_f)
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
UI_FONT = "Microsoft YaHei UI" if _find_font("msyh.ttc") else "Arial"


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
        self.img = Image.new("RGB", (width, height), "white")
        self.dr = ImageDraw.Draw(self.img)
        self._cache = {}

    def _f(self, size, cjk):
        key = (size, cjk)
        if key in self._cache:
            return self._cache[key]
        ImageFont = self._ImageFont
        fnt = None
        path = YAHEI_PATH if cjk else ARIAL_PATH
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

    def text(self, x, y, s, rgb, size=15, anchor="lt", rotate=0):
        fnt = self._f(size, _has_cjk(s))
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

    def text(self, x, y, s, rgb, size=15, anchor="lt", rotate=0):
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


_RL_STATE = {"ready": False, "arial": False, "yahei": False}


def _ensure_rl_fonts():
    """Register Arial + Microsoft YaHei with reportlab once (embeds subsets)."""
    if _RL_STATE["ready"]:
        return
    _RL_STATE["ready"] = True
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    if ARIAL_PATH:
        try:
            pdfmetrics.registerFont(TTFont("Arial", ARIAL_PATH))
            _RL_STATE["arial"] = True
        except Exception:  # noqa: BLE001
            pass
    if YAHEI_PATH:
        try:
            if YAHEI_PATH.lower().endswith(".ttc"):
                pdfmetrics.registerFont(TTFont("YaHei", YAHEI_PATH,
                                               subfontIndex=0))
            else:
                pdfmetrics.registerFont(TTFont("YaHei", YAHEI_PATH))
            _RL_STATE["yahei"] = True
        except Exception:  # noqa: BLE001
            pass


class _RLBackend:
    """Vector PDF backend via reportlab, embedding Arial + YaHei subsets."""

    def __init__(self, path, width, height):
        from reportlab.pdfgen import canvas
        _ensure_rl_fonts()
        self.c = canvas.Canvas(path, pagesize=(width, height))
        self.h = height

    def new_page(self):
        self.c.showPage()

    def _font_for(self, s):
        if _has_cjk(s) and _RL_STATE["yahei"]:
            return "YaHei"
        return "Arial" if _RL_STATE["arial"] else "Helvetica"

    def line(self, x0, y0, x1, y1, rgb, width=1):
        self.c.setStrokeColorRGB(*(v / 255.0 for v in rgb))
        self.c.setLineWidth(width)
        Y = self.h
        self.c.line(x0, Y - y0, x1, Y - y1)

    def text(self, x, y, s, rgb, size=15, anchor="lt", rotate=0):
        self.c.setFillColorRGB(*(v / 255.0 for v in rgb))
        font = self._font_for(s)
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


def _draw_spectrum(be, width, height, spec, relative=True,
                   scheme=DEFAULT_SCHEME, bar_width=2, transparency=0):
    """Draw a spectrum stick-plot onto a backend (shared by PNG and PDF)."""
    axis, grey, label = (68, 68, 68), (85, 85, 85), (192, 60, 43)
    pad_l, pad_r, pad_t, pad_b = 90, 30, 55, 70
    x0, y0 = pad_l, pad_t
    x1, y1 = width - pad_r, height - pad_b

    title = spec.name if spec.name != "(unnamed)" else \
        f"{spec.source} #{spec.num}"
    prec = spec.precursor
    if prec:
        title += f"    precursor m/z {prec}"
    be.text(x0, 20, title, (17, 17, 17), size=22, anchor="lm")

    be.line(x0, y1, x1, y1, axis, 1)
    be.line(x0, y0, x0, y1, axis, 1)

    peaks = spec.peaks or []
    if not peaks:
        be.text((x0 + x1) / 2, (y0 + y1) / 2, "No peaks", (150, 150, 150),
                size=15, anchor="cm")
        return

    mzs = [p[0] for p in peaks]
    lo, hi = min(mzs), max(mzs)
    span = max(hi - lo, 1.0)
    vmin, vmax = lo - span * 0.03, hi + span * 0.03
    base = max((it for _, it in peaks), default=1.0) or 1.0

    def mz2x(mz):
        return x0 + (mz - vmin) / (vmax - vmin) * (x1 - x0)

    for i in range(6):
        frac = i / 5.0
        y = y1 - frac * (y1 - y0)
        be.line(x0 - 5, y, x0, y, axis, 1)
        lab = f"{frac * 100:.0f}%" if relative else _fmt_si(frac * base)
        be.text(x0 - 10, y, lab, grey, size=14, anchor="rm")

    step = _nice_step((vmax - vmin) / 8.0)
    v = (int(vmin / step) + 1) * step
    while v < vmax:
        x = mz2x(v)
        be.line(x, y1, x, y1 + 5, axis, 1)
        be.text(x, y1 + 8, f"{v:g}", grey, size=14, anchor="ct")
        v += step
    be.text((x0 + x1) / 2, y1 + 42, "m/z", grey, size=15, anchor="ct")
    be.text(24, (y0 + y1) / 2,
            "Relative intensity" if relative else "Intensity",
            grey, size=14, anchor="cm", rotate=90)

    for mz, it in peaks:
        x = mz2x(mz)
        frac = it / base
        y = y1 - frac * (y1 - y0)
        be.line(x, y1, x, y,
                _blend_white(scheme_color(scheme, frac), transparency),
                bar_width)

    for mz, it in sorted(peaks, key=lambda p: p[1], reverse=True)[:10]:
        x = mz2x(mz)
        y = y1 - (it / base) * (y1 - y0)
        be.text(x, y - 6, f"{mz:.4f}".rstrip("0").rstrip("."),
                label, size=14, anchor="cb")


def render_spectrum_image(spec, relative=True, width=1200, height=700,
                          scheme=DEFAULT_SCHEME, bar_width=2, transparency=0):
    """Render a spectrum stick-plot to a PIL RGB image (raster)."""
    be = _PILBackend(width, height)
    _draw_spectrum(be, width, height, spec, relative, scheme,
                   bar_width, transparency)
    return be.img


def save_spectra_pdf(specs, path, relative=True, scheme=DEFAULT_SCHEME,
                     width=1200, height=700, bar_width=2, transparency=0):
    """Write one or more spectra to a multi-page vector PDF.

    Text is embedded as Arial (Latin) / Microsoft YaHei (CJK) subsets.
    """
    be = _make_pdf_backend(path, width, height)
    for i, s in enumerate(specs):
        if i > 0:
            be.new_page()
        _draw_spectrum(be, width, height, s, relative, scheme,
                       bar_width, transparency)
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

        self._build_ui()
        self._build_menu()
        self._setup_dnd()

        self.plot.set_scheme(self.scheme_name.get())
        self.plot.set_bar_style(self.bar_width.get(), self.bar_trans.get())
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
        setmenu.add_separator()
        setmenu.add_checkbutton(
            label="Show relative intensity (%)  (else absolute)",
            variable=self.relative_var, command=self._toggle_relative)
        setmenu.add_separator()
        setmenu.add_command(label="Metadata fields to show…",
                            command=self.choose_meta_fields)
        menubar.add_cascade(label="Settings", menu=setmenu)

        self.config(menu=menubar)
        self.bind("<Control-o>", lambda e: self.open_dialog())
        self.bind("<Control-s>", lambda e: self.save_spectra("msp"))
        self.bind("<Control-c>", lambda e: self.copy_peaks())

    def _on_scheme_change(self):
        self.plot.set_scheme(self.scheme_name.get())

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

    def _ask_image_size(self):
        """Ask for export image size: current preview, or a custom W×H.

        Returns (width, height) or None if cancelled.
        """
        pw = max(int(self.plot.winfo_width()), 100)
        ph = max(int(self.plot.winfo_height()), 100)
        dlg = tk.Toplevel(self)
        dlg.title("Image size")
        dlg.transient(self)
        dlg.resizable(False, False)
        result = {"size": None}
        mode = tk.StringVar(value="preview")
        wv = tk.StringVar(value="1200")
        hv = tk.StringVar(value="700")

        ttk.Radiobutton(dlg, text=f"Current preview size ({pw} × {ph})",
                        value="preview", variable=mode).grid(
            row=0, column=0, columnspan=3, sticky="w", padx=10, pady=(12, 2))
        ttk.Radiobutton(dlg, text="Custom size:", value="custom",
                        variable=mode).grid(row=1, column=0, sticky="w",
                                            padx=10)
        ttk.Entry(dlg, textvariable=wv, width=7).grid(row=1, column=1)
        ttk.Label(dlg, text="×").grid(row=1, column=2)
        ttk.Entry(dlg, textvariable=hv, width=7).grid(row=1, column=3,
                                                      padx=(0, 10))

        btns = ttk.Frame(dlg)
        btns.grid(row=2, column=0, columnspan=4, sticky="e", padx=10,
                  pady=(8, 12))

        def ok():
            if mode.get() == "preview":
                result["size"] = (pw, ph)
            else:
                w = _to_float(wv.get())
                h = _to_float(hv.get())
                if not w or not h or w < 50 or h < 50:
                    self._alert("Image size", "Enter valid width and height "
                                              "(≥ 50 px).")
                    return
                result["size"] = (int(w), int(h))
            dlg.destroy()
        ttk.Button(btns, text="OK", command=ok).pack(side="right")
        ttk.Button(btns, text="Cancel", command=dlg.destroy).pack(
            side="right", padx=(0, 6))
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
        self.mz_enable.trace_add("write", lambda *a: self.apply_filter())

        # target m/z (typed directly, or computed from formula+adduct)
        self.mz_target = tk.StringVar()
        self.mz_target.trace_add("write", lambda *a: self._on_mz_filter_change())
        ttk.Label(row3, text="m/z").pack(side="left", padx=(8, 2))
        ttk.Entry(row3, textvariable=self.mz_target, width=12).pack(side="left")

        ttk.Label(row3, text="±").pack(side="left", padx=(6, 2))
        self.mz_tol = tk.StringVar(value="10")
        self.mz_tol.trace_add("write", lambda *a: self._on_mz_filter_change())
        ttk.Entry(row3, textvariable=self.mz_tol, width=7).pack(side="left")
        self.mz_tol_unit = tk.StringVar(value="ppm")
        ttk.Combobox(row3, textvariable=self.mz_tol_unit, width=5,
                     state="readonly", values=["ppm", "Da"]).pack(
            side="left", padx=(2, 8))
        self.mz_tol_unit.trace_add("write",
                                   lambda *a: self._on_mz_filter_change())

        ttk.Label(row3, text="match").pack(side="left")
        self.mz_target_kind = tk.StringVar(value="precursor")
        ttk.Combobox(row3, textvariable=self.mz_target_kind, width=10,
                     state="readonly",
                     values=["precursor", "any peak"]).pack(side="left",
                                                            padx=(4, 12))
        self.mz_target_kind.trace_add("write",
                                      lambda *a: self._on_mz_filter_change())

        # formula (+ extra mass) + adduct -> computes the m/z above
        ttk.Label(row3, text="Formula").pack(side="left")
        self.mz_formula = tk.StringVar()
        self.mz_formula.trace_add("write", lambda *a: self._compute_mz())
        ttk.Entry(row3, textvariable=self.mz_formula, width=14).pack(
            side="left", padx=(4, 4))
        ttk.Label(row3, text="+Δmass").pack(side="left")
        self.mz_extra = tk.StringVar()
        self.mz_extra.trace_add("write", lambda *a: self._compute_mz())
        ttk.Entry(row3, textvariable=self.mz_extra, width=9).pack(
            side="left", padx=(4, 8))
        ttk.Label(row3, text="Adduct").pack(side="left")
        self.mz_adduct = tk.StringVar(value="[M+H]+")
        ttk.Combobox(row3, textvariable=self.mz_adduct, width=13,
                     state="readonly", values=list(ADDUCTS.keys())).pack(
            side="left", padx=(4, 0))
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
                         text="Drag = zoom range   •   Wheel = zoom   •   "
                              "Right-click = menu   •   Hover = peak info",
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

    def save_image(self, fmt, specs=None):
        if specs is None:
            specs = self._selected_spectra()
        if not specs:
            self._alert("Nothing to save",
                        "Load and select at least one spectrum first.")
            return
        size = self._ask_image_size()
        if size is None:
            return
        w, h = size
        ext = "." + fmt
        path = filedialog.asksaveasfilename(
            defaultextension=ext, title=f"Save plot as {ext}",
            filetypes=[(fmt.upper(), "*" + ext), ("All files", "*.*")])
        if not path:
            return
        rel = self.relative_var.get()
        scheme = self.scheme_name.get()
        bw = self.bar_width.get()
        tr = self.bar_trans.get()
        try:
            if fmt == "png":
                if len(specs) == 1:
                    render_spectrum_image(specs[0], rel, w, h, scheme,
                                          bw, tr).save(path)
                    n = 1
                else:
                    stem, _ = os.path.splitext(path)
                    for i, s in enumerate(specs, 1):
                        render_spectrum_image(s, rel, w, h, scheme,
                                              bw, tr).save(f"{stem}_{i}.png")
                    n = len(specs)
            else:  # pdf — one vector page per spectrum
                save_spectra_pdf(specs, path, relative=rel, scheme=scheme,
                                 width=w, height=h, bar_width=bw,
                                 transparency=tr)
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
