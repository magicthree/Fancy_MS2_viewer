"""Screenshot test on the real timsTOF mzML: m/z filter + colour/bar settings."""
import os
import time
import ctypes
from ctypes import wintypes

from PIL import Image
import fancyMS2viewer as m

FILE = r"E:\81_Methanogen\mzML\2016_1-MT.F-1 injekt_20% _GA1.d.mzML"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(OUT, exist_ok=True)

user32 = ctypes.windll.user32
gdi32 = ctypes.windll.gdi32
try:
    ctypes.windll.shcore.SetProcessDpiAwareness(2)
except Exception:
    pass


def capture(widget, path):
    widget.update_idletasks(); widget.update(); time.sleep(0.35); widget.update()
    hwnd = user32.GetAncestor(widget.winfo_id(), 2)
    rect = wintypes.RECT()
    user32.GetWindowRect(hwnd, ctypes.byref(rect))
    w, h = rect.right - rect.left, rect.bottom - rect.top
    hdc = user32.GetWindowDC(hwnd)
    memdc = gdi32.CreateCompatibleDC(hdc)
    bmp = gdi32.CreateCompatibleBitmap(hdc, w, h)
    gdi32.SelectObject(memdc, bmp)
    user32.PrintWindow(hwnd, memdc, 2)

    class BMIH(ctypes.Structure):
        _fields_ = [("biSize", wintypes.DWORD), ("biWidth", wintypes.LONG),
                    ("biHeight", wintypes.LONG), ("biPlanes", wintypes.WORD),
                    ("biBitCount", wintypes.WORD),
                    ("biCompression", wintypes.DWORD),
                    ("biSizeImage", wintypes.DWORD),
                    ("biXPelsPerMeter", wintypes.LONG),
                    ("biYPelsPerMeter", wintypes.LONG),
                    ("biClrUsed", wintypes.DWORD),
                    ("biClrImportant", wintypes.DWORD)]

    class BMI(ctypes.Structure):
        _fields_ = [("bmiHeader", BMIH), ("bmiColors", wintypes.DWORD * 3)]

    bmi = BMI()
    bmi.bmiHeader.biSize = ctypes.sizeof(BMIH)
    bmi.bmiHeader.biWidth = w
    bmi.bmiHeader.biHeight = -h
    bmi.bmiHeader.biPlanes = 1
    bmi.bmiHeader.biBitCount = 32
    buf = ctypes.create_string_buffer(w * h * 4)
    gdi32.GetDIBits(memdc, bmp, 0, h, buf, ctypes.byref(bmi), 0)
    img = Image.frombuffer("RGB", (w, h), buf, "raw", "BGRX", 0, 1)
    gdi32.DeleteObject(bmp); gdi32.DeleteDC(memdc); user32.ReleaseDC(hwnd, hdc)
    img.save(path)
    print("saved", os.path.basename(path), img.size)


app = m.App()
app.geometry("1180x720+80+40")
app.update()
t0 = time.time()
app.load_files([FILE])
while getattr(app, "_loading", False):
    app.update(); time.sleep(0.02)
print(f"loaded {len(app.spectra):,} MS2 spectra in {time.time()-t0:.1f}s")

# --- m/z filter: 653.681, 10 ppm ---
app.mz_target.set("653.681")
app.mz_tol.set("10")
app.mz_tol_unit.set("ppm")
app.mz_target_kind.set("precursor")
app.mz_enable.set(True)
app.apply_filter()
prec_n = app.listbox.size()
kind = "precursor"
if prec_n == 0:
    app.mz_target_kind.set("any peak")
    app.apply_filter()
    kind = "any peak"
n = app.listbox.size()
print(f"precursor matches: {prec_n} | using kind='{kind}' -> {n} spectra")

if n:
    app.listbox.selection_clear(0, "end")
    app.listbox.selection_set(0)
    app._show_spectrum(app.filtered[0])
app.lift(); app.attributes("-topmost", True)
capture(app, os.path.join(OUT, "01_search_653p681_10ppm.png"))
app.attributes("-topmost", False)

# --- colour + bar settings dialogs ---
import tkinter as tk
_orig = tk.Toplevel.wait_window
tk.Toplevel.wait_window = lambda self, *a: None

app.edit_custom_color()
app.update()
dlg = [w for w in app.winfo_children() if isinstance(w, tk.Toplevel)][-1]
dlg.lift(); dlg.attributes("-topmost", True)
capture(dlg, os.path.join(OUT, "02_peak_colour_gradient.png"))

app._pick_color((70, 190, 255))
picker = [w for w in app.winfo_children() if isinstance(w, tk.Toplevel)][-1]
picker.lift(); picker.attributes("-topmost", True)
capture(picker, os.path.join(OUT, "03_colour_wheel.png"))

for w in app.winfo_children():
    if isinstance(w, tk.Toplevel):
        w.destroy()

app.edit_bar_style()
app.update()
bar = [w for w in app.winfo_children() if isinstance(w, tk.Toplevel)][-1]
bar.lift(); bar.attributes("-topmost", True)
capture(bar, os.path.join(OUT, "04_bar_style.png"))
for w in app.winfo_children():
    if isinstance(w, tk.Toplevel):
        w.destroy()
tk.Toplevel.wait_window = _orig

# --- show a colour scheme + thicker, semi-transparent bars applied ---
app.scheme_name.set("Ocean")
app.plot.set_scheme("Ocean")
app.bar_width.set(5)
app.bar_trans.set(35)
app.plot.set_bar_style(5, 35)
app.lift(); app.attributes("-topmost", True)
capture(app, os.path.join(OUT, "05_colour_bar_applied.png"))
app.attributes("-topmost", False)

app.destroy()
print("done")
