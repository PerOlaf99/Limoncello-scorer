"""Manual genotyping for capillary-electrophoresis reads (built from scratch).

Click on a peak (or just next to it) and the best available peak-recognition
method locates the peak and shades its area.  The recognition algorithm is
selectable; a polymerase-A-addition (+A) and stutter peaks that trail/lead the
main peak by one repeat can be tagged automatically.

The four internal-standard peaks are captured through the CTC-CE duplex
pattern: peaks 1-2 are the two homoduplexes (they differ by the single SNP
base of the rs number) and peaks 3-4 are the two heteroduplexes made in the
PCR by Watson/Crick re-annealing (one mismatch base).  All four are the SAME
fragment (same bp), so a shared ``Fragment length (bp)`` applies to all four.

Peak position, channel/base, height, area and kind are collected in a table
that can be saved as CSV, Excel (.xlsx) or JSON — readable in Excel and usable
as ML training input.  For positions showing two peaks, the mutant (variant)
fraction is computed as small / (small + large).

The recognition logic lives in ``PeakPicker`` (headless, reusable from the main
window for batch picking); ``GenotypingEditor`` is the Tk widget that edits one
trace and ``GenotypingDialog`` the standalone Toplevel wrapper kept for
scripts/tests.
"""
import csv
import json
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from analyzer_core import load_trace

CHANNEL_ORDER = "ACGT"
CLICK_RADIUS = 30          # scans searched around a click
HET_WINDOW = 8             # scans in which two mains count as heterozygote
STUTTER_FRAC = 0.05        # min height of a stutter/+A peak vs the main peak
MIN_PEAK_FRAC = 0.05       # peak candidates must stand off the segment floor

PEAK_FINDERS = [
    ("best", "Best (prominence + area)"),
    ("max", "Simple local maxima"),
    ("gauss", "Gaussian fit"),
]


class _Record(dict):
    """One picked peak. Writable dict so rows feed straight into export."""


class PeakPicker:
    """Headless click-to-pick peak recognition for one trace.

    Holds the picked records, the CTC-CE internal-standard duplexes, the
    optional shared fragment length and all recognition/detection logic.  The
    UI reads the same attributes (records, std, length_bp, col_color) and
    calls the same methods, so one engine can drive a dialog, an embedded
    widget or the main window directly.
    """

    def __init__(self, doc, path, colors=None, base_order=None,
                 theme_mode="base", show=None, include_sh=True,
                 show_d2=True, finder="best"):
        self.doc = doc
        self.path = Path(path)
        self.show = show or (lambda col: True)
        self.colors = dict(colors or {"A": "#00AA00", "C": "#0000DD",
                                      "G": "#111111", "T": "#DD0000"})
        self.base_order = (base_order or "ACTG").upper()
        self.theme_mode = theme_mode
        self.include_sh = bool(include_sh)
        self.show_d2 = bool(show_d2)
        self.finder = finder or "best"

        self.records: list[_Record] = []
        self._gid = 0
        self._reject = None
        self.std = None
        self.length_bp = None
        self.col_color = {}
        self._d2_cache: dict[int, tuple] = {}
        for ci, base in enumerate(self.base_order[:4]):
            if base in CHANNEL_ORDER:
                if self.theme_mode == "channel":
                    self.col_color[CHANNEL_ORDER.index(base)] = \
                        self.colors.get(CHANNEL_ORDER[ci], "#444444")
                else:
                    self.col_color[CHANNEL_ORDER.index(base)] = \
                        self.colors.get(base, "#444444")

    # ------------------------------------------------------------- detection
    def pick(self, scan, vol=None):
        """Pick the peak you clicked on: nearest scan wins, then nearest
        voltage (tells channels apart when two share a scan), then tallest.
        A peak whose area is already picked on the same channel is never
        picked again — undo it first if you need to (neighbouring peaks, e.g.
        the two alleles of a heterozygote, stay pickable).  Returns the new
        record (plus any stutter/+A records appended) or None; on a refused
        re-pick, ``_reject`` is set to ``"area"`` for the UI message."""
        self._reject = None
        radius = max(CLICK_RADIUS, int(self._spacing() * 2.0))
        cands = []
        for col, color in self.col_color.items():
            if not self.show(col):
                continue
            pk = self._detect(col, scan, radius)
            if pk is None:
                continue
            pk["col"] = col
            pk["color"] = color
            cands.append(pk)
        if not cands:
            self._reject = "none"
            return None
        if vol is None:
            cands.sort(key=lambda p: (abs(p["apex"] - scan), -p["height"]))
        else:
            cands.sort(key=lambda p: (abs(p["apex"] - scan),
                                      abs(p["height"] - float(vol)),
                                      -p["height"]))
        best = cands[0]
        tol = max(1, int(round(self._spacing() * 0.25)))
        for rec in self.records:
            if rec["col"] == best["col"] and abs(rec["scan"] - best["apex"]) <= tol:
                self._reject = "area"
                return None
        self._gid += 1
        rec = _Record(
            file=str(self.path),
            well=self.doc.well,
            scan=best["apex"],
            channel=best["col"] + 1,
            base=CHANNEL_ORDER[best["col"]],
            kind="main",
            height=best["height"],
            area=best["area"],
            left=best["left"],
            right=best["right"],
            onset=best.get("onset", best["left"]),
            end=best.get("end", best["right"]),
            color=best["color"],
            col=best["col"],
            gid=self._gid,
        )
        if self.include_sh:
            for sib in self._shoulders(best, radius):
                rec2 = _Record(
                    file=str(self.path), well=self.doc.well,
                    scan=sib["apex"], channel=sib["col"] + 1,
                    base=CHANNEL_ORDER[sib["col"]], kind=sib["kind"],
                    height=sib["height"], area=sib["area"],
                    left=sib["left"], right=sib["right"], color=sib["color"],
                    col=sib["col"], onset=sib.get("onset", sib["left"]),
                    end=sib.get("end", sib["right"]),
                    gid=self._gid,
                )
                self.records.append(rec2)
        self.records.append(rec)
        return rec

    def _detect(self, col, scan, radius):
        """Run the chosen peak-recognition algorithm near scan on one column."""
        y = np.asarray(self.doc.acgt[:, col], dtype=float)
        n = y.size
        if n < 3:
            return None
        c = int(np.clip(scan, 0, n - 1))
        w = np.arange(max(0, c - radius), min(n, c + radius + 1))
        if w.size < 3:
            return None
        alg = dict(PEAK_FINDERS).get(self.finder, "best")
        if alg == "gauss":
            v = self._numeric_peak(y, w, col)
            if v is None:
                return v
            g = self._gauss_fit(y, v["left"], v["right"], v["apex"], v["height"])
            if g is not None:
                v.update(g)
            return v
        return self._numeric_peak(y, w, col)

    def _d2(self, col):
        """Smoothed second derivative of one channel, plus a robust noise
        estimate (median-absolute-deviation * 1.4826).  Cached per channel."""
        if col in self._d2_cache:
            return self._d2_cache[col]
        y = np.asarray(self.doc.acgt[:, col], dtype=float)
        try:
            from scipy.signal import savgol_filter
            yy = savgol_filter(y, window_length=7, polyorder=2, mode="interp")
        except Exception:
            k = np.ones(5) / 5.0
            yy = np.convolve(y, k, mode="same")
            yy = np.convolve(yy, k, mode="same")
        d2 = np.gradient(np.gradient(np.asarray(yy, dtype=float)))
        s2 = 1.4826 * np.median(np.abs(d2 - np.median(d2)))
        self._d2_cache[col] = (d2, float(s2))
        return self._d2_cache[col]

    def _onset_end(self, col, apex, left, right):
        """Start/stop of the peak from the second derivative: the concave-up
        lift-off on the rising flank (onset) and the concave-up return on the
        falling flank (end).  Falls back to the valley boundaries (left/right)
        when no curvature threshold stands out."""
        d2, s2 = self._d2(col)
        th = 3.0 * s2
        apex = int(apex)
        left = max(0, int(left))
        right = min(self.doc.acgt.shape[0] - 1, int(right))
        onset = None
        for i in range(left, max(left, apex - 1)):
            if d2[i] > th and d2[i + 1] > th:
                onset = i
                break
        end = None
        for i in range(right, max(apex + 1, right - 1), -1):
            if d2[i] > th and d2[i - 1] > th:
                end = i
                break
        return (left if onset is None else onset,
                right if end is None else end)

    def _numeric_peak(self, y, w, col):
        seg = y[w]
        base = float(np.nanmin(seg))
        top = float(np.nanmax(seg))
        rng = top - base
        if not np.isfinite(rng) or rng <= 1e-12:
            return None
        m = (seg[1:-1] >= seg[:-2]) & (seg[1:-1] > seg[2:])
        idx = np.flatnonzero(m) + 1
        if idx.size == 0:
            return None
        on = idx[seg[idx] - base >= MIN_PEAK_FRAC * rng]
        if on.size == 0:
            on = idx
        # strongest candidate nearest to the window centre
        c = w.size // 2
        dist = np.abs(w[on] - w[c])
        order = np.lexsort((seg[on], dist))
        best = int(on[order[0]])
        apex = int(w[best])
        left = int(w[:best + 1][np.argmin(seg[:best + 1])])
        tail = seg[best:]
        right = int(w[best + np.argmin(tail)])
        if right <= left:
            right = min(y.size - 1, left + 2)
        xs = np.arange(left, right + 1)
        bl = np.linspace(float(y[left]), float(y[right]), right - left + 1)
        area = float(np.sum(np.clip(y[xs] - bl, 0.0, None)))
        height = float(y[apex] - max(y[left], y[right]))
        if height <= 0 or area <= 0:
            return None
        onset, end = self._onset_end(col, apex, left, right)
        return {"apex": apex, "left": left, "right": right,
                "onset": onset, "end": end,
                "area": area, "height": height}

    def _gauss_fit(self, y, left, right, apex, height):
        from scipy.optimize import curve_fit
        xs = np.arange(left, right + 1)
        if xs.size < 5:
            return None

        def gauss(x, a, mu, s, b):
            return a * np.exp(-(x - mu) ** 2 / (2.0 * s * s)) + b

        p0 = (float(height), float(apex), max(2.0, (right - left) / 3.0),
              float(y[apex] - height))
        try:
            popt, _ = curve_fit(gauss, xs, y[xs], p0=p0, maxfev=4000)
        except Exception:
            return None
        a, mu, s, _b = popt
        if not (np.isfinite(a) and np.isfinite(mu) and np.isfinite(s)) or s <= 0:
            return None
        return {"apex": int(round(mu)),
                "area": float(a * abs(s) * np.sqrt(2.0 * np.pi)),
                "height": float(a)}

    def _shoulders(self, main, radius):
        """Tag the strongest satellite on each side of the main peak: the
        upstream one (stutter) and the downstream one (+A A-addition), found
        within ~±(0.4–1.6) × one repeat (≈ one base) of the main apex."""
        sp = self._spacing()
        y = np.asarray(self.doc.acgt[:, main["col"]], dtype=float)
        n = y.size
        out = []
        for k, kind in ((+1, "+A"), (-1, "stutter")):
            a, b = (main["apex"] + k * sp * 0.4,
                    main["apex"] + k * sp * 1.6)
            lo, hi = (int(min(a, b)), int(max(a, b)))
            cands = self._window_peaks(main["col"], lo, hi)
            cands = [(x, h) for x, h in cands
                     if STUTTER_FRAC * main["height"] <= h
                     <= 0.9 * main["height"]]
            if not cands:
                continue
            x, h = max(cands, key=lambda c: c[1])
            if k == +1 and x <= main["apex"] + 2:      # must trail the main
                continue
            if k == -1 and x >= main["apex"] - 2:      # must lead the main
                continue
            r = max(2, int(round(sp * 0.20)))
            left = int(np.argmin(y[max(0, x - r): x + 1])) + max(0, x - r)
            right = int(np.argmin(y[x: min(n, x + r + 1)])) + x
            if right <= left:
                right = min(n - 1, left + 2)
            xs = np.arange(left, right + 1)
            bl = np.linspace(y[left], y[right], right - left + 1)
            area = float(np.sum(np.clip(y[xs] - bl, 0.0, None)))
            hgt = float(y[x] - max(y[left], y[right]))
            if hgt <= 0 or area <= 0:
                continue
            onset, end = self._onset_end(main["col"], x, left, right)
            out.append({"apex": x, "left": left, "right": right,
                        "onset": onset, "end": end,
                        "area": area, "height": hgt,
                        "col": main["col"], "color": main["color"],
                        "kind": kind})
        return out

    def _window_peaks(self, col, lo, hi):
        """All local maxima in a scan window of one channel, as (scan, value)."""
        y = np.asarray(self.doc.acgt[:, col], dtype=float)
        lo = max(0, int(lo))
        hi = min(y.size - 1, int(hi))
        if hi - lo < 3:
            return []
        w = np.arange(lo, hi + 1)
        seg = y[w]
        rng = float(seg.max()) - float(seg.min())
        if not np.isfinite(rng) or rng <= 1e-12:
            return []
        m = (seg[1:-1] >= seg[:-2]) & (seg[1:-1] > seg[2:])
        idx = [int(i) + 1 for i in np.flatnonzero(m)
               if seg[int(i) + 1] - seg.min() >= MIN_PEAK_FRAC * rng]
        return [(int(w[i]), float(seg[i])) for i in idx]

    def _spacing(self):
        p = np.asarray(getattr(self.doc, "peak_positions", []) or [], dtype=float)
        if p.size > 1:
            d = float(np.median(np.diff(p)))
            if np.isfinite(d) and 2.0 < d < 80.0:
                return d
        return float(max(6.0, self.doc.n_scans * 0.004))

    # ------------------------------------- internal-standard CTC-CE duplexes
    def mark_std(self, length_bp=None):
        """Tag the standard mains (scan order) as the four CTC-CE duplexes:
        HOM1, HOM2 (homoduplexes, one rs SNP base apart), HET1, HET2
        (heteroduplexes, one mismatch base from Watson/Crick re-annealing).
        All four are the same fragment, so a shared length (bp) is optional.
        The standard channel is the one of the last main you picked.

        Raises ValueError with a plain-language reason, or returns a message."""
        mains = [r for r in self.records if r["kind"] == "main"]
        if not mains:
            raise ValueError("Pick the standard main peaks first.")
        last_col = mains[-1]["col"]
        col_mains = sorted((m for m in mains if m["col"] == last_col),
                           key=lambda m: m["scan"])
        if len(col_mains) < 2:
            raise ValueError("Need ≥ 2 standard main peaks on the same channel.")
        names = ["HOM1", "HOM2", "HET1", "HET2"][:min(len(col_mains), 4)]
        self.std = [(m["scan"], n) for m, n in zip(col_mains[:4], names)]
        self.length_bp = length_bp
        msg = "Standard set: " + ", ".join(f"{n}@{x}" for x, n in self.std)
        if self.length_bp is not None:
            msg += f"  (len {self.length_bp:g} bp)"
        return msg

    def clear_std(self):
        self.std = None
        self.length_bp = None

    def duplex_of(self, rec):
        """Duplex label (HOM1/HOM2/HET1/HET2) for a main peak that is one of
        the marked standard peaks; '' otherwise."""
        if not self.std or rec["kind"] != "main":
            return ""
        for x, n in self.std:
            if abs(rec["scan"] - x) <= 2:
                return n
        return ""

    def _duplex_of(self, rec):
        """Back-compat alias for duplex_of (used by the editor table/tests)."""
        return self.duplex_of(rec)

    def _clust_frac(self, rec):
        """Back-compat alias for clust_frac (used by the editor table/tests)."""
        return self.clust_frac(rec)

    def _het_window(self):
        """How many scans count as one allelic position: a repeat (~ one base,
        from the channel spacing) with a little slack, never below HET_WINDOW.
        On a CTC-CE run the two alleles of a heterozygote sit one base apart,
        which is typically far more than the old fixed 8-scan window."""
        return max(HET_WINDOW, int(round(self._spacing() * 1.6)))

    def clusters(self):
        """Main peaks grouped into allelic positions: same channel and within
        one repeat (≈ one base) of another main in the group."""
        mains = [r for r in self.records if r["kind"] == "main"]
        win = self._het_window()
        out = []
        for r in sorted(mains, key=lambda m: m["scan"]):
            placed = False
            for cl in out:
                if any(r["col"] == m["col"]
                       and abs(r["scan"] - m["scan"]) <= win for m in cl):
                    cl.append(r)
                    placed = True
                    break
            if not placed:
                out.append([r])
        return out

    def clust_frac(self, rec):
        """Variant (mutant) fraction small/(small+large) for a main peak in a
        position with two peaks; 0.0 otherwise."""
        for cl in self.clusters():
            if rec in cl and len(cl) >= 2:
                areas = sorted(m["area"] for m in cl)
                if sum(areas) > 0:
                    return areas[0] / sum(areas)
        return 0.0

    # ---------------------------------------------------------------- mutate
    def undo_last(self):
        """Remove the most recently picked peak (plus its stutter/+A tags)."""
        if not self.records:
            return False
        gid = max(r["gid"] for r in self.records)
        self.records = [r for r in self.records if r["gid"] != gid]
        return True

    def clear_all(self):
        self.records = []
        self._d2_cache = {}

    # ---------------------------------------------------------------- export
    def export_rows(self):
        """One writable dict per picked peak, ready for CSV/JSON/XLSX."""
        rows = []
        for r in self.records:
            rows.append({
                "file": r["file"], "well": r["well"], "scan": r["scan"],
                "channel": r["channel"], "base": r["base"], "kind": r["kind"],
                "start_scan": r.get("onset"), "end_scan": r.get("end"),
                "height_V": round(r["height"], 4),
                "area_Vscan": round(r["area"], 3),
                "duplex": self.duplex_of(r),
                "length_bp": self.length_bp,
                "fraction": round(self.clust_frac(r), 4)
                if r["kind"] == "main" else "",
            })
        return rows

    # ---------------------------------------------------------------- overlay
    def plot_overlay(self, ax):
        """Draw the picked records onto a matplotlib axes: shaded area, apex
        ring, 2nd-derivative start/end marks and a small label."""
        n = self.doc.acgt.shape[0]
        tick = float(np.nanmax(self.doc.acgt)) * 0.03
        for r in self.records:
            if not (0 <= r["left"] <= r["right"] < n):
                continue
            xs = np.arange(r["left"], r["right"] + 1)
            y = self.doc.acgt[xs, r["col"]]
            ax.fill_between(xs, 0, y, color=r["color"], alpha=0.25, zorder=1)
            ax.plot([r["scan"]], [self.doc.acgt[r["scan"], r["col"]]],
                    marker="o", ms=5, mfc="none", mec=r["color"], zorder=4)
            if self.show_d2:
                for xx in (r["onset"], r["end"]):
                    if xx is None or not (0 <= int(xx) < n):
                        continue
                    xx = int(xx)
                    yy = self.doc.acgt[xx, r["col"]]
                    ax.vlines(xx, max(0, yy - tick), yy + tick,
                              color=r["color"], lw=1.0, zorder=4)
                    ax.plot([xx], [yy], marker="s", ms=2.5, mfc=r["color"],
                            mec=r["color"], zorder=4)
            txt = f"{r['scan']}·{r['base']}"
            if r["kind"] != "main":
                txt = f"{r['kind']} " + txt
            ax.annotate(txt, (r["scan"], y.max()), textcoords="offset points",
                        xytext=(0, 6), fontsize=6, color=r["color"],
                        ha="center", zorder=5)


def _write_csv(path, rows):
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(list(rows[0].keys()))
        for r in rows:
            w.writerow(list(r.values()))


def _write_xlsx(path, rows):
    from openpyxl import Workbook
    from openpyxl.utils import get_column_letter
    wb = Workbook()
    ws = wb.active
    ws.title = "peaks"
    headers = list(rows[0].keys())
    ws.append(headers)
    for r in rows:
        ws.append([r[h] for h in headers])
    for i, h in enumerate(headers, 1):
        ws.column_dimensions[get_column_letter(i)].width = \
            max(8, min(28, 6 + len(h)))
    ws.freeze_panes = "A2"
    wb.save(path)


def save_table(path, rows):
    """Write peak rows to CSV, JSON or XLSX depending on the file suffix."""
    ext = Path(path).suffix.lower()
    if ext == ".json":
        Path(path).write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    elif ext == ".xlsx":
        _write_xlsx(path, rows)
    else:
        _write_csv(path, rows)


class GenotypingEditor(ttk.Frame):
    """Click-to-pick peak editor. Usable as a standalone Toplevel widget pack
    (GenotypingDialog wraps it); channels can be switched off one at a time.
    The recognition engine is a PeakPicker, so the editor and the main window
    share exactly the same picking behaviour."""

    def __init__(self, master, path, colors=None, base_order=None,
                 theme_mode="base", show=None, on_exit=None):
        super().__init__(master)
        self.show = show or (lambda col: True)
        self.path = Path(path)
        self.colors = dict(colors or {"A": "#00AA00", "C": "#0000DD",
                                      "G": "#111111", "T": "#DD0000"})
        self.base_order = (base_order or "ACTG").upper()
        self.theme_mode = theme_mode
        self.on_exit = on_exit

        try:
            self.doc = load_trace(self.path, base_order=self.base_order)
        except Exception as e:
            raise ValueError(f"Could not load trace: {e}")

        self.include_sh = tk.BooleanVar(value=True)
        self.show_d2 = tk.BooleanVar(value=True)
        self.finder = tk.StringVar(value="best")
        self.pk = PeakPicker(self.doc, self.path, colors=self.colors,
                             base_order=self.base_order,
                             theme_mode=self.theme_mode, show=self.show,
                             include_sh=self.include_sh.get(),
                             show_d2=self.show_d2.get())

        self._build()
        self.redraw()

    def __getattr__(self, name):
        pk = self.__dict__.get("pk")
        if pk is not None:
            return getattr(pk, name)
        raise AttributeError(name)

    @property
    def records(self):
        return self.pk.records

    # ------------------------------------------------------------------ UI
    def _build(self):
        pane = ttk.Panedwindow(self, orient=tk.HORIZONTAL)
        pane.pack(fill=tk.BOTH, expand=True)

        left = ttk.Frame(pane)
        pane.add(left, weight=3)
        self.fig = Figure(figsize=(9, 5.5), dpi=100, facecolor="#FFFFFF")
        self.fig.patch.set_facecolor("#FFFFFF")
        self.canvas = FigureCanvasTkAgg(self.fig, master=left)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, window=left)
        self.toolbar.update()
        self.canvas.mpl_connect("button_press_event", self._on_click)

        right = ttk.Frame(pane, width=360)
        pane.add(right, weight=1)
        info = ttk.LabelFrame(right, text="Run", padding=6)
        info.pack(fill=tk.X, padx=4, pady=4)
        ttk.Label(info, text=str(self.path), wraplength=330, justify=tk.LEFT).pack(anchor=tk.W)
        sub = (f"{self.doc.n_scans} scans, {self.doc.source.lower()} "
               f"· base order {self.doc.base_order}"
               + (f" · called {len(self.doc.sequence)} bp"
                  if self.doc.sequence else " · not called"))
        ttk.Label(info, text=sub, foreground="#555").pack(anchor=tk.W)

        det = ttk.LabelFrame(right, text="Peak recognition", padding=6)
        det.pack(fill=tk.X, padx=4, pady=4)
        self.finder.trace_add("write", lambda *a: self._sync_opts(redraw=False))
        ttk.Label(det, text="Algorithm").pack(anchor=tk.W)
        ttk.Combobox(det, textvariable=self.finder, state="readonly",
                     values=[v for _, v in PEAK_FINDERS], width=30).pack(anchor=tk.W)
        self.include_sh = tk.BooleanVar(value=True)
        ttk.Checkbutton(det, text="Add stutter & +A (A-addition) peaks",
                        variable=self.include_sh,
                        command=self._sync_opts).pack(anchor=tk.W)
        self.show_d2 = tk.BooleanVar(value=True)
        ttk.Checkbutton(det, text="Mark start/end from the 2nd derivative",
                        variable=self.show_d2,
                        command=self._sync_opts).pack(anchor=tk.W)
        self.het_label = ttk.Label(det, text="Click a peak (or near it) to pick it.",
                                   foreground="#0F3A6E", wraplength=320, justify=tk.LEFT)
        self.het_label.pack(anchor=tk.W, pady=(4, 0))

        std = ttk.LabelFrame(right, text="Internal standard (CTC-CE duplex pattern)",
                             padding=6)
        std.pack(fill=tk.X, padx=4, pady=4)
        fw = 330
        ttk.Label(std, text="The four standard peaks are the SAME fragment (same bp),\n"
                            "separated by cycling-temperature CE according to sequence:\n"
                            "peaks 1-2 = the two homoduplexes (they differ by the single\n"
                            "SNP base of the rs number); peaks 3-4 = the two heteroduplexes\n"
                            "made in the PCR by Watson/Crick re-annealing (one mismatch base).",
                  justify=tk.LEFT, foreground="#555", wraplength=fw).pack(anchor=tk.W, pady=(0, 4))
        sz = ttk.Frame(std)
        sz.pack(fill=tk.X)
        ttk.Label(sz, text="Fragment length (bp)").pack(side=tk.LEFT)
        self.len_entry = ttk.Entry(sz, width=8, justify=tk.RIGHT)
        self.len_entry.pack(side=tk.LEFT, padx=4)
        ttk.Label(sz, text="optional — same for all four",
                  foreground="#888").pack(side=tk.LEFT)
        sbtn = ttk.Frame(std)
        sbtn.pack(fill=tk.X, pady=(4, 0))
        ttk.Button(sbtn, text="Mark picked peaks as standard",
                   command=self._mark_std).pack(side=tk.LEFT)
        ttk.Button(sbtn, text="Clear", command=self._clear_std).pack(side=tk.LEFT, padx=4)
        ttk.Label(std, text="The four first main peaks (scan order) are tagged\n"
                            "HOM1, HOM2 (homoduplexes) then HET1, HET2 (heteroduplexes).",
                  foreground="#777", wraplength=fw, justify=tk.LEFT).pack(anchor=tk.W, pady=(4, 0))
        self.std_lbl = ttk.Label(std, text="No standard set.",
                                 foreground="#555", wraplength=fw, justify=tk.LEFT)
        self.std_lbl.pack(anchor=tk.W, pady=(4, 0))

        bars = ttk.Frame(right)
        bars.pack(fill=tk.X, padx=4, pady=2)
        ttk.Button(bars, text="Undo last",
                   command=self._undo_last).pack(side=tk.LEFT)
        ttk.Button(bars, text="Clear all",
                   command=self._clear_all).pack(side=tk.LEFT, padx=4)

        tblf = ttk.LabelFrame(right, text="Picked peaks",
                              padding=4)
        tblf.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)
        cols = ("#", "scan", "duplex", "ch", "kind", "height V", "area V·sc", "frac")
        self.tree = ttk.Treeview(tblf, columns=cols, show="headings", height=12)
        widths = {"#": 34, "scan": 54, "duplex": 58, "ch": 40, "kind": 66,
                  "height V": 70, "area V·sc": 78, "frac": 50}
        for c in cols:
            self.tree.heading(c, text=c)
            self.tree.column(c, width=widths[c], anchor="e" if c not in ("#", "kind") else "w",
                             stretch=(c in ("scan", "kind")))
        vs = ttk.Scrollbar(tblf, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.configure(yscrollcommand=vs.set)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vs.pack(side=tk.RIGHT, fill=tk.Y)

        http_bar = ttk.Frame(right)
        http_bar.pack(fill=tk.X, padx=4, pady=4)
        ttk.Button(http_bar, text="Save table…", command=self._save).pack(side=tk.LEFT)
        ttk.Button(http_bar, text="Close / back to viewer",
                   command=self._request_close).pack(side=tk.RIGHT)

        self.frac_var = tk.StringVar(value="No positions with two peaks yet.")
        fout = ttk.Label(right, textvariable=self.frac_var,
                         foreground="#555", wraplength=340, justify=tk.LEFT)
        fout.pack(fill=tk.X, padx=4, pady=(0, 6))

    # ------------------------------------------------------------- picking
    def _on_click(self, event):
        if event.xdata is None or event.inaxes is None:
            return
        if self.toolbar.mode != "":
            return
        self._pick_peak(int(round(event.xdata)), vol=event.ydata)

    def _pick_peak(self, scan, vol=None):
        rec = self.pk.pick(scan, vol=vol)
        if rec is None:
            if self.pk._reject == "area":
                self._status("That area is already picked — undo it first to "
                             "pick it again.")
            else:
                self._status("No peak found near that scan — try again closer "
                             "to a hump.")
            return
        self._status(f"Peak at scan {rec['scan']} · {rec['base']}"
                     f" height {rec['height']:.3f} V, area {rec['area']:.2f} V·scan")
        self.redraw()
        self._sync_table()

    def _sync_opts(self, redraw=True):
        self.pk.include_sh = bool(self.include_sh.get())
        self.pk.show_d2 = bool(self.show_d2.get())
        self.pk.finder = self.finder.get()
        if redraw:
            self.redraw()
            self._sync_table()

    def _request_close(self):
        """Leave peak picking: back to the main viewer when embedded."""
        if self.on_exit is not None:
            self.on_exit()
        else:
            self.destroy()

    def _close(self):
        self._request_close()

    def _status(self, msg):
        self.het_label.config(text=msg)

    # ------------------------------------- internal-standard CTC-CE duplexes
    def _mark_std(self):
        try:
            val = float(self.len_entry.get())
        except ValueError:
            val = None
        try:
            msg = self.pk.mark_std(val)
        except ValueError as e:
            self.std_lbl.config(text=str(e), foreground="#A33")
            return
        self._sync_table()

    def _clear_std(self):
        self.pk.clear_std()
        self.std_lbl.config(text="No standard set.", foreground="#555")
        self._sync_table()

    # ---------------------------------------------------------------- view
    def redraw(self):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        n = self.doc.acgt.shape[0]
        x = np.arange(n)
        for col, color in self.pk.col_color.items():
            if not self.show(col):
                continue
            ax.plot(x, self.doc.acgt[:, col], color=color, lw=0.7,
                    label=f"Ch{col + 1} {CHANNEL_ORDER[col]}")
        # base letters of an existing call, faint, for orientation
        if self.doc.sequence and self.doc.peak_positions:
            seq = self.doc.sequence
            pos = np.asarray(self.doc.peak_positions, dtype=int)
            top = float(np.nanmax(self.doc.acgt))
            for pi, b in enumerate(seq):
                if pi >= len(pos) or b not in CHANNEL_ORDER:
                    continue
                if not self.show(CHANNEL_ORDER.index(b)):
                    continue
                ax.text(pos[pi], top * 1.01, b, ha="center", va="bottom",
                        fontsize=5, color=self.colors.get(b, "#444"), alpha=0.85, zorder=2)
        self.pk.plot_overlay(ax)
        ax.set_xlim(0, n)
        ax.set_ylim(0, float(np.nanmax(self.doc.acgt)) * 1.08 or 1.0)
        ax.set_xlabel("scan")
        ax.set_ylabel("V")
        ax.set_title(f"{self.path.name} — click to pick a peak", fontsize=9)
        ax.grid(True, alpha=0.15)
        ax.legend(loc="upper right", fontsize=7, ncol=2, framealpha=0.6)
        self.fig.tight_layout()
        self.canvas.draw_idle()
        self._sync_fractions()

    def _sync_fractions(self):
        lines = []
        for cl in self.pk.clusters():
            if len(cl) >= 2 and sum(m["area"] for m in cl) > 0:
                areas = sorted(m["area"] for m in cl)
                frac = areas[0] / sum(areas)
                chans = "+".join(sorted(m["base"] for m in cl))
                lines.append(f"scan {cl[0]['scan']}  {chans}:  "
                             f"variant {frac:.3f}")
        self.frac_var.set("\n".join(lines) if lines
                          else "No position with two peaks yet (mutant "
                               "fraction appears here).")

    # ---------------------------------------------------------------- table
    def _sync_table(self):
        self.tree.delete(*self.tree.get_children())
        for i, r in enumerate(self.records, 1):
            fr = self._clust_frac(r) if r["kind"] == "main" else ""
            self.tree.insert("", tk.END, values=(
                i, r["scan"], self._duplex_of(r), r["base"], r["kind"],
                f"{r['height']:.3f}", f"{r['area']:.1f}",
                f"{fr:.3f}" if fr else ""))

    def _undo_last(self):
        if not self.pk.undo_last():
            return
        self.redraw()
        self._sync_table()
        self._status("Removed last picked peak.")

    def _clear_all(self):
        self.pk.clear_all()
        self.redraw()
        self._sync_table()
        self._status("Table cleared.")

    # -------------------------------------------------------------- export
    def _save(self):
        if not self.records:
            messagebox.showinfo("Save table", "Pick some peaks first.", parent=self)
            return
        types = [("CSV (Excel-compatible)", "*.csv"), ("JSON (ML)", "*.json")]
        try:
            import openpyxl  # noqa: F401
            types.insert(1, ("Excel workbook (.xlsx)", "*.xlsx"))
        except ImportError:
            pass
        path = filedialog.asksaveasfilename(parent=self, defaultextension=".csv",
                                            filetypes=types)
        if not path:
            return
        rows = self.pk.export_rows()
        try:
            save_table(path, rows)
        except Exception as e:
            messagebox.showerror("Save table", f"Could not write file:\n{e}", parent=self)
            return
        self._status(f"Saved {len(rows)} peak rows to {path}")


class GenotypingDialog(tk.Toplevel):
    """Standalone window wrapping a GenotypingEditor (kept for scripts/tests;
    the app picks peaks directly in the main window instead, batching several
    samples at a time)."""

    def __init__(self, parent, path, colors=None, base_order=None,
                 theme_mode="base"):
        super().__init__(parent)
        self.title("Manual genotyping — " + Path(path).name)
        self.geometry("1080x680")
        try:
            ed = GenotypingEditor(self, path, colors=colors, base_order=base_order,
                                  theme_mode=theme_mode)
        except Exception as e:
            messagebox.showerror("Genotyping", f"Could not load trace:\n{e}",
                                 parent=self)
            self.destroy()
            return
        self._editor = ed
        ed.pack(fill=tk.BOTH, expand=True)

    def __getattr__(self, name):
        ed = self.__dict__.get("_editor")
        if ed is None:
            raise AttributeError(name)
        return getattr(ed, name)