"""Manual genotyping for capillary-electrophoresis reads (built from scratch).

Click on a peak (or just next to it) and the best available peak-recognition
method locates the peak and shades its area.  The recognition algorithm is
selectable; a polymerase-A-addition (+A) and stutter peaks that trail/lead the
main peak by one repeat can be tagged automatically.

Peak position, channel/base, height, area and kind are collected in a table
that can be saved as CSV, Excel (.xlsx) or JSON — readable in Excel and usable
as ML training input.  For positions showing two peaks, the mutant (variant)
fraction is computed as  small / (small + large).
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


class GenotypingDialog(tk.Toplevel):
    """A standalone window: the selected well's traces, click-to-pick peaks,
    area shading, stutter/+A tagging, a results table and Excel/CSV/JSON out."""

    def __init__(self, parent, path, colors=None, base_order=None,
                 theme_mode="base"):
        super().__init__(parent)
        self.parent = parent
        self.path = Path(path)
        self.colors = dict(colors or {"A": "#00AA00", "C": "#0000DD",
                                      "G": "#111111", "T": "#DD0000"})
        self.base_order = (base_order or "ACTG").upper()
        self.theme_mode = theme_mode

        try:
            self.doc = load_trace(self.path, base_order=self.base_order)
        except Exception as e:
            messagebox.showerror("Genotyping", f"Could not load trace:\n{e}", parent=self)
            self.destroy()
            return

        self.records: list[_Record] = []
        self._gid = 0
        self.std = None
        self.length_bp = None
        self.col_color = {}
        for ci, base in enumerate(self.base_order[:4]):
            if base in CHANNEL_ORDER:
                if self.theme_mode == "channel":
                    self.col_color[CHANNEL_ORDER.index(base)] = \
                        self.colors.get(CHANNEL_ORDER[ci], "#444444")
                else:
                    self.col_color[CHANNEL_ORDER.index(base)] = \
                        self.colors.get(base, "#444444")

        self.title("Manual genotyping — " + self.path.name)
        self.geometry("1080x680")
        self._build()
        self.redraw()

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
        self.finder = tk.StringVar(value="best")
        ttk.Label(det, text="Algorithm").pack(anchor=tk.W)
        ttk.Combobox(det, textvariable=self.finder, state="readonly",
                     values=[v for _, v in PEAK_FINDERS], width=30).pack(anchor=tk.W)
        self.include_sh = tk.BooleanVar(value=True)
        ttk.Checkbutton(det, text="Add stutter & +A (A-addition) peaks",
                        variable=self.include_sh,
                        command=lambda: (self.redraw(), self._sync_table())).pack(anchor=tk.W)
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
        ttk.Button(http_bar, text="Close", command=self.destroy).pack(side=tk.RIGHT)

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
        """Pick the peak you clicked on: nearest scan wins, then nearest
        voltage (tells channels apart when two share a scan), then tallest."""
        radius = max(CLICK_RADIUS, int(self._spacing() * 2.0))
        cands = []
        for col, color in self.col_color.items():
            pk = self._detect(col, scan, radius)
            if pk is None:
                continue
            pk["col"] = col
            pk["color"] = color
            cands.append(pk)
        if not cands:
            self._status("No peak found near that scan — try again closer to a hump.")
            return
        if vol is None:
            cands.sort(key=lambda p: (abs(p["apex"] - scan), -p["height"]))
        else:
            cands.sort(key=lambda p: (abs(p["apex"] - scan),
                                      abs(p["height"] - float(vol)),
                                      -p["height"]))
        best = cands[0]
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
            color=best["color"],
            col=best["col"],
            gid=self._gid,
        )
        if self.include_sh.get():
            for sib in self._shoulders(best, radius):
                rec2 = _Record(
                    file=str(self.path), well=self.doc.well,
                    scan=sib["apex"], channel=sib["col"] + 1,
                    base=CHANNEL_ORDER[sib["col"]], kind=sib["kind"],
                    height=sib["height"], area=sib["area"],
                    left=sib["left"], right=sib["right"], color=sib["color"],
                    col=sib["col"],
                    gid=self._gid,
                )
                self.records.append(rec2)
        self.records.append(rec)
        self._status(f"Peak at scan {rec['scan']} · {rec['base']}"
                     f" height {rec['height']:.3f} V, area {rec['area']:.2f} V·scan")
        self.redraw()
        self._sync_table()

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
        alg = dict(PEAK_FINDERS)[self.finder.get()]
        if alg == "gauss":
            v = self._numeric_peak(y, w)
            if v is None:
                return v
            g = self._gauss_fit(y, v["left"], v["right"], v["apex"], v["height"])
            if g is not None:
                v.update(g)
            return v
        return self._numeric_peak(y, w)

    def _numeric_peak(self, y, w):
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
        return {"apex": apex, "left": left, "right": right,
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
            out.append({"apex": x, "left": left, "right": right,
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

    def _status(self, msg):
        self.het_label.config(text=msg)

    # ------------------------------------- internal-standard CTC-CE duplexes
    def _mark_std(self):
        """Tag the standard mains (scan order) as the four CTC-CE duplexes:
        HOM1, HOM2 (homoduplexes, one rs SNP base apart), HET1, HET2
        (heteroduplexes, one mismatch base from Watson/Crick re-annealing).
        All four are the same fragment, so a shared length (bp) is optional.
        The standard channel is the one of the last main you picked."""
        mains = [r for r in self.records if r["kind"] == "main"]
        if not mains:
            self.std_lbl.config(text="Pick the standard main peaks first.",
                                foreground="#A33")
            return
        last_col = mains[-1]["col"]
        col_mains = sorted((m for m in mains if m["col"] == last_col),
                           key=lambda m: m["scan"])
        if len(col_mains) < 2:
            self.std_lbl.config(text="Need ≥ 2 standard main peaks on the same channel.",
                                foreground="#A33")
            return
        names = ["HOM1", "HOM2", "HET1", "HET2"][:min(len(col_mains), 4)]
        pairs = [(m["scan"], n) for m, n in zip(col_mains[:4], names)]
        self.std = pairs
        try:
            self.length_bp = float(self.len_entry.get())
        except ValueError:
            self.length_bp = None
        msg = "Standard set: " + ", ".join(f"{n}@{x}" for x, n in pairs)
        if self.length_bp is not None:
            msg += f"  (len {self.length_bp:g} bp)"
        self.std_lbl.config(text=msg, foreground="#0A5")
        self._sync_table()

    def _clear_std(self):
        self.std = None
        self.length_bp = None
        self.std_lbl.config(text="No standard set.", foreground="#555")
        self._sync_table()

    def _duplex_of(self, rec):
        """Duplex label (HOM1/HOM2/HET1/HET2) for a main peak that is one of
        the marked standard peaks; '' otherwise."""
        if not self.std or rec["kind"] != "main":
            return ""
        for x, n in self.std:
            if abs(rec["scan"] - x) <= 2:
                return n
        return ""

    # ---------------------------------------------------------------- view
    def redraw(self):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        n = self.doc.acgt.shape[0]
        x = np.arange(n)
        for col, color in self.col_color.items():
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
                ax.text(pos[pi], top * 1.01, b, ha="center", va="bottom",
                        fontsize=5, color=self.colors.get(b, "#444"), alpha=0.85, zorder=2)
        # shade each picked peak and label it
        for r in self.records:
            xs = np.arange(r["left"], r["right"] + 1)
            y = self.doc.acgt[xs, r["col"]]
            ax.fill_between(xs, 0, y, color=r["color"], alpha=0.25, zorder=1)
            ax.plot([r["scan"]], [self.doc.acgt[r["scan"], r["col"]]],
                    marker="o", ms=5, mfc="none", mec=r["color"], zorder=4)
            txt = f"{r['scan']}·{r['base']}"
            if r["kind"] != "main":
                txt = f"{r['kind']} " + txt
            ax.annotate(txt, (r["scan"], y.max()), textcoords="offset points",
                        xytext=(0, 6), fontsize=6, color=r["color"],
                        ha="center", zorder=5)
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
        mains = [r for r in self.records if r["kind"] == "main"]
        clusters = []
        for r in sorted(mains, key=lambda m: m["scan"]):
            placed = False
            for cl in clusters:
                if any(abs(r["scan"] - m["scan"]) <= HET_WINDOW for m in cl):
                    cl.append(r)
                    placed = True
                    break
            if not placed:
                clusters.append([r])
        lines = []
        for cl in clusters:
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
        # store per-scan fractions so rows show their cluster value
        frac_by_scan = {}
        for r in self.records:
            if r["kind"] == "main":
                frac_by_scan.setdefault(r["scan"], self._clust_frac(r))
        for i, r in enumerate(self.records, 1):
            fr = self._clust_frac(r) if r["kind"] == "main" else ""
            self.tree.insert("", tk.END, values=(
                i, r["scan"], self._duplex_of(r), r["base"], r["kind"],
                f"{r['height']:.3f}", f"{r['area']:.1f}",
                f"{fr:.3f}" if fr else ""))

    def _clust_frac(self, rec):
        for cl in self._clusters():
            if rec in cl and len(cl) >= 2:
                areas = sorted(m["area"] for m in cl)
                if sum(areas) > 0:
                    return areas[0] / sum(areas)
        return 0.0

    def _clusters(self):
        mains = [r for r in self.records if r["kind"] == "main"]
        out = []
        for r in sorted(mains, key=lambda m: m["scan"]):
            placed = False
            for cl in out:
                if any(abs(r["scan"] - m["scan"]) <= HET_WINDOW for m in cl):
                    cl.append(r)
                    placed = True
                    break
            if not placed:
                out.append([r])
        return out

    def _undo_last(self):
        """"Remove the most recently picked peak (plus its stutter/+A tags)."""
        if not self.records:
            return
        gid = max(r["gid"] for r in self.records)
        self.records = [r for r in self.records if r["gid"] != gid]
        self.redraw()
        self._sync_table()
        self._status("Removed last picked peak.")

    def _clear_all(self):
        self.records = []
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
        rows = []
        for r in self.records:
            rows.append({
                "file": r["file"], "well": r["well"], "scan": r["scan"],
                "channel": r["channel"], "base": r["base"], "kind": r["kind"],
                "height_V": round(r["height"], 4),
                "area_Vscan": round(r["area"], 3),
                "duplex": self._duplex_of(r),
                "length_bp": self.length_bp,
                "fraction": round(self._clust_frac(r), 4)
                if r["kind"] == "main" else "",
            })
        try:
            ext = Path(path).suffix.lower()
            if ext == ".json":
                Path(path).write_text(
                    json.dumps(rows, indent=2) + "\n", encoding="utf-8")
            elif ext == ".xlsx":
                self._write_xlsx(path, rows)
            else:
                self._write_csv(path, rows)
        except Exception as e:
            messagebox.showerror("Save table", f"Could not write file:\n{e}", parent=self)
            return
        self._status(f"Saved {len(rows)} peak rows to {path}")

    def _write_csv(self, path, rows):
        with open(path, "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.writer(fh)
            w.writerow(list(rows[0].keys()))
            for r in rows:
                w.writerow(list(r.values()))

    def _write_xlsx(self, path, rows):
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