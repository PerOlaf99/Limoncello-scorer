"""Plate overlay: every well of a run in one graph, aligned on its standard.

Each well's sample (or standard) channel is put on the plate's common time
axis through the internal-standard bands the SNP engine places in it
(:func:`snp_genotyping.align_plate`), and coloured by the engine's call.  A
well that sits apart from its group -- a wrong call, a misplaced standard,
an odd product -- stands out at once.  Clicking a trace shows that well in
the main viewer.

:func:`overlay_data` is the headless part (tested); :class:`PlateOverlay` is
the Tk window, given the main window's loaders as callbacks so it does not
depend on it.
"""
from __future__ import annotations

from typing import Callable, Dict, Optional

import numpy as np

import snp_genotyping as sg

METHODS = (("piecewise", "Piecewise through the bands"),
           ("affine", "Shift + stretch (affine)"),
           ("scale", "Sum of band times (2005 CEQ)"),
           ("none", "Not aligned (as recorded)"))
CALL_COLOURS = {"hom1": "#3B82C4", "hom2": "#E08A2E", "het": "#3BA55C",
                "no-call": "#9AA0A6"}
MARGIN = 150            # scans shown either side of the standard


def _col_for(note: str, role: str, base_order: str) -> Optional[int]:
    """acgt column of the 'sample' or 'standard' channel named in a plate note."""
    key = "sample on Ch" if role == "sample" else "standard on Ch"
    i = note.find(key)
    if i < 0:
        return None
    try:
        ch = int(note[i + len(key)])
    except (ValueError, IndexError):
        return None
    return sg.acgt_index_for_channel(base_order, ch)


def overlay_data(docs: Dict[str, object], base_order: str = "ACTG",
                 method: str = "piecewise", channel: str = "sample",
                 rows: Optional[Dict[str, dict]] = None) -> dict:
    """Everything the overlay draws, for one plate.

    *docs* maps well -> trace document (``acgt``) or array.  *rows* are the
    :func:`snp_genotyping.genotype_plate` rows if already made.  Returns
    ``{"traces": {well: (x, y)}, "calls": {well: call}, "ref": [bands],
    "note": plate note, "rows": rows, "xlim": (lo, hi)}``; ``traces`` is empty
    when no standard was found."""
    if rows is None:
        rows = sg.genotype_plate(docs, base_order=base_order)
    note = next(iter(rows.values()))["plate_note"] if rows else ""
    arrs = {w: np.asarray(getattr(d, "acgt", d), float) for w, d in docs.items()}
    col = _col_for(note, channel, base_order)
    std = sg.std_from_rows(rows)
    ref = sg.plate_reference(std)
    out = {"traces": {}, "calls": {w: r["call"] for w, r in rows.items()},
           "ref": ref or [], "note": note, "rows": rows, "xlim": None}
    if col is None or not ref:
        return out
    if method == "none":
        aligned = {w: a for w, a in arrs.items() if std.get(w)}
        allstd = [s for v in std.values() for s in v]
        lo, hi = min(allstd) - MARGIN, max(allstd) + MARGIN
    else:
        aligned, ref = sg.align_plate(arrs, rows, method, ref=ref)
        lo, hi = ref[0] - MARGIN, ref[-1] + MARGIN
    lo, hi = int(max(0, lo)), int(hi)
    for w, a in aligned.items():
        y = sg.flatten(a[:, col])
        x = np.arange(len(y))
        keep = (x >= lo) & (x < hi)
        out["traces"][w] = (x[keep], y[keep])
    out["xlim"] = (lo, hi)
    return out


class PlateOverlay:
    """The Tk window.  Callbacks from the main window:

    *runs()* -> run names; *files_for(run)* -> ``{well: path}``;
    *load(path)* -> trace document; *order_for(doc)* -> dye order;
    *show(paths)* shows wells in the main viewer; *status(text)*."""

    def __init__(self, master, runs: Callable, files_for: Callable, load: Callable,
                 order_for: Callable, show: Callable, status: Callable = print):
        import tkinter as tk
        from tkinter import ttk
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        from matplotlib.figure import Figure

        self.runs, self.files_for, self.load = runs, files_for, load
        self.order_for, self.show, self.status = order_for, show, status
        self._rows: Dict[str, dict] = {}          # run -> engine rows (method-independent)
        self._docs: Dict[str, dict] = {}
        self._paths: Dict[str, dict] = {}
        self.win = tk.Toplevel(master)
        self.win.title("Plate overlay")
        self.win.transient(master)
        bar = ttk.Frame(self.win, padding=(8, 8, 8, 2))
        bar.pack(fill=tk.X)
        ttk.Label(bar, text="Run").pack(side=tk.LEFT)
        names = list(runs())
        self.run_var = tk.StringVar(value=names[0] if names else "")
        cb = ttk.Combobox(bar, textvariable=self.run_var, values=names, state="readonly", width=36)
        cb.pack(side=tk.LEFT, padx=(4, 10))
        cb.bind("<<ComboboxSelected>>", lambda e: self.refresh())
        ttk.Label(bar, text="Alignment").pack(side=tk.LEFT)
        self.method_var = tk.StringVar(value=METHODS[0][1])
        mb = ttk.Combobox(bar, textvariable=self.method_var, values=[m[1] for m in METHODS],
                          state="readonly", width=28)
        mb.pack(side=tk.LEFT, padx=(4, 10))
        mb.bind("<<ComboboxSelected>>", lambda e: self.refresh())
        ttk.Label(bar, text="Channel").pack(side=tk.LEFT)
        self.chan_var = tk.StringVar(value="sample")
        for v, t in (("sample", "Sample"), ("standard", "Standard")):
            ttk.Radiobutton(bar, text=t, value=v, variable=self.chan_var,
                            command=self.refresh).pack(side=tk.LEFT, padx=2)
        self.fig = Figure(figsize=(10, 4.8), dpi=96)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.win)
        tb = ttk.Frame(self.win)
        tb.pack(fill=tk.X, padx=8)
        NavigationToolbar2Tk(self.canvas, tb).update()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=8, pady=2)
        self.canvas.mpl_connect("pick_event", self._picked)
        self.note = ttk.Label(self.win, foreground="#555", padding=(8, 0, 8, 2), justify=tk.LEFT)
        self.note.pack(fill=tk.X)
        hint = ttk.Label(self.win, foreground="#555", padding=(8, 0, 8, 8), text=(
            "Every well on the plate's own time axis, aligned on its internal standard; "
            "colour = the engine's call. Click a trace to show that well."))
        hint.pack(fill=tk.X)
        for lbl in (self.note, hint):
            lbl.bind("<Configure>", lambda e, l=lbl: l.configure(wraplength=max(200, e.width - 16)))
        self.win.geometry("1000x620")
        self.refresh()

    def exists(self) -> bool:
        try:
            return bool(self.win.winfo_exists())
        except Exception:                                  # noqa: BLE001
            return False

    def lift(self):
        self.win.deiconify()
        self.win.lift()

    def _method(self) -> str:
        return next((k for k, t in METHODS if t == self.method_var.get()), "piecewise")

    def _plate(self, run):
        if run not in self._docs:
            paths = self.files_for(run)
            docs = {}
            for i, (w, p) in enumerate(sorted(paths.items())):
                self.status(f"Plate overlay: loading {run} {w} ({i + 1}/{len(paths)})…")
                self.win.update_idletasks()
                try:
                    docs[w] = self.load(p)
                except Exception:                          # noqa: BLE001
                    continue
            self._docs[run], self._paths[run] = docs, paths
        return self._docs[run]

    def refresh(self):
        run = self.run_var.get()
        if not run:
            return
        self.win.config(cursor="watch")
        self.win.update_idletasks()
        try:
            docs = self._plate(run)
            if not docs:
                self._draw_empty(f"No traces found for {run}.")
                return
            order = self.order_for(next(iter(docs.values())))
            if run not in self._rows:
                self.status(f"Plate overlay: finding the standard on {run}…")
                self.win.update_idletasks()
                self._rows[run] = sg.genotype_plate(docs, base_order=order, run_name=run)
            data = overlay_data(docs, order, self._method(), self.chan_var.get(),
                                rows=self._rows[run])
            self._draw(run, data)
        finally:
            self.win.config(cursor="")
            self.status("")

    def _draw_empty(self, text):
        self.ax.clear()
        self.ax.text(0.5, 0.5, text, ha="center", va="center", transform=self.ax.transAxes)
        self.canvas.draw_idle()
        self.note.config(text="")

    def _draw(self, run, data):
        self.ax.clear()
        if not data["traces"]:
            self._draw_empty("No internal standard found on this plate, so its wells "
                             "cannot be aligned.")
            self.note.config(text=data["note"])
            return
        counts = {}
        order = ("no-call", "hom1", "hom2", "het")     # no-calls underneath
        for call in order:
            for w, (x, y) in data["traces"].items():
                c = data["calls"].get(w, "no-call")
                if c != call:
                    continue
                counts[c] = counts.get(c, 0) + 1
                ln, = self.ax.plot(x, y, color=CALL_COLOURS.get(c, "#8E6CC7"), lw=0.7,
                                   alpha=0.45 if c == "no-call" else 0.75, picker=True,
                                   pickradius=4)
                ln.set_gid(w)
        if self._method() != "none":
            for r in data["ref"]:
                self.ax.axvline(r, color="k", ls=":", lw=0.8)
        self.ax.set_xlim(*data["xlim"])
        self.ax.set_ylim(bottom=0)
        self.ax.set_xlabel("scan on the plate's axis" if self._method() != "none" else "scan")
        from matplotlib.lines import Line2D
        handles = [Line2D([], [], color=CALL_COLOURS[c], lw=2) for c in order if c in counts]
        labels = [f"{c} ({counts[c]})" for c in order if c in counts]
        self.ax.legend(handles, labels, loc="upper right", fontsize=8)
        self.ax.set_title(f"{run}: {len(data['traces'])} wells, "
                          f"{self.chan_var.get()} channel", fontsize=10)
        self.fig.tight_layout()
        self.canvas.draw_idle()
        self.note.config(text=data["note"])

    def _picked(self, event):
        w = event.artist.get_gid()
        path = self._paths.get(self.run_var.get(), {}).get(w)
        if path is not None:
            self.show([path])
