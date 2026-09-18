#!/usr/bin/env python3
"""
Limoncello CE Analyzer (Python)
===============================
A desktop viewer and base caller for capillary-electrophoresis traces:

  • Multi-folder project load (RSD / SCF / ABI .ab1 / text traces)
  • Multi-instrument trace display with current (µA) overlay
  • 1–8 electropherogram panes, plus a Wrap view (one well in N rows)
  • Raw / processed / wrapped / base-called views
  • Quality profile (0–100, teal=good / brown=poor) + base letters
  • Basecaller versions (pos_bonus07, pos_profile, hz_soften, raw_peaks)
  • Auto-tour (slideshow) through the loaded wells
  • Sort wells by well / name / folder, jump-to-well filter
  • Export: FASTA, peak CSV, trace text (channels V + current µA)

Base calling uses our own tuned caller (best_basecaller spacing tracker),
not commercial instrument algorithms.

Genotyping (allele calling next to base calling) is planned for a later
release; the toolbar/menu layout already leaves room for it.

Run:
  python sequence_analyzer.py
  python sequence_analyzer.py --folder /path/to/rsd_dir
"""
from __future__ import annotations

import argparse
import json
import random
import re
import shutil
import subprocess
import sys
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from pathlib import Path
from typing import List, Optional

import numpy as np

# matplotlib embedded
import matplotlib

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter

from analyzer_core import (
    AnalysisSettings,
    BASECALLER_VERSIONS,
    SOURCE_LABELS,
    SUPPORTED_EXTS,
    TraceDocument,
    discover_files,
    find_run_folders,
    load_trace,
    run_basecall,
    display_trace,
)

CHANNEL_ORDER = "ACGT"
QUALITY_COLOR = "#C77B00"   # brownish (poor-quality colour, manual)
QUALITY_GOOD = "#00998A"    # blue-green (good-quality colour, manual)
CURRENT_COLOR = "#8C8C9E"

# Limoncello palette — the lime-yellow, slightly fluorescent colour of the
# Italian lemon liqueur that gives this app its name.
LIMON_BAR = "#D9F24A"       # fluorescent lime-yellow top strip
LIMON_DRK = "#22350E"       # dark green-lemon for text on the bar
LIMON_ACCENT = "#BEE027"    # accent green-yellow
DOC_CACHE_MAX = 300         # keep memory bounded on huge run collections
SCAN_RATE_HZ = 1.75         # instrument scan rate (~0.571 s per scan)
WELL_RE = re.compile(r"([A-Ha-h]\d{1,2})")

TRACE_THEMES = {
    "Classic": {"A": "#00AA00", "C": "#0000DD", "G": "#111111", "T": "#DD0000"},
    "Chromas": {"A": "#00AA00", "C": "#1E90FF", "G": "#444444", "T": "#FF0000"},
    "High-contrast": {"A": "#2E8B57", "C": "#1F4FC0", "G": "#000000", "T": "#E03030"},
    "Monochrome": {"A": "#777777", "C": "#555555", "G": "#333333", "T": "#888888"},
}

SRC_FG = {
    "rsd": "#000000",
    "scf": "#1E55C0",
    "abi": "#0F7A4D",
    "text": "#B45A00",
}

# well sorting modes: "letter first" = A01,A02,…; "number first" = A01,B01,…H01,A02…
SORT_OPTIONS = [
    ("well_row", "Well (letter first)"),
    ("well_col", "Well (number first)"),
    ("name", "File name"),
    ("folder", "Folder + name"),
]
SORT_CODES = {label: code for code, label in SORT_OPTIONS}
SORT_LABELS = {code: label for code, label in SORT_OPTIONS}


class ToolTip:
    """Small floating label that appears after hovering a widget.
    `getter` is called on hover and returns the text to show (or '')."""

    def __init__(self, widget, getter):
        self.widget = widget
        self.getter = getter
        self.tip = None
        self._after = None
        widget.bind("<Enter>", lambda e: self._schedule(), add="+")
        widget.bind("<Leave>", lambda e: self._cancel(), add="+")
        widget.bind("<Motion>", lambda e: self._schedule(), add="+")

    def _schedule(self):
        # Never stack tips: drop any visible one and re-arm.
        self._cancel()
        self._after = self.widget.after(350, self._show)

    def _cancel(self):
        if self._after:
            self.widget.after_cancel(self._after)
            self._after = None
        if self.tip is not None:
            self.tip.destroy()
            self.tip = None

    def _show(self):
        self._after = None
        text = self.getter()
        if not text:
            return
        if self.tip is not None:
            self.tip.destroy()
            self.tip = None
        self.tip = tk.Toplevel(self.widget, bg="#FFFFE1")
        self.tip.wm_overrideredirect(True)
        self.tip.wm_geometry(f"+{self.widget.winfo_pointerx() + 14}"
                             f"+{self.widget.winfo_pointery() + 14}")
        tk.Label(self.tip, text=text, bg="#FFFFE1", fg="#222222",
                 justify=tk.LEFT, font=("DejaVu Sans", 9),
                 padx=4, pady=2).pack()


class LimoncelloAnalyzerApp(tk.Tk):
    def __init__(self, initial_folders: Optional[List[Path]] = None):
        super().__init__()
        self.title("Limoncello CE Analyzer")
        self.minsize(1040, 720)
        self.configure(bg="#F2F4F7")
        # Undecorated but STILL WM-managed: 'splash' removes the gray OS
        # titlebar while mutter keeps managing the window (focus, placement).
        # Unlike override-redirect, this does not wedge the desktop.
        try:
            self.attributes("-type", "splash")
        except tk.TclError:
            pass
        self.geometry("1440x920")
        self._drag_off = (0, 0)
        self._zoomed = False
        self._pre_zoom = ""
        self._busy = False
        self._cancel = False
        self.bind("<Escape>", self._on_escape)
        self.bind("<F11>", lambda e: self._toggle_zoom())
        self.bind("<F1>", lambda e: self.show_help())
        self.bind("<Home>", lambda e: self._reset_zoom())

        self.folders: List[Path] = list(initial_folders or [])
        self.files: List[Path] = []
        self.docs: dict[str, TraceDocument] = {}
        self.selected: List[Path] = []
        self.excluded: set[str] = set()   # individual files hidden by the user
        self._bg_checked = False
        self._bg_imgs: List[np.ndarray] = []
        self.settings = AnalysisSettings()
        self.n_graphs = tk.IntVar(value=1)
        self.view_mode = tk.StringVar(value="processed")
        self.basecaller = tk.StringVar(value="pos_bonus07")
        self.status_var = tk.StringVar(value="Ready — add a data folder to begin")

        # Visualization state
        self.chan_show = [tk.BooleanVar(value=True) for _ in CHANNEL_ORDER]
        self.show_letters = tk.BooleanVar(value=True)
        self.show_qnum = tk.BooleanVar(value=False)
        self.show_qcurve = tk.BooleanVar(value=True)
        self.show_current = tk.BooleanVar(value=True)
        self.x_time = tk.BooleanVar(value=False)  # x axis: scans (off) vs time (on)
        self.theme = tk.StringVar(value="Classic")
        self.wrap_rows = tk.IntVar(value=5)

        # Shared X/Y view, stored as [first, last] fractions of the full data
        # range. One model drives every visible graph plus the axis bars.
        self._view_x = [0.0, 1.0]
        self._view_y = [0.0, 1.0]
        self._plot_axes: List = []
        self._full_xlim = None
        self._full_ylim = None

        # List sorting / filtering / tour
        self.sort_mode = tk.StringVar(value="well_row")
        self.filter_var = tk.StringVar()
        self.cycle_ms = tk.IntVar(value=3000)
        self.touring = False
        self.cycle_pos = 0
        self._tour_job = None

        self._style()
        self._build_menu()
        self._build_layout()
        roots = list(self.folders)
        self.folders = []
        for root in roots:
            self._absorb_runs(root)
        self.refresh_file_list()
        # Draw the empty (background) view once the window is up, otherwise
        # the plot area stays blank until the first selection.
        self.after(50, self.redraw)

    # ------------------------------------------------------------------ style
    def _style(self):
        style = ttk.Style(self)
        # NOTE: stay on the native theme — the 'clam' theme repaints the
        # scale thumb on every hover motion and flickers on X11.
        bg = "#F2F4F7"
        style.configure("TFrame", background=bg)
        style.configure("TLabel", background=bg)
        style.configure("TLabelframe", background=bg)
        style.configure("TLabelframe.Label", background=bg, foreground="#0F3A6E",
                        font=("", 9, "bold"))
        style.configure("TListbox", background="#FFFFFF", fieldbackground="#FFFFFF",
                        foreground="#1A1C22")
        style.configure("TCheckbutton", background=bg)

    # ------------------------------------------------------ custom titlebar
    def _start_move(self, event):
        self._drag_off = (event.x_root - self.winfo_x(),
                          event.y_root - self.winfo_y())

    def _on_move(self, event):
        x = event.x_root - self._drag_off[0]
        y = event.y_root - self._drag_off[1]
        sw, sh = self.winfo_screenwidth(), self.winfo_screenheight()
        # keep at least a sliver on-screen so the bar can always be grabbed
        x = max(200 - self.winfo_width(), min(x, sw - 200))
        y = max(0, min(y, sh - 60))
        self.geometry(f"+{x}+{y}")

    def _toggle_zoom(self):
        if self._zoomed:
            self.geometry(self._pre_zoom)
            self._zoomed = False
        else:
            self._pre_zoom = self.geometry()
            sw, sh = self.winfo_screenwidth(), self.winfo_screenheight()
            self.geometry(f"{sw}x{sh}+0+0")
            self._zoomed = True

    def _close(self):
        self.destroy()

    def _popup_menu(self, menu, mb):
        """Open a top menu bar dropdown (stock Menubutton posting fails here)."""
        try:
            menu.tk_popup(mb.winfo_rootx(), mb.winfo_rooty() + mb.winfo_height())
        finally:
            menu.grab_release()
        return "break"

    def _on_escape(self, event=None):
        if self._busy:
            self._cancel = True
            self.status_var.set("Cancelling…")
        elif self._zoomed:
            self._toggle_zoom()

    # ------------------------------------------------------------------ UI
    def _build_menu(self):
        file_m = tk.Menu(self, tearoff=0)
        file_m.add_command(label="Add data folder…", command=self.add_folder)
        file_m.add_command(label="Clear folders", command=self.clear_folders)
        file_m.add_separator()
        file_m.add_command(label="Export sequence (FASTA)…", command=self.export_fasta)
        file_m.add_command(label="Export peak table (CSV)…", command=self.export_peaks)
        file_m.add_command(label="Export trace text (V + µA)…", command=self.export_text)
        file_m.add_separator()
        file_m.add_command(label="Save graph image…", command=self.save_figure)
        file_m.add_separator()
        file_m.add_command(label="Save settings JSON…", command=self.save_settings)
        file_m.add_command(label="Load settings JSON…", command=self.load_settings)
        file_m.add_separator()
        file_m.add_command(label="Exit", command=self.destroy)

        view_m = tk.Menu(self, tearoff=0)
        for mode, label in [
            ("raw", "Raw traces"),
            ("processed", "Processed (ACGT)"),
            ("called", "Base-called (with peaks)"),
            ("wrap", "Wrap (one well, N rows)"),
        ]:
            view_m.add_radiobutton(
                label=label, variable=self.view_mode, value=mode, command=self.redraw)
        view_m.add_separator()
        view_m.add_checkbutton(
            label="Show base letters on peaks", variable=self.show_letters,
            command=self.redraw)
        view_m.add_checkbutton(
            label="Show quality numbers (zoomed)", variable=self.show_qnum,
            command=self.redraw)
        view_m.add_checkbutton(
            label="Show quality profile (0-100)", variable=self.show_qcurve,
            command=self.redraw)
        view_m.add_checkbutton(
            label="Show instrument current (µA)", variable=self.show_current,
            command=lambda: (self._sync_current_btn(), self.redraw()))
        view_m.add_separator()

        sort_m = tk.Menu(self, tearoff=0)
        for code, label in SORT_OPTIONS:
            sort_m.add_radiobutton(label=label, variable=self.sort_mode, value=code,
                                   command=self.refresh_file_list)
        view_m.add_cascade(label="Sort files", menu=sort_m)

        wrap_m = tk.Menu(self, tearoff=0)
        wrap_var = tk.StringVar(value=str(self.wrap_rows.get()))
        for r in (2, 3, 4, 5, 6, 8):
            wrap_m.add_radiobutton(
                label=f"{r} rows", variable=wrap_var, value=str(r),
                command=lambda r=r: (self.wrap_rows.set(r), self.redraw()))
        view_m.add_cascade(label="Wrap rows", menu=wrap_m)

        tour_m = tk.Menu(self, tearoff=0)
        for ms, label in [(2000, "2 s"), (3000, "3 s"), (4000, "4 s"), (5000, "5 s")]:
            tour_m.add_radiobutton(label=label, value=ms, variable=self.cycle_ms)
        view_m.add_cascade(label="Tour interval", menu=tour_m)

        view_m.add_separator()
        view_m.add_command(label="Page backward", command=lambda: self._page_by(-1))
        view_m.add_command(label="Page forward", command=lambda: self._page_by(1))
        view_m.add_command(label="Start/stop auto-tour", command=self.toggle_tour)
        view_m.add_command(label="Reset view (X & Y)", command=self._reset_zoom)
        view_m.add_command(label="Redraw", command=self.redraw)

        chan_m = tk.Menu(self, tearoff=0)
        for i, base in enumerate(CHANNEL_ORDER):
            chan_m.add_checkbutton(
                label=f"Channel {base}", variable=self.chan_show[i],
                command=self.redraw)
        view_m.add_cascade(label="Channels", menu=chan_m)

        theme_m = tk.Menu(self, tearoff=0)
        for name in TRACE_THEMES:
            theme_m.add_radiobutton(
                label=name, variable=self.theme, value=name, command=self.redraw)
        view_m.add_cascade(label="Trace colors", menu=theme_m)

        analysis_m = tk.Menu(self, tearoff=0)
        analysis_m.add_command(label="Basecall selected", command=self.basecall_selected)
        analysis_m.add_command(label="Basecall all in list", command=self.basecall_all)
        analysis_m.add_separator()
        analysis_m.add_command(label="⚙ Basecall settings…", command=self.open_settings)
        analysis_m.add_separator()
        analysis_m.add_command(label="Genotyping (planned)", state=tk.DISABLED)

        help_m = tk.Menu(self, tearoff=0)
        help_m.add_command(label="User manual…", command=self.show_help)
        help_m.add_command(label="About", command=self.show_about)

        self._menus = {"File": file_m, "View": view_m,
                       "Analysis": analysis_m, "Help": help_m}
    def _build_layout(self):
        # Yellow title bar — the window is undecorated ('splash'), so this is
        # the top of the app and carries the name (centred) + window buttons.
        tbar = tk.Frame(self, bg=LIMON_BAR, bd=0)
        tbar.pack(side=tk.TOP, fill=tk.X)
        tbar.columnconfigure(0, weight=1)
        tbar.columnconfigure(2, weight=1)

        left_pad = tk.Frame(tbar, bg=LIMON_BAR)
        left_pad.grid(row=0, column=0, sticky="nsew")

        mid = tk.Frame(tbar, bg=LIMON_BAR)
        mid.grid(row=0, column=1, pady=4)
        title_lbl = ttk.Label(mid, text="Limoncello CE Analyzer",
                              font=("DejaVu Sans", 16, "bold"),
                              background=LIMON_BAR, foreground=LIMON_DRK)
        title_lbl.pack()
        sub_lbl = ttk.Label(mid, text="multi-instrument trace viewer",
                            background=LIMON_BAR, foreground="#4A5A18")
        sub_lbl.pack()

        btn_style = dict(relief=tk.FLAT, bd=0, font=("DejaVu Sans", 12),
                         bg=LIMON_BAR, fg=LIMON_DRK,
                         activebackground=LIMON_ACCENT, activeforeground=LIMON_DRK,
                         takefocus=0, cursor="hand2", padx=10, pady=3)
        btns = tk.Frame(tbar, bg=LIMON_BAR)
        btns.grid(row=0, column=2, sticky="e")
        tk.Button(btns, text="✕", command=self._close, **btn_style).pack(side=tk.RIGHT)
        tk.Button(btns, text="□", command=self._toggle_zoom, **btn_style).pack(side=tk.RIGHT)

        for widget in (tbar, left_pad, mid, title_lbl, sub_lbl, btns):
            widget.bind("<ButtonPress-1>", self._start_move)
            widget.bind("<B1-Motion>", self._on_move)
            widget.bind("<Double-Button-1>", lambda e: self._toggle_zoom())
        acc = tk.Frame(self, bg=LIMON_ACCENT, height=2)
        acc.pack(side=tk.TOP, fill=tk.X)
        acc.pack_propagate(False)

        # Menu bar (gray strip below the yellow banner)
        mbar = tk.Frame(self, bg="#E8EBEF", bd=0, highlightthickness=0)
        mbar.pack(side=tk.TOP, fill=tk.X)
        for label, menu in self._menus.items():
            mb = tk.Menubutton(
                mbar, text=f"  {label}  ", menu=menu, relief=tk.FLAT,
                font=("DejaVu Sans", 10), fg="#1A1C22", bg="#E8EBEF",
                activebackground="#D6DCE3", activeforeground="#0F3A6E",
                bd=0, padx=0, pady=2, takefocus=0, cursor="hand2")
            mb.pack(side=tk.LEFT, padx=2, pady=1)
            # Stock Menubutton posting misbehaves under this WM, and a menu
            # opened on press is closed again by the same click's release.
            # Open on release instead, with the default handling suppressed.
            mb.bind("<ButtonPress-1>", lambda e: "break")
            mb.bind("<ButtonRelease-1>",
                    lambda e, m=menu, b=mb: self._popup_menu(m, b))

        # Status bar
        ttk.Label(self, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W).pack(
            side=tk.BOTTOM, fill=tk.X
        )

        body = ttk.Panedwindow(self, orient=tk.HORIZONTAL)
        body.pack(fill=tk.BOTH, expand=True)

        # ----- Left: folders + files -----
        left = ttk.Frame(body, width=380)
        body.add(left, weight=1)

        ttk.Label(left, text="⛁ Data folders", font=("", 10, "bold"),
                  foreground="#0F3A6E").pack(anchor=tk.W, padx=4, pady=(4, 0))
        self.folder_list = tk.Listbox(left, height=5, exportselection=False,
                                      selectmode=tk.EXTENDED,
                                      bg="#FFFFFF", highlightthickness=1,
                                      highlightbackground="#B9C2D2",
                                      activestyle="none")
        self.folder_list.pack(fill=tk.X, padx=4, pady=2)
        self.folder_list.bind("<<ListboxSelect>>", self._on_folder_select)
        self.folder_list.bind("<Double-Button-1>", self._on_folder_activate)
        fscroll = ttk.Scrollbar(left, orient=tk.HORIZONTAL,
                                command=self.folder_list.xview)
        fscroll.pack(fill=tk.X, padx=4)
        self.folder_list.configure(xscrollcommand=fscroll.set)
        self._folder_paths: List[Path] = []
        ToolTip(self.folder_list, self._folder_tip_text)
        fb = ttk.Frame(left)
        fb.pack(fill=tk.X, padx=4)
        ttk.Button(fb, text="⊕ Add…", command=self.add_folder).pack(side=tk.LEFT)
        ttk.Button(fb, text="– Remove", command=self.remove_folder).pack(side=tk.LEFT, padx=4)

        self.samples_label = ttk.Label(left, text="◧ Samples / files",
                                       font=("", 10, "bold"),
                                       foreground="#0F3A6E")
        self.samples_label.pack(anchor=tk.W, padx=4, pady=(8, 0))

        # sort + filter row
        srow = ttk.Frame(left)
        srow.pack(fill=tk.X, padx=4)
        ttk.Label(srow, text="Sort:").pack(side=tk.LEFT)
        self.sort_labels = {c: l for c, l in SORT_OPTIONS}
        self.sort_codes = {l: c for c, l in SORT_OPTIONS}
        self.sort_label_var = tk.StringVar(value=self.sort_labels["well_row"])
        self.sort_cb = ttk.Combobox(srow, textvariable=self.sort_label_var, width=18,
                                    values=list(SORT_LABELS.values()), state="readonly")
        self.sort_cb.pack(side=tk.LEFT, padx=(2, 4))
        self.sort_cb.bind("<<ComboboxSelected>>",
                          lambda e: (self.sort_mode.set(self.sort_codes[self.sort_label_var.get()]),
                                     self.refresh_file_list()))
        ttk.Label(srow, text="Jump/well:").pack(side=tk.LEFT)
        self.well_entry = ttk.Entry(srow, textvariable=self.filter_var, width=6)
        self.well_entry.pack(side=tk.LEFT, padx=2)
        self.filter_var.trace_add("write", lambda *a: self.refresh_file_list())

        self.file_list = tk.Listbox(left, selectmode=tk.EXTENDED, exportselection=False,
                                    bg=LIMON_BAR, fg=LIMON_DRK,
                                    selectbackground=LIMON_ACCENT, selectforeground=LIMON_DRK,
                                    highlightthickness=1,
                                    highlightbackground="#B9C2D2",
                                    activestyle="none")
        self.file_list.pack(fill=tk.BOTH, expand=True, padx=4, pady=2)
        self.file_list.bind("<<ListboxSelect>>", self.on_file_select)
        self.file_list.bind("<space>", lambda e: self.toggle_tour())
        self.file_list.bind("<Down>", self._page_key)
        self.file_list.bind("<Up>", self._page_key)
        self.file_list.bind("<Page_Down>", self._page_key)
        self.file_list.bind("<Page_Up>", self._page_key)
        self.file_list.bind("<Delete>", lambda e: self.remove_samples())
        self.bind("<Down>", self._page_key_root)
        self.bind("<Up>", self._page_key_root)
        self.bind("<Page_Down>", self._page_key_root)
        self.bind("<Page_Up>", self._page_key_root)

        rrow = ttk.Frame(left)
        rrow.pack(fill=tk.X, padx=4, pady=(0, 2))
        ttk.Button(rrow, text="– Remove file(s)",
                   command=self.remove_samples).pack(side=tk.LEFT)

        frow = ttk.Frame(left)
        frow.pack(fill=tk.X, padx=4, pady=(0, 2))
        ttk.Label(frow, text="Graphs").pack(side=tk.LEFT)
        tk.Spinbox(
            frow, from_=1, to=8, textvariable=self.n_graphs, width=3,
            command=self.redraw, highlightthickness=0
        ).pack(side=tk.LEFT, padx=2)
        ttk.Label(frow, text="Tour (ms)").pack(side=tk.LEFT, padx=(8, 0))
        tk.Spinbox(
            frow, from_=500, to=20000, increment=500, textvariable=self.cycle_ms,
            width=5, highlightthickness=0
        ).pack(side=tk.LEFT, padx=2)
        ttk.Button(frow, text="▶", width=3, command=self.toggle_tour).pack(side=tk.LEFT, padx=2)

        # ----- Center: plots -----
        center = ttk.Frame(body)
        body.add(center, weight=5)

        chan_bar = ttk.Frame(center)
        chan_bar.pack(fill=tk.X, padx=4, pady=(4, 0))
        ttk.Label(chan_bar, text="Channels:").pack(side=tk.LEFT)
        for i, base in enumerate(CHANNEL_ORDER):
            cb = tk.Checkbutton(
                chan_bar, text=f"  {base}  ", variable=self.chan_show[i],
                command=self.redraw, bg="#F2F4F7",
                fg=TRACE_THEMES["Classic"][base], activebackground="#F2F4F7",
                selectcolor="white")
            cb.pack(side=tk.LEFT, padx=2)
        self._cur_btn = tk.Checkbutton(
            chan_bar, text="  µA  ", variable=self.show_current,
            bg="#DFF6DC", activebackground="#DFF6DC", fg="#14532D",
            selectcolor="white",
            command=lambda: (self._sync_current_btn(), self.redraw()))
        self._cur_btn.pack(side=tk.LEFT, padx=(8, 2))
        self._sync_current_btn()
        ttk.Label(chan_bar, text="   Signal: Volts").pack(side=tk.LEFT, padx=4)
        self._time_btn = tk.Checkbutton(
            chan_bar, text="  Time (s)  ", variable=self.x_time,
            bg="#F2F4F7", activebackground="#F2F4F7", selectcolor="white",
            command=self.redraw)
        self._time_btn.pack(side=tk.LEFT, padx=(8, 2))
        ttk.Button(chan_bar, text="⟲ Reset view",
                   command=self._reset_zoom).pack(side=tk.RIGHT, padx=4)

        # Axis bars: drag to pan, mouse-wheel to zoom. One bar per axis and it
        # drives every visible graph at once (shared X / shared Y).
        plotf = ttk.Frame(center)
        plotf.pack(fill=tk.BOTH, expand=True)
        plotf.rowconfigure(0, weight=1)
        plotf.columnconfigure(0, weight=1)

        self.fig = Figure(figsize=(9, 7), dpi=100, facecolor="#FFFFFF")
        self.canvas = FigureCanvasTkAgg(self.fig, master=plotf)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        self.canvas.get_tk_widget().bind("<Button-3>", self._on_right_click)
        self.xbar = ttk.Scrollbar(plotf, orient=tk.HORIZONTAL)
        self.xbar.grid(row=1, column=0, sticky="ew")
        self.ybar = ttk.Scrollbar(plotf, orient=tk.VERTICAL)
        self.ybar.grid(row=0, column=1, sticky="ns")
        self.xbar.configure(command=lambda *a: self._bar_cmd("x", a))
        self.ybar.configure(command=lambda *a: self._bar_cmd("y", a))
        for bar, axis in ((self.xbar, "x"), (self.ybar, "y")):
            bar.bind("<Button-4>", lambda e, ax=axis: self._wheel_zoom(ax, 1))
            bar.bind("<Button-5>", lambda e, ax=axis: self._wheel_zoom(ax, -1))
            bar.bind("<MouseWheel>",
                     lambda e, ax=axis: self._wheel_zoom(ax, 1 if e.delta > 0 else -1))
            bar.bind("<Double-Button-1>", lambda e, ax=axis: self._reset_axis(ax))
            bar.bind("<Button-3>", self._on_right_click)
            bar.configure(cursor="hand2")

        # Sequence readout
        hdr = ttk.Frame(center)
        hdr.pack(fill=tk.X, padx=4)
        ttk.Label(hdr, text="Called sequence").pack(side=tk.LEFT)
        self.seq_text = scrolledtext.ScrolledText(center, height=4, wrap=tk.CHAR,
                                                  font=("Courier", 9))
        self.seq_text.pack(fill=tk.X, padx=4, pady=2)

        # ----- Advanced base-call parameters live in a menu dialog ----
        # (no right-hand knob panel → more room for plots, no flicker)
        self._build_settings_dialog()

    def _build_settings_dialog(self):
        """Advanced base-call parameters in a separate (hidden) window."""
        win = tk.Toplevel(self)
        win.title("⚙ Basecall settings (advanced)")
        win.geometry("380x620")
        win.minsize(340, 420)
        win.withdraw()
        self.settings_win = win

        frm = ttk.Frame(win)
        frm.pack(fill=tk.BOTH, expand=True)
        canvas = tk.Canvas(frm, highlightthickness=0, bg="#F2F4F7")
        sb = ttk.Scrollbar(frm, orient=tk.VERTICAL, command=canvas.yview)
        inner = ttk.Frame(canvas)
        inner.bind("<Configure>",
                   lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=inner, anchor=tk.NW)
        canvas.configure(yscrollcommand=sb.set, bg="#F2F4F7")
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb.pack(side=tk.RIGHT, fill=tk.Y)

        ttk.Label(inner,
                  text="Tuned once per caller — normally you won't need to change "
                       "these.\nApply re-runs the selected wells with them.",
                  foreground="#666", wraplength=340).pack(anchor=tk.W, padx=8, pady=6)

        def section(title: str) -> ttk.LabelFrame:
            f = ttk.LabelFrame(inner, text=title, padding=4)
            f.pack(fill=tk.X, padx=6, pady=4)
            return f

        def check(parent, text, var, cmd=None):
            tk.Checkbutton(parent, text=text, variable=var, command=cmd,
                           bg="#F2F4F7", activebackground="#F2F4F7",
                           selectcolor="#FFFFFF").pack(anchor=tk.W)

        def scale(parent, var, lo, hi, res):
            tk.Scale(parent, from_=lo, to=hi, variable=var, orient=tk.HORIZONTAL,
                     resolution=res, showvalue=0, highlightthickness=0,
                     relief=tk.FLAT, length=280, troughcolor="#DCE4EF",
                     bg="#F2F4F7", activebackground="#DCE4EF", sliderrelief=tk.RAISED
                     ).pack(fill=tk.X)

        f = section("Basecaller version")
        for key, desc in BASECALLER_VERSIONS.items():
            tk.Radiobutton(f, text=key, variable=self.basecaller, value=key,
                           bg="#F2F4F7", activebackground="#F2F4F7",
                           selectcolor="#F2F4F7").pack(anchor=tk.W)
            ttk.Label(f, text=f"  {desc}", foreground="#555").pack(anchor=tk.W)

        f = section("Dye / channel")
        self.base_order_var = tk.StringVar(value=self.settings.base_order)
        ttk.Label(f, text="Base order (instrument)").pack(anchor=tk.W)
        ttk.Combobox(f, textvariable=self.base_order_var,
                     values=["TGCA", "ACGT", "GATC", "CTAG"], width=12).pack(anchor=tk.W)
        ttk.Label(f, text="Applies to RSD / text traces; ABI uses its own FWO order.",
                  foreground="#777", wraplength=250).pack(anchor=tk.W)

        f = section("Baseline")
        self.bl_enable = tk.BooleanVar(value=True)
        check(f, "Enable baseline", self.bl_enable)
        self.bl_win = tk.IntVar(value=self.settings.baseline_window)
        ttk.Label(f, text="Window").pack(anchor=tk.W)
        scale(f, self.bl_win, 51, 401, 1)

        f = section("Spectral & mobility")
        self.spec_en = tk.BooleanVar(value=True)
        self.spec_adapt = tk.BooleanVar(value=True)
        self.mob_en = tk.BooleanVar(value=True)
        check(f, "Spectral separation", self.spec_en)
        check(f, "Position-adaptive matrix", self.spec_adapt)
        check(f, "Mobility correction", self.mob_en)

        f = section("Band filter (deconv)")
        self.gauss_en = tk.BooleanVar(value=True)
        self.mp_wiener = tk.BooleanVar(value=False)
        self.gauss_seg = tk.IntVar(value=384)
        self.gauss_reg = tk.DoubleVar(value=0.05)
        check(f, "Gaussian reconstruction", self.gauss_en)
        check(f, "Multi-pass Wiener (2048/1900)", self.mp_wiener)
        ttk.Label(f, text="Segment size").pack(anchor=tk.W)
        scale(f, self.gauss_seg, 128, 1024, 1)
        ttk.Label(f, text="Noise reg").pack(anchor=tk.W)
        scale(f, self.gauss_reg, 0.01, 0.2, 0.001)

        f = section("Spacing tracker")
        self.bonus = tk.DoubleVar(value=0.7)
        self.pullback = tk.DoubleVar(value=0.008)
        self.ema = tk.DoubleVar(value=0.08)
        self.wlo = tk.DoubleVar(value=0.70)
        self.whi = tk.DoubleVar(value=1.30)
        ttk.Label(f, text="Channel peak bonus").pack(anchor=tk.W)
        scale(f, self.bonus, 0.0, 1.5, 0.1)
        ttk.Label(f, text="Pullback weight").pack(anchor=tk.W)
        scale(f, self.pullback, 0.0, 0.03, 0.001)
        ttk.Label(f, text="EMA alpha").pack(anchor=tk.W)
        scale(f, self.ema, 0.02, 0.25, 0.01)
        ttk.Label(f, text="Window frac lo / hi").pack(anchor=tk.W)
        scale(f, self.wlo, 0.5, 0.95, 0.01)
        scale(f, self.whi, 1.05, 1.6, 0.01)

        f = section("Scan region (display)")
        self.scan_start = tk.IntVar(value=0)
        self.scan_end = tk.IntVar(value=0)
        ttk.Label(f, text="Start (0=auto)").pack(anchor=tk.W)
        ttk.Entry(f, textvariable=self.scan_start, width=10).pack(anchor=tk.W)
        ttk.Label(f, text="End (0=full)").pack(anchor=tk.W)
        ttk.Entry(f, textvariable=self.scan_end, width=10).pack(anchor=tk.W)

        bf = ttk.Frame(inner)
        bf.pack(fill=tk.X, padx=6, pady=8)
        ttk.Button(bf, text="Apply + redraw",
                   command=self.apply_knobs_redraw).pack(fill=tk.X, pady=2)
        ttk.Button(bf, text="Basecall selected with these",
                   command=self.basecall_selected).pack(fill=tk.X, pady=2)
        ttk.Button(bf, text="Close", command=win.withdraw).pack(fill=tk.X, pady=(4, 2))

    def open_settings(self):
        self.settings_win.deiconify()
        self.settings_win.lift()

    # ------------------------------------------------------------------ data
    def _absorb_runs(self, root: Path):
        """Expand a root folder into its run folders and add them to the list.
        Returns (n_found, n_added)."""
        runs = find_run_folders(root)
        known = {str(p.resolve()) for p in self.folders}
        added = 0
        for p in runs:
            if str(p.resolve()) in known:
                continue
            self.folders.append(p)
            self._folder_paths.append(p)
            self.folder_list.insert(tk.END, p.name)
            added += 1
        return len(runs), added

    def _pick_folders(self) -> List[Path]:
        """Pick one or more folders. Prefer zenity (a proper picker where you
        can Ctrl-click several folders); fall back to the Tk chooser."""
        start = str(self._folder_paths[-1].parent) if self._folder_paths else ""
        if shutil.which("zenity"):
            cmd = ["zenity", "--file-selection", "--directory", "--multiple",
                   "--separator=\n",
                   "--title=Add data folder(s): select a run, or Ctrl-click several"]
            if start:
                cmd.append(f"--filename={start}/")
            try:
                r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
                if r.returncode != 0:
                    return []
                return [Path(p) for p in r.stdout.splitlines() if p.strip()]
            except Exception:
                pass
        d = filedialog.askdirectory(parent=self,
                                    title="Select a folder with run folders "
                                          "(rsd / scf / ab1 / text)")
        return [Path(d)] if d else []

    def add_folder(self):
        # A folder that holds many runs (e.g. an OY collection) loads every
        # run; selecting a single run folder loads just that run. Several
        # folders can be picked at once.
        roots = self._pick_folders()
        if not roots:
            return
        total_added = 0
        nofiles = []
        for root in roots:
            n_found, added = self._absorb_runs(root)
            total_added += added
            if n_found == 0:
                nofiles.append(root)
        # show everything again (a stale folder selection would hide new runs)
        self.folder_list.selection_clear(0, tk.END)
        self.refresh_file_list()
        if nofiles:
            messagebox.showwarning(
                "No trace files",
                "No trace files found under:\n"
                + "\n".join(str(p) for p in nofiles), parent=self)
        if total_added:
            self.status_var.set(f"Added {total_added} run folder(s)")

    def _folder_tip_text(self) -> str:
        if not self._folder_paths:
            return ""
        i = self.folder_list.nearest(self.winfo_pointery() - self.folder_list.winfo_rooty())
        if not 0 <= i < len(self._folder_paths):
            return ""
        return str(self._folder_paths[i])

    def _on_folder_select(self, _evt=None):
        sel = self.folder_list.curselection()
        if sel:
            idx = sel[0]
            if 0 <= idx < len(self._folder_paths):
                self.status_var.set(f"Folder: {self._folder_paths[idx]}")
        self.refresh_file_list()

    def _on_folder_activate(self, _evt=None):
        # Double-click a folder row: show only that folder's files and open
        # the first one. Ctrl/Shift-click several rows first to show a union.
        self.refresh_file_list()
        if self.files:
            self.file_list.selection_clear(0, tk.END)
            self.file_list.selection_set(0)
            self.file_list.activate(0)
            self.on_file_select(None)

    def _selected_folder_paths(self) -> List[Path]:
        return [self._folder_paths[i] for i in self.folder_list.curselection()
                if 0 <= i < len(self._folder_paths)]

    def remove_folder(self):
        sel = list(self.folder_list.curselection())
        if not sel:
            self.status_var.set("Select a folder/run in the list above first")
            return
        removed = []
        for i in sorted(sel, reverse=True):
            if i < len(self._folder_paths):
                path = self._folder_paths.pop(i)
                self.folder_list.delete(i)
                removed.append(path)
                try:
                    self.folders.remove(path)
                except ValueError:
                    pass
        # any file-level hides under a removed run no longer apply
        for p in removed:
            rp = str(p.resolve())
            self.excluded = {e for e in self.excluded
                             if not e.startswith(rp + "/") and e != rp}
        gone = {str(p.resolve()) for p in removed}
        self.selected = [p for p in self.selected
                         if str(p.resolve()) not in gone]
        self.folder_list.selection_clear(0, tk.END)
        self.refresh_file_list()
        self.redraw()
        self.status_var.set(f"Removed {len(removed)} run folder(s)")

    def remove_samples(self):
        idxs = self.file_list.curselection()
        if not idxs:
            self.status_var.set("Select sample(s) in the list below first")
            return
        n = 0
        for i in idxs:
            if i < len(self.files):
                self.excluded.add(str(self.files[i].resolve()))
                n += 1
        self.file_list.selection_clear(0, tk.END)
        self.selected = []
        self.refresh_file_list()
        self.redraw()
        self.status_var.set(f"Removed {n} sample file(s) from the list")

    def clear_folders(self):
        self.folders.clear()
        self._folder_paths.clear()
        self.excluded.clear()
        self.selected = []
        self.folder_list.delete(0, tk.END)
        self.refresh_file_list()
        self.redraw()

    def _well_parts(self, name: str):
        m = WELL_RE.search(name)
        if not m:
            return "", 9999
        w = m.group(1).upper()
        letter = w[:1] if w[:1] in "ABCDEFGH" else "Z"
        nums = re.findall(r"\d+", w)
        num = int(nums[0]) if nums else 9999
        return letter, num

    def _sort_key(self, f: Path):
        letter, num = self._well_parts(f.stem)
        mode = self.sort_mode.get()
        if mode == "well_col":
            return (num, letter, f.parent.name, f.name)
        if mode == "well_row":
            return (letter, num, f.parent.name, f.name)
        if mode == "name":
            return (f.name.lower(), f.parent.name)
        return (f.parent.name, f.name)

    def refresh_file_list(self):
        selected = self._selected_folder_paths()
        # No folder selected -> every added run (cross-run compare). Selecting
        # one or more folders narrows the sample list to those folders' files.
        folders = selected if selected else self.folders
        files = [f for f in discover_files(folders)
                 if str(f.resolve()) not in self.excluded]
        filt = self.filter_var.get().strip().lower()
        if filt:
            files = [f for f in files
                     if filt in f.stem.lower()
                     or WELL_RE.search(f.stem) and filt in WELL_RE.search(f.stem).group(1).lower()
                     or filt in f.parent.name.lower()]
        self.files = sorted(files, key=self._sort_key)
        self.file_list.delete(0, tk.END)
        if self.files:
            self.file_list.insert(tk.END, *(f"{f.parent.name}/{f.name}"
                                            for f in self.files))
        colorize = len(self.files) <= 5000  # keep huge runs fast
        vis: dict = {}
        for i, f in enumerate(self.files):
            if colorize:
                ext = f.suffix.lower().lstrip(".")
                self.file_list.itemconfig(i, fg=SRC_FG.get(ext, "#000"))
            vis[str(f.resolve())] = i
        # keep a selection if possible
        sel_idxs = []
        for p in self.selected:
            i = vis.get(str(p.resolve()))
            if i is not None:
                sel_idxs.append(i)
        if sel_idxs:
            for i in sel_idxs:
                self.file_list.selection_set(i)
        exts = " ".join(sorted({f".{f.suffix.lower().lstrip('.')}" for f in self.files}))
        scope = (f"{len(selected)} selected folder(s)" if selected
                 else f"{len(self.folders)} folder(s)")
        if hasattr(self, "samples_label"):
            self.samples_label.configure(
                text=f"◧ Samples / files  ({len(self.files)})")
        self.status_var.set(f"{len(self.files)} files ({exts or 'none'}) in {scope}"
                            + ("  [filtered]" if filt else ""))

    def on_file_select(self, _evt=None):
        idxs = self.file_list.curselection()
        self.selected = [self.files[i] for i in idxs if i < len(self.files)]
        if self.selected:
            self.cycle_pos = 0  # manual selection resets any tour position
        self.redraw()

    # ------------------------------------------------------------------ page
    def _current_page_start(self) -> int:
        if self.selected:
            key = str(self.selected[0].resolve())
            for i, p in enumerate(self.files):
                if str(p.resolve()) == key:
                    return i
        return self.cycle_pos

    def _page_by(self, delta: int):
        """Manual paging: advance/retreat the selection window by n graphs."""
        m = len(self.files)
        if m == 0:
            return
        n = max(1, self.n_graphs.get())
        if m <= n:
            idxs = list(range(m))
        else:
            start = self._current_page_start() + delta * n
            if start >= m:
                start = 0
            elif start < 0:
                start = m - n
            if start + n > m:
                start = m - n
            idxs = list(range(start, start + n))
        self.cycle_pos = idxs[0] if idxs else 0
        self.selected = [self.files[i] for i in idxs]
        self.file_list.selection_clear(0, tk.END)
        for i in idxs:
            self.file_list.selection_set(i)
        if idxs:
            self.file_list.see(idxs[-1])
        if idxs:
            self.status_var.set("Page: " + "  ".join(
                f"{p.parent.name}/{p.name}" for p in self.selected))
        self.redraw()

    def _page_key(self, event):
        self._page_by(1 if event.keysym in ("Down", "Page_Down") else -1)
        return "break"

    def _page_key_root(self, event):
        if event.keysym not in ("Down", "Up", "Page_Down", "Page_Up"):
            return None
        w = self.focus_get()
        if w is not None and w.winfo_class() in (
                "Listbox", "Canvas", "FigureCanvasTkAgg"):
            self._page_by(1 if event.keysym in ("Down", "Page_Down") else -1)
        return None

    # ------------------------------------------------------------------ tour
    def _sync_current_btn(self):
        self.cur_btn_on = self.show_current.get()
        if getattr(self, "_cur_btn", None) is not None:
            if self.cur_btn_on:
                self._cur_btn.configure(bg="#DFF6DC", activebackground="#DFF6DC",
                                        fg="#14532D")
            else:
                self._cur_btn.configure(bg="#F2F4F7", activebackground="#F2F4F7",
                                        fg="#777777")

    def toggle_current(self):
        self.show_current.set(not self.show_current.get())
        self._sync_current_btn()
        self.status_var.set("Current trace (µA): " +
                            ("ON" if self.cur_btn_on else "OFF"))
        self.redraw()

    def toggle_tour(self):
        self.touring = not self.touring
        if self.touring:
            self.cycle_pos = 0
            if not self.selected:
                self.cycle_pos = 0
            self._tour_tick()
        else:
            if self._tour_job:
                self.after_cancel(self._tour_job)
                self._tour_job = None

    def _tour_tick(self):
        if not self.touring:
            return
        m = len(self.files)
        if m == 0:
            self.toggle_tour()
            return
        n = max(1, self.n_graphs.get())
        paths = [self.files[(self.cycle_pos + i) % m] for i in range(n)]
        self.selected = paths[:n]
        idxs = []
        byp = {str(p.resolve()): i for i, p in enumerate(self.files)}
        for p in paths:
            i = byp.get(str(p.resolve()))
            if i is not None:
                idxs.append(i)
        self.file_list.selection_clear(0, tk.END)
        for i in idxs:
            self.file_list.selection_set(i)
        if idxs:
            self.file_list.see(idxs[-1])
        self.cycle_pos = (self.cycle_pos + n) % m
        self.redraw()
        self._tour_job = self.after(max(300, self.cycle_ms.get()), self._tour_tick)

    # ------------------------------------------------------------------ analysis
    def _settings_from_ui(self) -> AnalysisSettings:
        s = AnalysisSettings(
            basecaller=self.basecaller.get(),
            base_order=self.base_order_var.get(),
            baseline_window=int(self.bl_win.get()),
            spectral_enable=self.spec_en.get(),
            position_adaptive_spectral=self.spec_adapt.get(),
            mobility_enable=self.mob_en.get(),
            use_gaussian_reconstruction=self.gauss_en.get(),
            use_multipass_wiener=self.mp_wiener.get(),
            gaussian_recon_segment_size=int(self.gauss_seg.get()),
            gaussian_recon_noise_reg=float(self.gauss_reg.get()),
            channel_peak_bonus=float(self.bonus.get()),
            pullback_weight=float(self.pullback.get()),
            ema_alpha=float(self.ema.get()),
            window_frac_lo=float(self.wlo.get()),
            window_frac_hi=float(self.whi.get()),
            view_mode=self.view_mode.get(),
            signal_start=int(self.scan_start.get()),
            signal_end=int(self.scan_end.get()),
        )
        return s

    def apply_knobs_redraw(self):
        self.settings = self._settings_from_ui()
        # re-load docs so the dye order / base_order takes effect on fresh loads
        self.redraw()

    def _ensure_doc(self, path: Path) -> TraceDocument:
        key = str(path.resolve())
        if key not in self.docs:
            self.status_var.set(f"Loading {path.name}…")
            self.update_idletasks()
            try:
                self.docs[key] = load_trace(path, base_order=self.base_order_var.get())
            except RuntimeError as e:
                messagebox.showerror("Load error", f"{path.name}:\n{e}")
                raise
            if len(self.docs) > DOC_CACHE_MAX:
                sel = {str(p.resolve()) for p in self.selected}
                for k in list(self.docs):
                    if k not in sel:
                        del self.docs[k]
                        if len(self.docs) <= DOC_CACHE_MAX:
                            break
        return self.docs[key]

    def basecall_selected(self):
        if not self.selected:
            messagebox.showinfo("Basecall", "Select one or more wells in the list.")
            return
        if self._busy:
            return
        self._busy = True
        self._cancel = False
        self.settings = self._settings_from_ui()
        total = len(self.selected)
        done = 0
        try:
            for path in self.selected:
                if self._cancel or not self.winfo_exists():
                    break
                try:
                    doc = self._ensure_doc(path)
                except Exception:
                    break
                try:
                    if doc.source == "rsd":
                        self.status_var.set(f"Basecalling {done + 1}/{total}: "
                                            f"{path.name} ({self.settings.basecaller})…")
                        run_basecall(doc, self.settings)
                    else:
                        self.status_var.set(f"{path.name}: no {self.basecaller.get()} "
                                            f"params for {doc.source} files — "
                                            f"showing container bases")
                except Exception as e:
                    messagebox.showerror("Basecall error", f"{path.name}:\n{e}")
                    break
                done += 1
                # keep the window (and the whole desktop) responsive on big jobs
                if total > 50 or done % 10 == 0:
                    if self.winfo_exists():
                        self.update()
            if self.winfo_exists():
                self.status_var.set(
                    f"Cancelled after {done} well(s)" if self._cancel
                    else f"Basecalled {total} well(s)")
                self.redraw()
                self._show_sequence()
        finally:
            self._busy = False

    def basecall_all(self):
        if not self.files:
            return
        n = len(self.files)
        if n > 96 and not messagebox.askyesno(
                "Basecall all",
                f"Basecall all {n} wells?\n\nThis runs for a long time and uses "
                "a lot of CPU/memory. You can press Esc to cancel.", parent=self):
            return
        self.selected = list(self.files)
        self.basecall_selected()

    def _show_sequence(self):
        self.seq_text.delete("1.0", tk.END)
        for path in self.selected[: self.n_graphs.get()]:
            key = str(path.resolve())
            doc = self.docs.get(key)
            if not doc:
                continue
            if not doc.sequence:
                self.seq_text.insert(
                    tk.END, f">{doc.well} ({doc.source}) — use ▶ Call to base-call\n\n")
                continue
            self.seq_text.insert(
                tk.END,
                f">{doc.well} {self._well_stats(doc)}  [{doc.source}]\n{doc.sequence}\n\n")

    # ------------------------------------------------------------------ plot
    def _theme_colors(self) -> dict:
        return dict(TRACE_THEMES.get(self.theme.get(), "Classic"))

    def _well_stats(self, doc) -> str:
        import textwrap
        seq, q = doc.sequence, doc.qualities
        n = len(seq)
        if n == 0:
            return "not basecalled"
        qa = np.asarray(q[:n], dtype=float) if q else np.array([], dtype=float)
        qmean = float(qa.mean()) if qa.size else 0.0
        qmin = float(qa.min()) if qa.size else 0.0
        nN = seq.count("N")
        gc = 100.0 * (seq.count("G") + seq.count("C")) / n if n else 0.0
        pos = np.asarray(doc.peak_positions[:n], dtype=float)
        sp = float(np.median(np.diff(pos))) if pos.size > 1 else 0.0
        return (f"len={n} Qmean={qmean:.1f} Qmin={qmin:.1f} "
                f"GC={gc:.1f}% N={nN} sp={sp:.2f} scans  [{doc.source}]")

    def _q100(self, q: np.ndarray) -> np.ndarray:
        """Map band-height qualities to a 0–100 score (manual-style profile)."""
        if q.size == 0:
            return q
        qq = np.asarray(q, dtype=float)
        scale = float(np.percentile(qq, 95))
        if scale <= 0:
            scale = float(qq.max()) or 1.0
        return np.clip(qq / scale * 100.0, 0.0, 100.0)

    def _draw_wrap(self, paths, settings, colors):
        rows = max(1, min(8, self.wrap_rows.get()))
        shows = [p for p in paths]
        total = len(shows) * rows
        cap = 20
        if total > cap:
            shows = shows[: max(1, cap // rows)]
            total = len(shows) * rows
        axes = self.fig.subplots(total, 1, sharex=True)
        if total == 1:
            axes = [axes]
        self._plot_axes = list(axes)
        gy0 = gy1 = None
        gx1 = 0
        for k, path in enumerate(shows):
            try:
                doc = self._ensure_doc(path)
            except Exception:
                ax = axes[k * rows]
                ax.text(0.5, 0.5, f"Load error: {path.name}", ha="center", transform=ax.transAxes)
                continue
            tr = display_trace(doc, settings)
            n = tr.shape[0]
            gx1 = max(gx1, n)
            segs = np.array_split(np.arange(n), rows)
            for r in range(rows):
                ax = axes[k * rows + r]
                sel = segs[r]
                for ci, base in enumerate(CHANNEL_ORDER):
                    if self.chan_show[ci].get():
                        y = tr[sel, ci]
                        ax.plot(sel, y, color=colors[base], lw=0.6, label=base if r == 0 else None)
                        if y.size:
                            m0 = float(np.nanmin(y)); m1 = float(np.nanmax(y))
                            gy0 = m0 if gy0 is None else min(gy0, m0)
                            gy1 = m1 if gy1 is None else max(gy1, m1)
                # current overlay (µA) in row 0 only
                if self.show_current.get() and r == 0 and doc.current_ua is not None:
                    cu = doc.current_ua
                    axc = ax.twinx()
                    axc.plot(np.arange(n), cu, color=CURRENT_COLOR, lw=0.8,
                             alpha=0.7, linestyle=":", label="I (µA)")
                    axc.set_ylabel("µA", color=CURRENT_COLOR, fontsize=6)
                    axc.tick_params(axis="y", labelcolor=CURRENT_COLOR, labelsize=6)
                ax.set_ylabel("", fontsize=7)
                if r == 0:
                    ax.text(0.004, 0.995, f"{doc.path.parent.name}/{doc.path.name}",
                            transform=ax.transAxes, ha="left", va="top",
                            fontsize=6, color="#333", zorder=6)
                ax.tick_params(labelsize=6)
                self._style_x_axis(ax, r == rows - 1,
                                   label=(r == rows - 1 and k == len(shows) - 1))
                self._draw_letters(ax, doc, sel[0], sel[-1] + 1, colors, force=True)
        for ax in axes:
            ax.grid(True, alpha=0.15)
        self.fig.suptitle(f"Wrap view — {len(shows)} well(s) × {rows} rows", fontsize=9)
        if gy0 is not None and gx1 > 0:
            pad = 0.02 * (gy1 - gy0) or 1.0
            self._full_xlim = (0.0, float(gx1))
            self._full_ylim = (float(gy0 - pad), float(gy1 + pad))
        return total

    def _draw_letters(self, ax, doc, s0, s1, colors, force=False):
        if not doc.sequence or not doc.peak_positions:
            return
        seq = doc.sequence
        pos = np.asarray(doc.peak_positions, dtype=float)
        q = doc.qualities or []
        med_sp = float(np.median(np.diff(pos))) if pos.size > 1 else 0.0
        tr = doc.acgt
        ytop = float(np.nanmax(tr[max(0, s0):min(tr.shape[0], s1)])) if s1 > s0 else 1.0
        tight = (s1 - s0) > 1400
        draw_let = force or (not tight) or med_sp >= 10
        for pi, base_ in enumerate(seq):
            if pi >= len(pos) or base_ not in colors:
                continue
            p = int(pos[pi])
            if not (s0 <= p < s1):
                continue
            ax.axvline(p, color=colors[base_], alpha=0.18, lw=0.4, zorder=1)
            if draw_let:
                ax.text(p, ytop + 0.015 * ytop, base_, color=colors[base_],
                        ha="center", va="bottom", fontsize=6.5, zorder=5, clip_on=True)
        return med_sp

    def _draw_quality(self, ax, doc, s0, s1):
        if not doc.sequence or not doc.qualities:
            return
        pos = np.asarray(doc.peak_positions, dtype=float)
        q = np.asarray(doc.qualities[: len(pos)], dtype=float)
        if q.size and q.max() > 0:
            q = self._q100(q)
        m = (pos >= s0) & (pos < s1)
        if not m.any():
            return
        x, y = pos[m], q[m]
        # poor (brown) everywhere, good (teal) overlaid on runs above threshold
        axb = ax.twinx()
        axb.plot(x, y, color=QUALITY_COLOR, lw=0.8, alpha=0.55, zorder=2)
        axb.fill_between(x, 0, y, color=QUALITY_COLOR, alpha=0.06, zorder=1)
        thresh = 32.0
        good = y >= thresh
        cross = np.zeros_like(good)
        cross[1:] = good[:-1]
        starts = np.where(good & ~cross)[0]
        ends = np.where(~good & (cross))[0]
        for st in starts:
            en = ends[ends > st]
            en = en[0] if len(en) else len(x)
            xt = x[st:en + 1]
            yt = y[st:en + 1]
            if len(xt) > 1:
                axb.plot(xt, yt, color=QUALITY_GOOD, lw=1.0, alpha=0.9, zorder=3)
        axb.set_ylim(0, 100)
        axb.set_ylabel("Q (0–100)", color="#555", fontsize=6)
        axb.tick_params(axis="y", labelsize=6, colors="#555")

    def _background_image(self):
        """Lazy-load every Background*/BG* picture shipped with the app (its own
        folder) or sitting next to it, and return one at random for the empty
        view.  Own-folder copies win, so a shared/zipped copy is self-contained."""
        if not self._bg_checked:
            self._bg_checked = True
            try:
                from PIL import Image, ImageOps
            except Exception:
                return None
            here = Path(__file__).resolve().parent
            seen_names = set()
            for root in (here, here.parent):
                if not root.is_dir():
                    continue
                for p in sorted(root.iterdir()):
                    if not p.is_file():
                        continue
                    stem = p.stem.lower()
                    if not (stem.startswith("background") or stem.startswith("bg")):
                        continue
                    if p.suffix.lower() not in (".jpg", ".jpeg", ".png", ".gif", ".bmp", ".webp"):
                        continue
                    key = p.name.lower()
                    if key in seen_names:
                        continue
                    seen_names.add(key)
                    try:
                        img = ImageOps.exif_transpose(Image.open(p)).convert("RGB")
                        img.thumbnail((1600, 1600))
                        self._bg_imgs.append(np.array(img))
                    except Exception:
                        continue
        if not self._bg_imgs:
            return None
        return random.choice(self._bg_imgs)

    def _on_right_click(self, _evt=None):
        """Right-click on the plot/axis bars resets the shared view."""
        self._reset_zoom()
        return "break"

    def _reset_axis(self, axis):
        (self._view_x if axis == "x" else self._view_y)[:] = [0.0, 1.0]
        self._apply_zoom()

    def _reset_zoom(self):
        self._view_x[:] = [0.0, 1.0]
        self._view_y[:] = [0.0, 1.0]
        if self._full_xlim is None:
            # Empty/background view: rebuild it so the picture is restored.
            self.redraw()
        else:
            self._apply_zoom()

    def _wheel_zoom(self, axis, direction):
        """Mouse-wheel over an axis bar: zoom that axis around its centre."""
        v = self._view_x if axis == "x" else self._view_y
        lo, hi = v
        width = max(1e-9, hi - lo)
        centre = 0.5 * (lo + hi)
        new_w = min(1.0, max(width / 1000.0, width * (0.8 if direction > 0 else 1.25)))
        lo2 = min(max(0.0, centre - new_w / 2.0), 1.0 - new_w)
        v[:] = [lo2, lo2 + new_w]
        self._apply_zoom()

    def _bar_cmd(self, axis, args):
        """ttk scrollbar -> pan that axis (grab the thumb and slide)."""
        if not args:
            return
        v = self._view_x if axis == "x" else self._view_y
        lo, hi = v
        width = hi - lo
        cmd = args[0]
        if cmd == "moveto" and len(args) > 1:
            new_lo = min(max(0.0, float(args[1])), max(0.0, 1.0 - width))
            v[:] = [new_lo, new_lo + width]
        elif cmd == "scroll" and len(args) > 1:
            n = float(args[1])
            step = width * (0.1 if len(args) > 2 and args[2] == "units" else 1.0)
            new_lo = min(max(0.0, lo + n * step), max(0.0, 1.0 - width))
            v[:] = [new_lo, new_lo + width]
        self._apply_zoom()

    def _apply_zoom(self, *_):
        """Push the shared X/Y view onto every visible graph + both bars."""
        if self._full_xlim is None or not self._plot_axes:
            # Empty/background view: make sure nothing (e.g. a toolbar zoom)
            # has pushed the picture out of frame.
            if self._full_xlim is None and self.fig.axes:
                for ax in self.fig.axes:
                    ax.set_xlim(0, 1)
                    ax.set_ylim(0, 1)
                    ax.set_autoscale_on(False)
                if hasattr(self, "xbar"):
                    self.xbar.set(0.0, 1.0)
                    self.ybar.set(0.0, 1.0)
                self.canvas.draw_idle()
            return
        x0, x1 = self._full_xlim
        y0, y1 = self._full_ylim
        fx0, fx1 = self._view_x
        fy0, fy1 = self._view_y
        w = (x1 - x0) * max(1e-9, fx1 - fx0)
        h = (y1 - y0) * max(1e-9, fy1 - fy0)
        xlo = x0 + (x1 - x0) * fx0
        ylo = y0 + (y1 - y0) * fy0
        for ax in self._plot_axes:
            ax.set_xlim(xlo, xlo + w)
            ax.set_ylim(ylo, ylo + h)
        if hasattr(self, "xbar"):
            self.xbar.set(fx0, fx1)
            self.ybar.set(fy0, fy1)
        self.canvas.draw_idle()

    def _fmt_time(self, scan, _pos=None):
        """Scan number -> time label (seconds, or m:ss past a minute)."""
        secs = float(scan) / SCAN_RATE_HZ
        if secs >= 60:
            m = int(secs // 60)
            return f"{m}:{secs - 60 * m:04.1f}"
        return f"{secs:.1f}"

    def _style_x_axis(self, ax, is_bottom, label=True):
        """Show x tick labels on the bottom pane only; label scans or time."""
        ax.tick_params(axis="x", labelbottom=bool(is_bottom))
        if not is_bottom:
            return
        if self.x_time.get():
            ax.xaxis.set_major_formatter(FuncFormatter(self._fmt_time))
        if label:
            ax.set_xlabel("Time (s)" if self.x_time.get() else "Scan")

    def redraw(self):
        self.fig.clear()
        self._plot_axes = []
        self._full_xlim = None
        self._full_ylim = None
        colors = self._theme_colors()
        n = max(1, min(8, self.n_graphs.get()))
        paths = self.selected[:n] if self.selected else []
        if not paths:
            ax = self.fig.add_subplot(111)
            bg = self._background_image()
            if bg is not None:
                ax.set_position([0, 0, 1, 1])
                ax.imshow(bg, extent=[0, 1, 0, 1], aspect="auto")
                ax.text(0.5, 0.96, "Select wells from the list",
                        ha="center", va="top", fontsize=15, color="#1A1A1A")
            else:
                ax.text(0.5, 0.5, "Select wells from the list\n(View → Sort files / Jump well to "
                                  "compare the same well across runs)",
                        ha="center", va="center")
            ax.set_axis_off()
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_autoscale_on(False)
            self.canvas.draw_idle()
            self.seq_text.delete("1.0", tk.END)
            return

        settings = self._settings_from_ui()
        wrap = settings.view_mode == "wrap"

        if wrap:
            self._draw_wrap(paths, settings, colors)
            self.fig.tight_layout()
            self._apply_zoom()
            self.canvas.draw_idle()
            self._show_sequence()
            return

        gx0 = gx1 = gy0 = gy1 = None
        for i, path in enumerate(paths):
            ax = self.fig.add_subplot(len(paths), 1, i + 1)
            try:
                doc = self._ensure_doc(path)
            except Exception as e:
                ax.text(0.5, 0.5, f"Load error: {e}", ha="center", transform=ax.transAxes)
                continue
            tr = display_trace(doc, settings)
            s0 = max(0, settings.signal_start)
            s1 = settings.signal_end if settings.signal_end > 0 else tr.shape[0]
            s1 = min(tr.shape[0], s1)
            x = np.arange(s0, s1)
            plotted = False
            for ci, base in enumerate(CHANNEL_ORDER):
                if not self.chan_show[ci].get():
                    continue
                seg = tr[s0:s1, ci]
                ax.plot(x, seg, color=colors[base], lw=0.7, label=base)
                if seg.size:
                    m0 = float(np.nanmin(seg)); m1 = float(np.nanmax(seg))
                    gy0 = m0 if gy0 is None else min(gy0, m0)
                    gy1 = m1 if gy1 is None else max(gy1, m1)
                plotted = True
            if s1 > s0:
                gx0 = s0 if gx0 is None else min(gx0, s0)
                gx1 = s1 if gx1 is None else max(gx1, s1)
            self._plot_axes.append(ax)
            if not plotted:
                ax.text(0.5, 0.5, "(all channels hidden)",
                        ha="center", va="center", transform=ax.transAxes, color="#888")

            if self.show_current.get() and doc.current_ua is not None:
                axc = ax.twinx()
                cu = doc.current_ua
                axc.plot(np.arange(s0, s1), cu[s0:s1], color=CURRENT_COLOR, lw=0.8,
                         alpha=0.8, linestyle=":", label="I (µA)")
                axc.set_ylabel("µA", color=CURRENT_COLOR, fontsize=6)
                axc.tick_params(axis="y", labelcolor=CURRENT_COLOR, labelsize=6)

            ax.set_ylabel("[Signal, V]" if i == 0 else "", fontsize=8)
            ax.tick_params(labelsize=7)
            ax.text(0.004, 0.995, f"{doc.path.parent.name}/{doc.path.name}",
                    transform=ax.transAxes, ha="left", va="top",
                    fontsize=6, color="#333", zorder=6)
            self._style_x_axis(ax, i == len(paths) - 1)

            if doc.sequence and settings.view_mode == "called" and doc.peak_positions:
                self._draw_letters(ax, doc, s0, s1, colors, force=False)
                ytop = float(np.nanmax(tr[s0:s1])) if s1 > s0 else 1.0
                if self.show_qnum.get():
                    q = doc.qualities or []
                    pos = np.asarray(doc.peak_positions, dtype=float)
                    for pi, base_ in enumerate(doc.sequence):
                        if pi >= len(pos):
                            break
                        p = int(pos[pi])
                        if not (s0 <= p < s1) or pi >= len(q):
                            continue
                        ax.text(p, -0.02 * ytop, f"{q[pi]:.1f}", color=QUALITY_COLOR,
                                ha="center", va="top", fontsize=6, zorder=5, clip_on=True)

            if (self.show_qcurve.get() and doc.sequence and doc.qualities
                    and settings.view_mode != "raw"):
                self._draw_quality(ax, doc, s0, s1)
        if gx0 is not None and gy0 is not None and gx1 > gx0:
            pad = 0.02 * (gy1 - gy0) or 1.0
            self._full_xlim = (float(gx0), float(gx1))
            self._full_ylim = (float(gy0 - pad), float(gy1 + pad))
        self.fig.tight_layout()
        self._apply_zoom()
        self.canvas.draw_idle()
        self._show_sequence()

    # ------------------------------------------------------------------ I/O
    def _docs_for_export(self):
        out = []
        for p in self.selected:
            doc = self.docs.get(str(p.resolve()))
            if doc:
                out.append(doc)
        return out

    def export_fasta(self):
        docs = self._docs_for_export()
        if not docs:
            messagebox.showinfo("Export", "Select loaded wells first.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".fasta", filetypes=[("FASTA", "*.fasta")])
        if not path:
            return
        lines = []
        for doc in docs:
            if doc.sequence:
                lines.append(f">{doc.well} len={len(doc.sequence)} [{doc.source}]\n{doc.sequence}")
        Path(path).write_text("\n".join(lines) + "\n")
        self.status_var.set(f"Wrote {path}")

    def export_peaks(self):
        """Export basecalled peaks to CSV: well, base, scan position, quality."""
        import csv
        docs = self._docs_for_export()
        if not docs:
            messagebox.showinfo("Export", "Select loaded wells first.")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if not path:
            return
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["well", "index", "base", "scan", "quality", "source"])
            for doc in docs:
                if not doc.sequence:
                    continue
                for i, (b, pos, q) in enumerate(zip(
                        doc.sequence, doc.peak_positions, doc.qualities or [])):
                    w.writerow([doc.well, i, b, int(pos), float(q), doc.source])
        self.status_var.set(f"Wrote {path}")

    def export_text(self):
        """Export trace text (manual 5.2): channels (V), current (µA), params."""
        import csv
        docs = self._docs_for_export()
        if not docs:
            messagebox.showinfo("Export", "Select loaded wells first.")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".txt", filetypes=[("Text", "*.txt")])
        if not path:
            return
        with open(path, "w", newline="") as fh:
            for doc in docs:
                fh.write(f"# WELL {doc.well}  SOURCE {doc.source}  {doc.meta}\n")
                fh.write("# Instrument parameters: source file "
                         f"{doc.path.name} ({doc.n_scans} scans)\n")
                if doc.sequence:
                    fh.write(f"# Bases called: {len(doc.sequence)}  start scan "
                             f"{doc.peak_positions[0] if doc.peak_positions else '-'} "
                             f"end scan "
                             f"{doc.peak_positions[-1] if doc.peak_positions else '-'}\n")
                w = csv.writer(fh, delimiter="\t", lineterminator="\n")
                hdr = ["scan"] + [f"{b}(V)" for b in CHANNEL_ORDER]
                if doc.current_ua is not None:
                    hdr.append("current(uA)")
                w.writerow(hdr)
                c = doc.current_ua
                for s in range(doc.n_scans):
                    row = [s] + [f"{doc.acgt[s, ci]:.4g}" for ci in range(4)]
                    if c is not None:
                        row.append(f"{c[s]:.4g}")
                    w.writerow(row)
                fh.write("\n")
        self.status_var.set(f"Wrote {path} ({len(docs)} well(s))")

    def save_figure(self):
        """Save the current plot area as an image (replaces the old toolbar)."""
        path = filedialog.asksaveasfilename(
            title="Save graph image",
            defaultextension=".png",
            filetypes=[("PNG image", "*.png"), ("PDF", "*.pdf"),
                       ("SVG", "*.svg"), ("All files", "*.*")])
        if not path:
            return
        try:
            self.fig.savefig(path, dpi=150, facecolor=self.fig.get_facecolor())
        except Exception as e:
            messagebox.showerror("Save graph image", str(e), parent=self)
            return
        self.status_var.set(f"Saved graph image: {path}")

    def save_settings(self):
        path = filedialog.asksaveasfilename(defaultextension=".json", filetypes=[("JSON", "*.json")])
        if not path:
            return
        s = self._settings_from_ui()
        Path(path).write_text(json.dumps(s.__dict__, indent=2))

    def load_settings(self):
        path = filedialog.askopenfilename(filetypes=[("JSON", "*.json")])
        if not path:
            return
        data = json.loads(Path(path).read_text())
        for k, v in data.items():
            if hasattr(self.settings, k):
                setattr(self.settings, k, v)
        self.basecaller.set(data.get("basecaller", "pos_bonus07"))
        self.base_order_var.set(data.get("base_order", "TGCA"))
        self.bonus.set(data.get("channel_peak_bonus", 0.7))
        self.pullback.set(data.get("pullback_weight", 0.008))
        self.status_var.set(f"Loaded settings from {path}")

    def show_help(self):
        """Scrollable in-app user manual (Help -> User manual, or F1)."""
        text = (
            "Limoncello CE Analyzer — user manual\n"
            "=====================================\n\n"
            "WHAT IT IS\n"
            "  A viewer and base caller for capillary-electrophoresis traces.\n"
            "  Supports .rsd, .scf, ABI .ab1 and text/CSV traces.\n"
            "  Fluorescence is shown in Volts; instrument current in µA.\n\n"
            "1. OPEN YOUR DATA\n"
            "  File ▸ Add data folder…  (or the  ⊕ Add…  button)\n"
            "    • Pick ONE run folder to load just that run.\n"
            "    • Ctrl-click SEVERAL folders in the picker to load them at once.\n"
            "    • Pick a parent folder that contains many run folders\n"
            "      (e.g. a whole OY run collection) to load every run under it.\n"
            "  The chosen runs appear in the  Data folders  list, top-left.\n\n"
            "2. THE DATA FOLDERS LIST\n"
            "  • Click a row to show only that run's files in Samples / files.\n"
            "  • Ctrl-click / Shift-click several rows to show their files together.\n"
            "  • Double-click a row to show that run and open its first sample.\n"
            "  • Select row(s) and click  – Remove  to unload those runs.\n"
            "  • File ▸ Clear folders unloads everything.\n\n"
            "3. THE SAMPLES / FILES LIST\n"
            "  Each entry is one trace (shown as  run/file).\n"
            "  • Click to plot it; Ctrl/Shift-click to plot several side by side.\n"
            "  • Set the  Graphs  spinner (bottom-left) to stack up to 8 plots.\n"
            "  • Type in  Jump/well  to filter by well name, file name or run name.\n"
            "  • Sort:  View ▸ Sort files  (well row, column, name, run, …).\n"
            "  • Remove selected samples with  – Remove file(s)  or the Delete key.\n\n"
            "4. MOVING AROUND THE PLOT\n"
            "  • Bottom bar = X (scan) axis, right bar = Y (signal) axis.\n"
            "  • X tick numbers appear only under the bottom pane, so stacked\n"
            "    plots keep their height. Tick  Time (s)  to read the X axis in\n"
            "    seconds instead of scan numbers (1.75 scans/s, ~0.57 s/scan).\n"
            "  • Grab a bar's thumb and slide to pan that axis.\n"
            "  • Roll the mouse wheel over the X bar to zoom X; over the Y bar to\n"
            "    zoom Y. Double-click a bar to reset that axis.\n"
            "  • Both bars move every visible plot together, so comparisons stay\n"
            "    aligned.\n"
            "  • Lost in the zoom?  Right-click the plot or a bar, press Home, or\n"
            "    use the  ⟲ Reset view  button (View ▸ Reset view) to start over.\n\n"
            "5. VIEW MODES  (View menu)\n"
            "  • Raw traces            — detector signal as recorded.\n"
            "  • Processed (ACGT)      — colour-separated channels.\n"
            "  • Base-called           — peaks with base letters and quality.\n"
            "  • Wrap                  — one well split over N rows; set\n"
            "                            View ▸ Wrap rows (2-8).\n"
            "  • Channels submenu toggles A/C/G/T; Trace colors picks a palette.\n"
            "  • Page forward / Page backward step through wells; auto-tour plays\n"
            "    them automatically (Space starts/stops it).\n\n"
            "6. BASECALLING  (Analysis menu)\n"
            "  • Basecall selected — calls the wells currently plotted.\n"
            "  • Basecall all in list — calls everything listed (asks first if big).\n"
            "  • ⚙ Basecall settings… — advanced caller parameters.\n"
            "  • Called sequence(s) appear in the box below the plot.\n"
            "  • Genotyping is planned for a later release.\n\n"
            "7. EXPORT  (File menu)\n"
            "  • Export sequence (FASTA)…   called bases per well.\n"
            "  • Export peak table (CSV)…   well, base, scan position, quality.\n"
            "  • Export trace text (V + µA)… raw values.\n"
            "  • Save graph image…          the plot area as PNG/PDF/SVG.\n"
            "  • Save / Load settings JSON…  remembers your caller setup.\n\n"
            "8. KEYBOARD SHORTCUTS\n"
            "  ↑ / ↓ / PgUp / PgDn   page through wells\n"
            "  Space                 start / stop auto-tour\n"
            "  Delete                remove selected sample(s)\n"
            "  Home / right-click    reset the X & Y view\n"
            "  F1                    this manual\n"
            "  F11                   zoom the window to full screen / restore\n"
            "  Esc                   cancel a long job, or leave full screen\n\n"
            "TIPS\n"
            "  • Comparing the same well across runs? Sort by well and use\n"
            "    Jump/well, then select the matching samples.\n"
            "  • Large collections take a moment to list; plotting is lazy.\n"
            "  • The empty graph area shows the project picture until you pick\n"
            "    a sample.\n"
        )
        win = tk.Toplevel(self)
        win.title("Limoncello CE Analyzer — User manual")
        win.geometry("780x660")
        win.minsize(520, 400)
        box = scrolledtext.ScrolledText(win, wrap=tk.WORD,
                                        font=("DejaVu Sans Mono", 10),
                                        padx=14, pady=10, bg="#FBFBF4")
        box.pack(fill=tk.BOTH, expand=True)
        box.insert("1.0", text)
        box.configure(state=tk.DISABLED)
        ttk.Button(win, text="Close", command=win.destroy).pack(pady=6)
        win.transient(self)

    def show_about(self):
        messagebox.showinfo(
            "About",
            "Limoncello CE Analyzer  —  CE trace viewer & base caller\n\n"
            "Serves all instruments that produce CE traces:\n"
            "Formats: .rsd, .scf, ABI .ab1, text/CSV traces\n"
            "Current trace: µA (RSD raw ÷10); fluorescence in Volts.\n"
            "Basecallers: pos_bonus07, pos_profile, hz_soften, raw_peaks\n"
            "Multi-folder load, multi-graph, auto-tour, well sort/filter.\n"
            "Genotyping is planned for a later release.\n\n"
            "Uses our own tuned spacing tracker.",
        )


def main():
    ap = argparse.ArgumentParser(description="Limoncello CE Analyzer")
    ap.add_argument("--folder", action="append", type=Path, help="Data folder (repeatable)")
    args = ap.parse_args()
    app = LimoncelloAnalyzerApp(initial_folders=args.folder)
    app.mainloop()


if __name__ == "__main__":
    main()