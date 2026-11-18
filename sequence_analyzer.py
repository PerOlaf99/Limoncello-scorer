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
  • Basecaller versions (mb1000/mb4000 × accuracy/length, legacy aliases)
  • Auto-tour (slideshow) through the loaded wells
  • Sort wells by well / name / folder, jump-to-well filter
  • Export: FASTA, peak CSV, trace text (channels V + current µA)

Base calling uses our own tuned caller (best_basecaller spacing tracker),
not commercial instrument algorithms.

Genotyping sits next to base calling in the Genotyping menu: manual peak
picking straight on the stacked viewer (batch 4-6 wells at a time, click to
pick, CTC-CE duplex internal standard, CSV/Excel/JSON export for ML).

Run:
  python sequence_analyzer.py
  python sequence_analyzer.py --folder /path/to/rsd_dir
"""
from __future__ import annotations

import sys  # noqa: I001 -- needed before the GUI imports below

# ``--check`` must work even where tkinter/matplotlib are not installed, so
# answer it before the GUI imports below are executed.
if __name__ == "__main__" and "--check" in sys.argv[1:]:
    from analyzer_core import environment_report
    sys.exit(environment_report())

import argparse
from types import SimpleNamespace
import json
import random
import re
import shutil
import subprocess
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext, simpledialog
from pathlib import Path
from typing import List, Optional

import numpy as np

import mark_library

# matplotlib embedded
import matplotlib

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter

from analyzer_core import (
    AnalysisSettings,
    BASECALLER_VERSIONS,
    TraceDocument,
    discover_files,
    find_run_folders,
    declared_dye_order,
    load_trace,
    run_basecall,
    display_trace,
    run_params_summary,
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
    # Classic = the MegaBACE software look, colours FIXED per base letter:
    # A green, C blue, T red, G black.  Traces are drawn Ch1..Ch4 in the run's
    # base order (default A·C·T·G), so on a default run Ch1=A (green),
    # Ch2=C (blue), Ch3=T (red), Ch4=G (black) — exactly as on the MegaBACE
    # sequence analyser.  The base order only shifts which channel shows
    # which colour; colours stay with the letters.
    "Classic": {"A": "#00AA00", "C": "#0000DD", "G": "#111111", "T": "#DD0000"},
    # Sequencing (MegaBACE SA) colours are FIXED per base, exactly as the SA
    # software shows them: A green, C blue, T red, G black.  The base order
    # changes per chemistry (it only decides which base sits on which physical
    # channel), so the traces keep these colours at all times.  The dye-kit
    # names below therefore match the equivalent of these colours on MegaBACE.
    "Seq DYEnamic (T·G·C·A)": {"A": "#00AA00", "C": "#0000DD", "G": "#111111", "T": "#DD0000"},
    "Seq ET primer (A·C·T·G)": {"A": "#00AA00", "C": "#0000DD", "G": "#111111", "T": "#DD0000"},
    "Chromas": {"A": "#00AA00", "C": "#1E90FF", "G": "#444444", "T": "#FF0000"},
    "High-contrast": {"A": "#2E8B57", "C": "#1F4FC0", "G": "#000000", "T": "#E03030"},
    "Monochrome": {"A": "#777777", "C": "#555555", "G": "#333333", "T": "#888888"},
    # MegaBACE genotyping dye sets are FIXED per CHANNEL (no base colours),
    # applied physically: Channel1,2,3,4.  Dye set 1 = Ch1 Red (ET-ROX),
    # Ch2 Blue (FAM), Ch3 Black (HEX/NED), Ch4 Green (TET/HEX).  Many users
    # are used to Red=Ch3 / Black=Ch4, hence the second variant.
    "Genotyping (R·B·Blk·G)": {"A": "#E03030", "C": "#2557E8", "G": "#111111", "T": "#00B052"},
    "Genotyping (G·B·R·Blk)": {"A": "#00B052", "C": "#2557E8", "G": "#E03030", "T": "#111111"},
}

# "channel" themes paint each trace by position (Channel1..4), ignoring which
# base letter the channel carries (genotyping dye sets).  Everything else,
# including Classic, paints by the colour of the base letter (MegaBACE SA).
TRACE_THEME_MODE = {
    "Genotyping (R·B·Blk·G)": "channel",
    "Genotyping (G·B·R·Blk)": "channel",
}

SRC_FG = {
    "rsd": "#000000",
    "scf": "#1E55C0",
    "abi": "#0F7A4D",
    "text": "#B45A00",
}

# well sorting modes: "letter first" = A01,A02,…; "number first" = A01,B01,…H01,A02…
def _event_mods(event):
    """(ctrl, shift) held during a matplotlib mouse event."""
    key = (getattr(event, "key", None) or "").lower()
    ge = getattr(event, "guiEvent", None)
    state = getattr(ge, "state", 0) if ge is not None else 0
    try:
        state = int(state)
    except (TypeError, ValueError):
        state = 0
    ctrl = "control" in key or "ctrl" in key or bool(state & 0x0004)
    shift = "shift" in key or bool(state & 0x0001)
    return ctrl, shift


# Most graphs stacked at once.  Every place that reads the Graphs box goes
# through _graph_count(), so the box, the viewer and the picking modes agree.
MAX_GRAPHS = 12

# Line styles that tell runs apart when "Overlay runs" puts the same well from
# several runs into one graph (colour stays the dye's).
RUN_STYLES = ["-", "--", ":", "-."]

SORT_OPTIONS = [
    ("well_row", "Well, row order"),
    ("well_col", "Well, column order"),
    ("name", "File name"),
    ("folder", "Run, then file"),
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
        # The system title bar stays: minimise, maximise, edge-resize and window
        # snapping come from the window manager on every platform.  (A
        # 'splash' window with a drawn title bar lost them on Linux and was
        # ignored on Windows/macOS, which then showed two title bars.)
        self.geometry("1440x920")
        self.protocol("WM_DELETE_WINDOW", self._on_close)
        self._zoomed = False
        self._pre_zoom = ""
        self._busy = False
        self._cancel = False
        self.bind("<Escape>", self._on_escape)
        self.bind("<F11>", lambda e: self._toggle_zoom())
        self.bind("<F1>", lambda e: self.show_help())
        self.bind("<Home>", lambda e: self._reset_zoom())
        # align the graph under the mouse: Alt+←/→ one scan, Alt+Shift ten
        self.bind("<Alt-Left>", lambda e: self._nudge(-1))
        self.bind("<Alt-Right>", lambda e: self._nudge(1))
        self.bind("<Alt-Shift-Left>", lambda e: self._nudge(-10))
        self.bind("<Alt-Shift-Right>", lambda e: self._nudge(10))
        self.bind("<Alt-Home>", lambda e: self._reset_alignment("hover"))
        self.bind("<Control-z>", lambda e: self._gen_undo())
        self.bind("<Control-Z>", lambda e: self._gen_undo())
        # Hold R and click a picked peak to remove it (its +A tag too).  Kept
        # at the window level and checked in the plot click, because matplotlib
        # does not carry a plain letter through a mouse event the way it does
        # Shift/Ctrl/Alt.
        self._remove_mod = False
        self.bind("<KeyPress-r>", self._on_remove_key)
        self.bind("<KeyPress-R>", self._on_remove_key)
        self.bind("<KeyRelease-r>", self._on_remove_release)
        self.bind("<KeyRelease-R>", self._on_remove_release)

        self.folders: List[Path] = list(initial_folders or [])
        self.files: List[Path] = []
        self.docs: dict[str, TraceDocument] = {}
        self._comment_cache: dict[str, str] = {}
        self.selected: List[Path] = []
        self.excluded: set[str] = set()   # individual files hidden by the user
        self._bg_checked = False
        self._bg_imgs: List[np.ndarray] = []
        self.settings = AnalysisSettings()
        self.base_order_var = tk.StringVar(value=self.settings.base_order)
        self.n_graphs = tk.IntVar(value=1)
        # View options (View menu + the view bar), remembered across sessions.
        _vp = mark_library.load_prefs().get("view", {})
        self.y_scale = tk.StringVar(
            value=_vp.get("y_scale", "shared") if _vp.get("y_scale") in
            ("shared", "each") else "shared")
        self.fit_y = tk.BooleanVar(value=bool(_vp.get("fit_y", True)))
        self.always_table = tk.BooleanVar(value=bool(_vp.get("always_table", False)))
        self.overlay_runs = tk.BooleanVar(value=bool(_vp.get("overlay_runs", False)))
        self.view_mode = tk.StringVar(value="processed")
        self.basecaller = tk.StringVar(value="mb1000_accuracy")
        self.status_var = tk.StringVar(value="Ready — add a data folder to begin")

        # Visualization state
        self.chan_show = [tk.BooleanVar(value=True) for _ in CHANNEL_ORDER]
        self.show_letters = tk.BooleanVar(value=True)
        self.show_qnum = tk.BooleanVar(value=False)
        self.show_qcurve = tk.BooleanVar(value=False)
        self.show_current = tk.BooleanVar(value=True)
        self.x_time = tk.BooleanVar(value=False)  # x axis: scans (off) vs time (on)
        self.theme = tk.StringVar(value="Classic")
        self.wrap_rows = tk.IntVar(value=5)

        # Main-window manual genotyping: the user pages the current batch
        # (n_graphs wells at a time), zooms in and clicks peaks on the stacked
        # viewer.  One PeakPicker per displayed well accumulates the records.
        self.genotyping_active = False
        # area_mode reuses this whole mode for the plainer drag-to-measure
        # tool: same stacked viewer, same table, same batching, but a
        # press-drag-release span instead of a click that detects a peak.
        self.area_mode = False
        self._gen_pickers: dict = {}
        self._area_pickers: dict = {}
        self._area_drag = None
        self._area_cursors: list = []
        self._gen_active_path: Optional[Path] = None
        # Auto-genotyping: which physical channel carries the internal
        # standard and which carries the sample.  These are per-run assay
        # roles, not dye names -- they start at the T9 kit's Ch3/Ch2 but the
        # user can move them, and the pairing is checked before a single well
        # is scored, because swapping them does not fail loudly: it scores the
        # sample's own peaks as the standard and returns confident nonsense.
        self._auto_is_channel = tk.IntVar(value=3)
        self._auto_sample_channel = tk.IntVar(value=2)
        # The operator's floor on the minor allele before a well may be called
        # a het: 0 keeps the engine's own behaviour, where an off-balance het
        # is called and only flagged.  It is remembered with the marks, because
        # a setting that shaped a plate should not quietly vanish on restart.
        self._auto_min_het_frac = tk.DoubleVar(
            value=float(mark_library.load_prefs().get("min_het_frac", 0.0)))
        self._auto_rows: list = []
        # The plate model learned by the semi-automatic pass, if one has been
        # run.  Auto-genotyping uses it when present and refuses to guess a
        # standard when it is not: there is no built-in geometry to fall back
        # on, because one plate's internal standard does not predict another's.
        self._auto_is_model = None
        self._auto_is_model_note = ""
        # One standard shape per run: a window can hold several runs at once
        # and each has its own drift, so the model is looked up by run rather
        # than assumed to be the only one there is.
        self._auto_models: dict = {}
        # A standard fitted on one plate and loaded here to seed another, used
        # only for a run that has no marks of its own to fit from.  Cleared by
        # "Forget learned standard" like everything else.
        self._auto_saved_model = None
        self._auto_saved_run = ""
        # Saved marks on disk (library/), read per run on first use.
        self._library_cache: dict = {}
        self._auto_tree = None
        self._auto_table = None
        # What sits under the graphs: the genotyping table by default, the
        # basecalled sequence when the toggle asks for it.  Manual picking and
        # area measuring take the slot over while they are on.
        self._below_mode = "geno"
        self._geno_rows: list = []          # the table's rows, in display order
        self._auto_iids: dict = {}          # (run, well) -> tree item id
        self._auto_cols: tuple = ()         # the table's columns, in order
        # Popup for editing a call cell off the auto-genotype table; the window
        # is re-raised, not recreated, while one edit is in progress.
        self._auto_edit_win = None
        self._pick_cid = None
        self._gen_motion_cid = None
        self._area_release_cid = None
        self._gen_cursors: list = []
        self.pick_table = None
        self.pick_tree = None
        self._gen_sh = tk.BooleanVar(value=False)
        self._gen_d2 = tk.BooleanVar(value=False)
        # The two pick modes are radio buttons in effect: at most one can own
        # the canvas, and the menu shows which one is live.
        self._mode_pick = tk.BooleanVar(value=False)
        self._mode_area = tk.BooleanVar(value=False)
        self._gen_lock_var = tk.StringVar(value="any")

        # Shared X/Y view, stored as [first, last] fractions of the full data
        # range. One model drives every visible graph plus the axis bars.
        self._view_x = [0.0, 1.0]
        self._view_y = [0.0, 1.0]
        self._plot_axes: List = []
        self._full_xlim = None
        self._full_ylim = None
        self._wrap_active = False      # wrap rows: each row its own X slice
        self._plot_drag = None         # box-zoom / pan gesture in progress
        # Per-graph alignment: file -> (a, b) with  reference = a * scan + b.
        # Display only: every graph keeps its true scans (picks, tables and
        # marks never move); alignment changes WHICH scans a graph shows.
        self._align: dict = {}
        self._axes_keys: List = []     # file key per graph (None: not alignable)
        self._axes_titles: List = []   # (text artist, base label) per graph
        self._data_full = None         # union of the graphs' true scan ranges
        self._hover_idx = None         # graph under the mouse (for nudging)

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

    def _toggle_zoom(self):
        """F11: full screen and back, through the window manager."""
        try:
            self._zoomed = not bool(self.attributes("-fullscreen"))
            self.attributes("-fullscreen", self._zoomed)
            return
        except tk.TclError:
            pass
        if self._zoomed:                       # no fullscreen support: fill the screen
            self.geometry(self._pre_zoom)
            self._zoomed = False
        else:
            self._pre_zoom = self.geometry()
            self.geometry(f"{self.winfo_screenwidth()}x{self.winfo_screenheight()}+0+0")
            self._zoomed = True

    def _on_close(self):
        """Write anything still waiting (batched alignment saves) and quit."""
        try:
            if getattr(self, "_align_save_job", None) is not None:
                self.after_cancel(self._align_save_job)
                self._flush_alignment()
        except Exception:                                  # noqa: BLE001
            pass
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
        elif getattr(self, "genotyping_active", False):
            self.exit_genotyping_picking()

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
        file_m.add_command(label="Exit", command=self._on_close)

        view_m = tk.Menu(self, tearoff=0)
        for mode, label in [
            ("raw", "Raw traces"),
            ("processed", "Processed (ACGT)"),
            ("called", "Base-called (with peaks)"),
            ("esd", "Sequencing trace (ESD peaks)"),
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
        yscale_m = tk.Menu(self, tearoff=0)
        for code, label in (("shared", "Shared by all graphs"),
                            ("each", "Each graph its own")):
            yscale_m.add_radiobutton(
                label=label, variable=self.y_scale, value=code,
                command=lambda: (self._view_y.__setitem__(slice(None), [0.0, 1.0]),
                                 self._sync_y_scale_cb(), self._view_option_changed()))
        view_m.add_cascade(label="Y scale", menu=yscale_m)
        view_m.add_checkbutton(label="Fit Y to view", variable=self.fit_y,
                               command=self._view_option_changed)
        view_m.add_checkbutton(label="Overlay runs of the same well",
                               variable=self.overlay_runs,
                               command=lambda: self._view_option_changed(redraw=True))
        view_m.add_command(label="Plate map…", command=self.show_plate_map)
        view_m.add_cascade(label="Alignment", menu=self._alignment_menu(view_m))
        view_m.add_checkbutton(label="Always show table under graphs",
                               variable=self.always_table,
                               command=lambda: (self._save_view_prefs(),
                                                self._below_refresh()))
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

        basecall_m = tk.Menu(self, tearoff=0)
        basecall_m.add_command(label="Basecall selected", command=self.basecall_selected)
        basecall_m.add_command(label="Basecall all in list", command=self.basecall_all)
        basecall_m.add_separator()
        basecall_m.add_command(label="⚙ Basecall settings…", command=self.open_settings)
        basecall_m.add_separator()
        basecall_m.add_radiobutton(
            label="Sequencing trace (ESD peaks)",
            variable=self.view_mode, value="esd", command=self.redraw)
        basecall_m.add_radiobutton(
            label="Plain channels (processed)",
            variable=self.view_mode, value="processed", command=self.redraw)
        basecall_m.add_separator()
        basecall_m.add_command(label="Clear base calls (undo)", command=self.clear_basecalls)

        # The picking workflow and the automated one used to share one menu,
        # with Auto-genotyping buried under a cascade; an operator doing the
        # automated pass never found it.  They are two jobs, so they get two
        # top-level menus.
        genotyping_m = tk.Menu(self, tearoff=0)
        genotyping_m.add_command(label="Manual peak picking…",
                                 command=self.toggle_genotyping_picking)
        genotyping_m.add_command(label="Peak area by drag…",
                                 command=self.toggle_area_picking)
        genotyping_m.add_separator()
        genotyping_m.add_command(label="Mark peaks as standard…",
                                 command=self._gen_mark_std)
        genotyping_m.add_command(label="Mark IS in all shown wells",
                                 command=self._gen_mark_std_batch)
        genotyping_m.add_command(label="Marked IS wells…",
                                 command=self.report_marked_wells)
        genotyping_m.add_separator()
        genotyping_m.add_command(label="Undo last pick / measurement",
                                 command=self._gen_undo)
        genotyping_m.add_command(label="Clear picks / measurements",
                                 command=self._gen_clear)
        genotyping_m.add_command(label="Tag duplex species for MF…",
                                 command=self._gen_mark_duplex)
        genotyping_m.add_command(label="Save peaks table…", command=self._gen_save)
        genotyping_m.add_command(label="Fragment-length sizing (ladder)…",
                                 command=self.size_selected_fragments)
        genotyping_m.add_separator()
        genotyping_m.add_cascade(label="Picking options…", menu=self._build_picking_opts(genotyping_m))
        genotyping_m.add_separator()
        genotyping_m.add_command(label="Next batch →",
                                 command=lambda: self._page_by(1))
        genotyping_m.add_command(label="← Previous batch",
                                 command=lambda: self._page_by(-1))
        genotyping_m.add_separator()
        genotyping_m.add_command(label="Exit peak picking",
                                 command=self.exit_genotyping_picking)

        auto_m = tk.Menu(self, tearoff=0)
        auto_m.add_command(label="Semi automatic standard…",
                           command=self.semi_auto_is)
        auto_m.add_command(label="Use a saved standard…",
                           command=self.auto_use_saved_model)
        auto_m.add_command(label="Forget learned standard",
                           command=self.forget_plate_model)
        auto_m.add_separator()
        auto_m.add_command(label="Auto-genotype sample window…",
                           command=self.auto_genotype_wells)
        auto_m.add_command(label="Channel roles (standard / sample)…",
                           command=self.auto_genotype_channels)
        auto_m.add_command(label="Minimum minor allele for a het…",
                           command=self.auto_genotype_min_frac)
        auto_m.add_separator()
        auto_m.add_command(label="Save auto-genotype table…",
                           command=self.auto_genotype_save)
        auto_m.add_command(label="Export library for ML…",
                           command=self.library_export_ml)
        auto_m.add_separator()
        auto_m.add_command(label="CTCE analysis…", command=self.show_ctce)

        help_m = tk.Menu(self, tearoff=0)
        help_m.add_command(label="User manual…", command=self.show_help)
        help_m.add_command(label="About", command=self.show_about)

        comments_m = tk.Menu(self, tearoff=0)
        comments_m.add_command(label="Dye order of this run…",
                               command=self.set_run_dye_order)
        comments_m.add_command(label="Run comments…", command=self.edit_comments)
        comments_m.add_command(label="Run info…", command=self.show_run_info)

        self._menus = {"File": file_m, "View": view_m,
                       "Base calling": basecall_m,
                       "Manual genotyping": genotyping_m,
                       "Auto genotyping": auto_m,
                       "Comments": comments_m, "Help": help_m}

    def _build_picking_opts(self, parent):
        m = tk.Menu(parent, tearoff=0)
        m.add_checkbutton(
            label="Add +A (A-addition) peak",
            variable=self._gen_sh, command=self._sync_gen_opts)
        m.add_checkbutton(
            label="Mark start/end from the 2nd derivative",
            variable=self._gen_d2, command=self._sync_gen_opts)
        m.add_separator()
        m.add_radiobutton(
            label="Pick on any channel (auto)",
            variable=self._gen_lock_var, value="any",
            command=self._sync_gen_opts)
        for _ci in (1, 2, 3, 4):
            m.add_radiobutton(
                label=f"Lock picks to Ch{_ci}",
                variable=self._gen_lock_var, value=str(_ci),
                command=self._sync_gen_opts)
        return m

    def _build_layout(self):
        # Yellow banner with the app's name; the window buttons are the
        # system's own.
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
        ttk.Button(fb, text="– Remove run", command=self.remove_folder).pack(side=tk.LEFT, padx=4)

        self.samples_label = ttk.Label(left, text="◧ Samples / files",
                                       font=("", 10, "bold"),
                                       foreground="#0F3A6E")
        self.samples_label.pack(anchor=tk.W, padx=4, pady=(8, 0))

        # sort + filter row
        srow = ttk.Frame(left)
        srow.pack(fill=tk.X, padx=4)
        ttk.Label(srow, text="Sort:").pack(side=tk.LEFT)
        self.sort_label_var = tk.StringVar(
            value=SORT_LABELS.get(self.sort_mode.get(), SORT_LABELS["well_row"]))
        self.sort_mode.trace_add(
            "write", lambda *_a: self.sort_label_var.set(
                SORT_LABELS.get(self.sort_mode.get(), self.sort_label_var.get())))
        self.sort_cb = ttk.Combobox(srow, textvariable=self.sort_label_var, width=18,
                                    values=list(SORT_LABELS.values()), state="readonly")
        self.sort_cb.pack(side=tk.LEFT, padx=(2, 4))
        self.sort_cb.bind("<<ComboboxSelected>>",
                          lambda e: (self.sort_mode.set(SORT_CODES[self.sort_label_var.get()]),
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
        self.file_list.bind("<Left>", self._page_key)
        self.file_list.bind("<Right>", self._page_key)
        self.file_list.bind("<Page_Down>", self._page_key)
        self.file_list.bind("<Page_Up>", self._page_key)
        # Help and README promise Delete removes the selected samples.
        self.file_list.bind("<Delete>", lambda e: self.remove_samples())

        self.bind("<Down>", self._page_key_root)
        self.bind("<Up>", self._page_key_root)
        self.bind("<Left>", self._page_key_root)
        self.bind("<Right>", self._page_key_root)
        self.bind("<Page_Down>", self._page_key_root)
        self.bind("<Page_Up>", self._page_key_root)

        rrow = ttk.Frame(left)
        rrow.pack(fill=tk.X, padx=4, pady=(0, 2))
        ttk.Button(rrow, text="– Hide file(s)",
                   command=self.remove_samples).pack(side=tk.LEFT)

        frow = ttk.Frame(left)
        frow.pack(fill=tk.X, padx=4, pady=(0, 2))
        ttk.Label(frow, text="Tour (ms)").pack(side=tk.LEFT)
        tk.Spinbox(
            frow, from_=500, to=20000, increment=500, textvariable=self.cycle_ms,
            width=5, highlightthickness=0
        ).pack(side=tk.LEFT, padx=2)
        ttk.Button(frow, text="▶", width=3, command=self.toggle_tour).pack(side=tk.LEFT, padx=2)

        # ----- Center: plots -----
        center = ttk.Frame(body)
        self.center = center
        body.add(center, weight=5)

        chan_bar = ttk.Frame(center)
        chan_bar.pack(fill=tk.X, padx=4, pady=(4, 0))
        ttk.Label(chan_bar, text="Channels:").pack(side=tk.LEFT)
        self._chan_cbs = []
        for i in range(len(CHANNEL_ORDER)):
            # widget i drives chan_show[i], and chan_show is indexed by
            # physical channel -- so widget i IS Ch i+1 and the widgets are
            # already in order.  Only the label needs the dye order.
            cb = tk.Checkbutton(
                chan_bar, text=f"  Ch{i + 1} {CHANNEL_ORDER[i]}  ", variable=self.chan_show[i],
                command=self._chan_changed, bg="#F2F4F7",
                fg="#000000", activebackground="#F2F4F7",
                selectcolor="white")
            cb.pack(side=tk.LEFT, padx=2)
            self._chan_cbs.append(cb)
        # the bar starts out indexed A,C,G,T; relabel to the run's dye order
        self._refresh_channel_labels()
        self.base_order_var.trace_add("write", self._on_base_order_changed)
        self._refresh_channel_buttons()
        self._cur_btn = tk.Checkbutton(
            chan_bar, text="  µA  ", variable=self.show_current,
            bg="#DFF6DC", activebackground="#DFF6DC", fg="#14532D",
            selectcolor="white",
            command=lambda: (self._sync_current_btn(), self.redraw()))
        self._cur_btn.pack(side=tk.LEFT, padx=(8, 2))
        self._sync_current_btn()
        self._time_btn = tk.Checkbutton(
            chan_bar, text="  Time (min)  ", variable=self.x_time,
            bg="#F2F4F7", activebackground="#F2F4F7", selectcolor="white",
            command=self.redraw)
        self._time_btn.pack(side=tk.LEFT, padx=(8, 2))
        self.theme_cb = ttk.Combobox(
            chan_bar, textvariable=self.theme, width=24, state="readonly",
            values=list(TRACE_THEMES))
        self.theme_cb.pack(side=tk.RIGHT, padx=(4, 8))
        self.theme_cb.bind("<<ComboboxSelected>>", lambda e: self.redraw())
        ttk.Label(chan_bar, text="Trace colors:").pack(side=tk.RIGHT, pady=4)
        self._chan_bar = chan_bar

        self._build_view_bar(center)
        self._build_mode_banner(center)

        # Graphs above, the table for the current mode below, with a sash the
        # user can drag; the lower pane is taken away when it has nothing to
        # show (see _below_refresh).
        self.vsplit = ttk.Panedwindow(center, orient=tk.VERTICAL)
        self.vsplit.pack(fill=tk.BOTH, expand=True)
        plot_pane = ttk.Frame(self.vsplit)
        self.below = ttk.Frame(self.vsplit)
        self.vsplit.add(plot_pane, weight=4)
        self._below_shown = False
        self._below_sash = None

        # Axis bars: drag to pan, mouse-wheel to zoom. One bar per axis and it
        # drives every visible graph at once (shared X / shared Y).
        plotf = self.plotf = ttk.Frame(plot_pane)
        plotf.pack(fill=tk.BOTH, expand=True)
        plotf.rowconfigure(0, weight=1)
        plotf.columnconfigure(0, weight=1)

        self.fig = Figure(figsize=(9, 7), dpi=100, facecolor="#FFFFFF")
        self.canvas = FigureCanvasTkAgg(self.fig, master=plotf)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        self.canvas.get_tk_widget().bind("<Button-3>", self._on_right_click)
        # On the graphs themselves: wheel zooms at the cursor (Ctrl: Y,
        # Shift: pan), drag zooms a box (Ctrl+drag while picking), middle-drag
        # pans, double-click resets.
        self.canvas.mpl_connect("resize_event", self._on_canvas_resize)
        self.canvas.mpl_connect("scroll_event", self._on_plot_scroll)
        self.canvas.mpl_connect("button_press_event", self._on_plot_press)
        self.canvas.mpl_connect("motion_notify_event", self._on_plot_motion)
        self.canvas.mpl_connect("button_release_event", self._on_plot_release)
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

        # Below the graphs: the genotyping table for the wells being shown, so
        # a call and its standard can be checked while the traces are on
        # screen.  The basecalled sequence is one toggle away, and manual peak
        # picking / area measuring swap the same slot for their own table.
        hdr = ttk.Frame(self.below)
        hdr.pack(fill=tk.X, padx=4)
        self.seq_hdr = ttk.Label(hdr, text="Genotyping")
        self.seq_hdr.pack(side=tk.LEFT)
        self.seq_btn = ttk.Button(hdr, text="Show sequence",
                                  command=self._toggle_below, width=14)
        self.seq_btn.pack(side=tk.RIGHT)
        self.seq_text = scrolledtext.ScrolledText(self.below, height=4, wrap=tk.CHAR,
                                                  font=("Courier", 9))
        # not packed here: _below_refresh() decides which pane goes where
        self._build_auto_table()

        # ----- Advanced base-call parameters live in a menu dialog ----
        # (no right-hand knob panel → more room for plots, no flicker)
        self._build_settings_dialog()

    def _build_view_bar(self, parent):
        """Second toolbar row: what is on screen and how it is scaled."""
        bar = ttk.Frame(parent)
        bar.pack(fill=tk.X, padx=4, pady=(2, 0))
        self.view_bar = bar

        ttk.Label(bar, text="Graphs").pack(side=tk.LEFT)
        spin = tk.Spinbox(bar, from_=1, to=MAX_GRAPHS, textvariable=self.n_graphs,
                          width=3, command=self.redraw, highlightthickness=0)
        spin.pack(side=tk.LEFT, padx=(2, 4))
        spin.bind("<Return>", lambda e: (self.redraw(), "break")[1])
        spin.bind("<FocusOut>", lambda e: self.redraw())
        self._page_lbl = ttk.Label(bar, text="", foreground="#1E4E8C")
        self._page_lbl.pack(side=tk.LEFT, padx=(0, 10))

        ttk.Separator(bar, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=4)
        ttk.Label(bar, text="Zoom").pack(side=tk.LEFT)
        ttk.Button(bar, text="−", width=2,
                   command=lambda: self._wheel_zoom("x", -1)).pack(side=tk.LEFT, padx=1)
        ttk.Button(bar, text="+", width=2,
                   command=lambda: self._wheel_zoom("x", 1)).pack(side=tk.LEFT, padx=1)
        ttk.Label(bar, text="Scans").pack(side=tk.LEFT, padx=(8, 2))
        self.scan_from = tk.StringVar()
        self.scan_to = tk.StringVar()
        self._scan_from_entry = ttk.Entry(bar, textvariable=self.scan_from, width=6)
        self._scan_from_entry.pack(side=tk.LEFT)
        ttk.Label(bar, text="–").pack(side=tk.LEFT)
        self._scan_to_entry = ttk.Entry(bar, textvariable=self.scan_to, width=6)
        self._scan_to_entry.pack(side=tk.LEFT)
        for e in (self._scan_from_entry, self._scan_to_entry):
            e.bind("<Return>", self._scan_entry_go)
            e.bind("<KP_Enter>", self._scan_entry_go)
        ttk.Button(bar, text="Go", width=3,
                   command=self._scan_entry_go).pack(side=tk.LEFT, padx=(2, 0))

        ttk.Separator(bar, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=8)
        ttk.Label(bar, text="Y scale").pack(side=tk.LEFT)
        ycb = ttk.Combobox(bar, state="readonly", width=12,
                           values=["Shared", "Each graph"])
        ycb.set("Each graph" if self.y_scale.get() == "each" else "Shared")
        ycb.pack(side=tk.LEFT, padx=(2, 6))

        def _y_pick(_e=None):
            self.y_scale.set("each" if ycb.get() == "Each graph" else "shared")
            self._view_y[:] = [0.0, 1.0]
            self._view_option_changed()
        ycb.bind("<<ComboboxSelected>>", _y_pick)
        self._y_scale_cb = ycb
        self.y_scale.trace_add("write", lambda *_a: self._sync_y_scale_cb())
        ttk.Checkbutton(bar, text="Fit Y to view", variable=self.fit_y,
                        command=self._view_option_changed).pack(side=tk.LEFT, padx=4)
        ttk.Checkbutton(bar, text="Overlay runs", variable=self.overlay_runs,
                        command=lambda: self._view_option_changed(redraw=True)
                        ).pack(side=tk.LEFT, padx=4)
        ttk.Button(bar, text="⟲ Reset view",
                   command=self._reset_zoom).pack(side=tk.RIGHT, padx=4)
        # Plate map and Align sit in the channel bar, which has room to spare;
        # the view bar is full at ordinary window widths.
        top = getattr(self, "_chan_bar", bar)
        ttk.Button(top, text="▦ Plate map",
                   command=self.show_plate_map).pack(side=tk.RIGHT, padx=(4, 10))
        amb = ttk.Menubutton(top, text="⇄ Align")
        amb["menu"] = self._alignment_menu(amb)
        amb.pack(side=tk.RIGHT, padx=4)
        ToolTip(spin, lambda: f"Graphs stacked at once (1–{MAX_GRAPHS})")
        ToolTip(ycb, lambda: "Shared: one height scale for all graphs (compare wells).\n"
                             "Each graph: every graph scaled to its own peaks.")
        ToolTip(self._scan_from_entry, lambda: "Type a scan range and press Enter")
        hint = ttk.Label(bar, text=" ⓘ ", foreground="#1E4E8C", cursor="question_arrow")
        hint.pack(side=tk.LEFT, padx=(2, 0))
        ToolTip(hint, lambda: (
            "On the graphs:\n"
            "  wheel = zoom at the cursor    Ctrl+wheel = zoom Y    Shift+wheel = pan\n"
            "  drag = zoom to a box (Ctrl+drag while picking peaks)\n"
            "  middle-drag = pan    double-click / right-click / Home = reset"))

    def _sync_y_scale_cb(self):
        cb = getattr(self, "_y_scale_cb", None)
        if cb is not None:
            cb.set("Each graph" if self.y_scale.get() == "each" else "Shared")

    def _alignment_menu(self, parent):
        m = tk.Menu(parent, tearoff=0)
        m.add_command(label="Align shown graphs to the top one (scans in view)",
                      command=self.align_to_top)
        m.add_command(label="Align on internal standard", command=self.align_on_standard)
        m.add_separator()
        m.add_command(label="Reset this graph (Alt+Home)",
                      command=lambda: self._reset_alignment("hover"))
        m.add_command(label="Reset all alignments", command=lambda: self._reset_alignment("all"))
        return m

    def _build_mode_banner(self, parent):
        """A coloured strip saying the mouse now records peaks, with the
        actions that only make sense in that mode.  Hidden while viewing."""
        ban = tk.Frame(parent, bg="#FFE9A8", highlightthickness=0)
        self._mode_banner = ban
        self._mode_banner_lbl = tk.Label(ban, bg="#FFE9A8", fg="#5A3E00", anchor="w",
                                         font=("DejaVu Sans", 10, "bold"))
        self._mode_banner_lbl.pack(side=tk.LEFT, padx=8, pady=3)
        self._mode_banner_hint = tk.Label(ban, bg="#FFE9A8", fg="#5A3E00", anchor="w")
        self._mode_banner_hint.pack(side=tk.LEFT, padx=4)
        ttk.Button(ban, text="Exit (Esc)",
                   command=self.exit_genotyping_picking).pack(side=tk.RIGHT, padx=6, pady=2)
        ttk.Button(ban, text="Clear shown wells",
                   command=self._gen_clear_shown).pack(side=tk.RIGHT, padx=2, pady=2)
        self._mode_banner_is = ttk.Button(ban, text="Mark IS in shown wells",
                                          command=self._gen_mark_std_batch)
        self._mode_banner_is.pack(side=tk.RIGHT, padx=2, pady=2)

    def _update_mode_banner(self):
        ban = getattr(self, "_mode_banner", None)
        if ban is None:
            return
        if not getattr(self, "genotyping_active", False):
            ban.pack_forget()
            return
        area = getattr(self, "area_mode", False)
        colour = "#D9ECFF" if area else "#FFE9A8"
        fg = "#0F3A6E" if area else "#5A3E00"
        title = ("● Peak area by drag" if area else "● Manual peak picking")
        hint = ("drag between two points to measure  ·  Ctrl+drag zooms  ·  Ctrl+Z undoes"
                if area else
                "click a peak to record it  ·  R+click removes  ·  Ctrl+drag zooms  ·  Ctrl+Z undoes")
        for w in (ban, self._mode_banner_lbl, self._mode_banner_hint):
            w.configure(bg=colour)
        self._mode_banner_lbl.configure(text=title, fg=fg)
        self._mode_banner_hint.configure(text=hint, fg=fg)
        if area:
            self._mode_banner_is.pack_forget()
        elif not self._mode_banner_is.winfo_ismapped():
            self._mode_banner_is.pack(side=tk.RIGHT, padx=2, pady=2)
        if not ban.winfo_ismapped():
            ban.pack(fill=tk.X, padx=4, pady=(2, 0), before=self.vsplit)

    # ----------------------------------------------------------- plate map
    PLATE_COLOURS = {
        "hom1": ("#3B82C4", "1"), "hom2": ("#E08A2E", "2"), "het": ("#3BA55C", "H"),
        "nocall": ("#9AA0A6", "–"), "other": ("#8E6CC7", "?"),
        "unscored": ("#FFFFFF", ""), "missing": ("#ECEDEF", ""),
    }

    @staticmethod
    def _call_class(call):
        c = str(call or "").lower().replace("_", "-")
        if not c:
            return "unscored"
        if "no" in c and "call" in c:
            return "nocall"
        if c.startswith("het"):
            return "het"
        if c in ("hom-1", "hom1"):
            return "hom1"
        if c in ("hom-2", "hom2"):
            return "hom2"
        return "other"

    @staticmethod
    def _well_rc(name):
        """('A', 1) for 'A01' / 'a1' / '..._A01'; None if no well id."""
        m = re.fullmatch(r"([A-Pa-p])(\d{1,2})", name)
        if m:                                   # a bare well id, 96 or 384
            return m.group(1).upper(), int(m.group(2))
        m = WELL_RE.search(name)                # a well inside a longer name
        if m:
            w = m.group(1)
            return w[0].upper(), int(w[1:])
        return None

    def _plate_runs(self):
        runs = []
        for p in list(self.selected or []) + list(self.files or []):
            r = Path(p).parent.name
            if r not in runs:
                runs.append(r)
        return runs

    def _plate_files(self, run):
        """{(row, col): path} for every trace of *run* in the loaded folders."""
        cache = getattr(self, "_plate_files_cache", {})
        if run not in cache:
            files = [p for p in discover_files(self.folders) if p.parent.name == run]
            out = {}
            for p in files:
                rc = self._well_rc(p.stem)
                if rc and rc not in out:
                    out[rc] = p
            cache[run] = out
            self._plate_files_cache = cache
        return cache[run]

    def show_plate_map(self):
        """A 96-well (or 384-well) map of one run: each well coloured by its
        call, the wells on screen outlined; click to show wells."""
        win = getattr(self, "_plate_win", None)
        try:
            if win is not None and win.winfo_exists():
                win.deiconify()
                win.lift()
                self._plate_map_refresh(force=True)
                return
        except tk.TclError:
            pass
        runs = self._plate_runs()
        if not runs:
            messagebox.showinfo("Plate map", "Add a data folder first.")
            return
        win = tk.Toplevel(self)
        self._plate_win = win
        win.title("Plate map")
        win.transient(self)
        top = ttk.Frame(win, padding=(8, 8, 8, 2))
        top.pack(fill=tk.X)
        ttk.Label(top, text="Run").pack(side=tk.LEFT)
        self._plate_run = tk.StringVar(value=runs[0])
        cb = ttk.Combobox(top, textvariable=self._plate_run, values=runs,
                          state="readonly", width=40)
        cb.pack(side=tk.LEFT, padx=6)
        cb.bind("<<ComboboxSelected>>", lambda e: self._plate_map_refresh(force=True))
        self._plate_canvas = tk.Canvas(win, bg="#FFFFFF", highlightthickness=0)
        self._plate_canvas.pack(fill=tk.BOTH, expand=True, padx=8, pady=4)
        legend = ttk.Frame(win, padding=(8, 0, 8, 2))
        legend.pack(fill=tk.X)
        for key, text in (("hom1", "hom-1"), ("hom2", "hom-2"), ("het", "het"),
                          ("nocall", "no-call"), ("unscored", "not scored"),
                          ("missing", "no file")):
            col, _ = self.PLATE_COLOURS[key]
            sw = tk.Canvas(legend, width=14, height=14, highlightthickness=0, bg=win.cget("bg"))
            sw.create_rectangle(1, 1, 13, 13, fill=col, outline="#888")
            sw.pack(side=tk.LEFT, padx=(6, 2))
            ttk.Label(legend, text=text).pack(side=tk.LEFT)
        hint = ttk.Label(win, foreground="#555", padding=(8, 0, 8, 8), justify=tk.LEFT, text=(
            "Click a well to show it · Ctrl-click to add/remove · click a row letter "
            "or column number for the whole row/column · red outline = on screen"))
        hint.pack(fill=tk.X)
        hint.bind("<Configure>", lambda e: hint.configure(wraplength=max(200, e.width - 16)))
        self._plate_canvas.bind("<Configure>", lambda e: self._plate_map_refresh(force=True))
        self._plate_canvas.bind("<Button-1>", lambda e: self._plate_click(e, add=False))
        self._plate_canvas.bind("<Control-Button-1>", lambda e: self._plate_click(e, add=True))
        win.geometry("640x420")
        self._plate_map_refresh(force=True)

    def _plate_geometry(self, files):
        rows = "ABCDEFGH"
        ncol = 12
        if any(r > "H" for r, _c in files) or any(c > 12 for _r, c in files):
            rows, ncol = "ABCDEFGHIJKLMNOP", 24
        cv = self._plate_canvas
        w = max(200, cv.winfo_width())
        h = max(150, cv.winfo_height())
        cell = max(12, min((w - 30) / ncol, (h - 24) / len(rows)))
        return rows, ncol, cell, 26, 20

    def _plate_map_refresh(self, force=False):
        win = getattr(self, "_plate_win", None)
        try:
            if win is None or not win.winfo_exists():
                return
        except tk.TclError:
            return
        run = self._plate_run.get()
        files = self._plate_files(run)
        calls = {}
        for r in (getattr(self, "_geno_rows", None) or []):
            if r.get("run") == run and r.get("well"):
                rc = self._well_rc(str(r["well"]))
                if rc:
                    calls[rc] = r.get("call", "")
        shown = {self._well_rc(Path(p).stem) for p in (self.selected or [])[: self._graph_count() * 4]
                 if Path(p).parent.name == run}
        cv = self._plate_canvas
        cv.delete("all")
        rows, ncol, cell, x0, y0 = self._plate_geometry(files)
        self._plate_layout = (rows, ncol, cell, x0, y0)
        small = cell < 26
        for c in range(1, ncol + 1):
            cv.create_text(x0 + (c - 0.5) * cell, y0 / 2, text=str(c),
                           font=("DejaVu Sans", 7 if small else 9), tags=("hdr", f"col{c}"))
        for ri, r in enumerate(rows):
            cv.create_text(x0 / 2, y0 + (ri + 0.5) * cell, text=r,
                           font=("DejaVu Sans", 7 if small else 9, "bold"), tags=("hdr", f"row{r}"))
            for c in range(1, ncol + 1):
                rc = (r, c)
                if rc not in files:
                    key = "missing"
                else:
                    key = self._call_class(calls.get(rc, ""))
                fill, mark = self.PLATE_COLOURS[key]
                xa, ya = x0 + (c - 1) * cell + 1, y0 + ri * cell + 1
                on = rc in shown
                cv.create_rectangle(xa, ya, xa + cell - 2, ya + cell - 2, fill=fill,
                                    outline="#D62828" if on else "#B8BDC4",
                                    width=3 if on else 1)
                if mark and not small:
                    cv.create_text(xa + cell / 2 - 1, ya + cell / 2 - 1, text=mark,
                                   fill="#FFFFFF" if key != "unscored" else "#333",
                                   font=("DejaVu Sans", 9, "bold"))

    def _plate_click(self, event, add=False):
        layout = getattr(self, "_plate_layout", None)
        if not layout:
            return "break"
        rows, ncol, cell, x0, y0 = layout
        run = self._plate_run.get()
        files = self._plate_files(run)
        col = int((event.x - x0) // cell) + 1
        ri = int((event.y - y0) // cell)
        targets = []
        if event.x < x0 and 0 <= ri < len(rows):              # row header
            targets = [files[(rows[ri], c)] for c in range(1, ncol + 1) if (rows[ri], c) in files]
        elif event.y < y0 and 1 <= col <= ncol:                 # column header
            targets = [files[(r, col)] for r in rows if (r, col) in files]
        elif 0 <= ri < len(rows) and 1 <= col <= ncol and (rows[ri], col) in files:
            targets = [files[(rows[ri], col)]]
        if not targets:
            return "break"
        if self.overlay_runs.get():
            # the same wells in every run of this plate
            plate = self._plate_of(targets[0])
            wells = {self._well_label(p) for p in targets}
            targets = [p for p in self.files
                       if self._plate_of(p) == plate and self._well_label(p) in wells] or targets
        if add:
            sel = list(self.selected or [])
            keys = {str(Path(p).resolve()) for p in sel}
            for p in targets:
                k = str(Path(p).resolve())
                if k in keys:
                    sel = [s for s in sel if str(Path(s).resolve()) != k]
                else:
                    sel.append(p)
            targets = sel
        if len(targets) > 1 and not add:
            per = len({self._well_label(p) for p in targets}) if self.overlay_runs.get() else len(targets)
            self.n_graphs.set(min(MAX_GRAPHS, max(self._graph_count(), per)))
        self.selected = list(targets)
        try:
            idx = {str(Path(p).resolve()): i for i, p in enumerate(self.files)}
            self.file_list.selection_clear(0, tk.END)
            for p in self.selected:
                i = idx.get(str(Path(p).resolve()))
                if i is not None:
                    self.file_list.selection_set(i)
                    self.file_list.see(i)
        except tk.TclError:
            pass
        self.redraw()
        return "break"

    def _update_page_label(self, shown):
        """'Wells 5–8 of 192' -- and say so when the selection is bigger than
        what fits, instead of silently drawing only the first few."""
        lbl = getattr(self, "_page_lbl", None)
        if lbl is None:
            return
        sel = list(self.selected or [])
        if not sel:
            lbl.config(text="")
            return
        unit = "wells" if self.overlay_runs.get() and not getattr(
            self, "genotyping_active", False) else "graphs"
        if len(sel) > len(shown) and unit == "graphs":
            lbl.config(text=f"showing {len(shown)} of {len(sel)} selected — "
                            "raise Graphs or page with ↓/↑", foreground="#A33")
            return
        wells = [self._well_label(p) for p in shown]
        uniq = sorted(set(wells), key=lambda w: (w[:1], w[1:]))
        runs = len({Path(p).parent.name for p in shown})
        if not uniq:
            lbl.config(text="")
            return
        span = uniq[0] if len(uniq) == 1 else f"{uniq[0]}–{uniq[-1]}"
        extra = f" · {runs} runs" if runs > 1 else ""
        lbl.config(text=f"{span}{extra} · {len(shown)} of {len(self.files)} files",
                   foreground="#1E4E8C")

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
        if not hasattr(self, "base_order_var"):        # keep the one traced
            self.base_order_var = tk.StringVar(value=self.settings.base_order)
        ttk.Label(f, text="Base order (instrument)").pack(anchor=tk.W)
        ttk.Combobox(f, textvariable=self.base_order_var, state="readonly",
                     values=["ACTG", "TGCA", "GATC", "CTAG"], width=12).pack(anchor=tk.W)
        ttk.Label(f, text="MegaBACE default is ACTG (Ch1=A, Ch2=C, Ch3=T, Ch4=G); "
                          "DYEnamic runs are TGCA. ABI uses its own FWO order.",
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
        self.gauss_seg = tk.IntVar(value=384)
        self.gauss_reg = tk.DoubleVar(value=0.05)
        check(f, "Gaussian reconstruction", self.gauss_en)
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
        # The genotyping table shows the library's calls for the runs in
        # play, so refill it when that set of runs changes (a folder added at
        # start-up or later, or cleared) -- not on every filter keystroke.
        runs = frozenset(f.parent.name for f in self.files)
        if runs != getattr(self, "_geno_runs_shown", None):
            self._geno_runs_shown = runs
            self._plate_files_cache = {}        # folders changed: re-scan
            self._sync_geno_table()

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
        n = self._graph_count()
        if self.overlay_runs.get() and not getattr(self, "genotyping_active", False):
            idxs = self._overlay_page(delta, n)
        elif m <= n:
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

    def _overlay_page(self, delta, n):
        """List indices of the next/previous n wells (every run of each), for
        paging while "Overlay runs" draws one graph per well of a plate."""
        m = len(self.files)
        key = lambda i: (self._plate_of(self.files[i]), self._well_label(self.files[i]))
        byp = {str(p.resolve()): i for i, p in enumerate(self.files)}
        shown = sorted(byp[str(p.resolve())] for p in (self.selected or [])
                       if str(p.resolve()) in byp)
        if delta > 0:
            start = (shown[-1] + 1) if shown else 0
            if start >= m:
                start = 0
        else:
            j = (shown[0] - 1) if shown else m - 1
            if j < 0:
                j = m - 1
            wells = []
            while j >= 0:
                w = key(j)
                if w not in wells:
                    if len(wells) == n:
                        break
                    wells.append(w)
                j -= 1
            start = j + 1
        idxs, wells = [], []
        for i in range(start, m):
            w = key(i)
            if w not in wells:
                if len(wells) == n:
                    break
                wells.append(w)
            idxs.append(i)
        return idxs

    def _page_key(self, event):
        if int(getattr(event, "state", 0) or 0) & 0x0008:     # Alt: alignment nudge
            return None
        keys = event.keysym
        if keys in ("Down", "Right", "Page_Down"):
            self._page_by(1)
        else:
            self._page_by(-1)
        return "break"

    def _page_key_root(self, event):
        if int(getattr(event, "state", 0) or 0) & 0x0008:     # Alt: alignment nudge
            return None
        keys = event.keysym
        if keys not in ("Down", "Up", "Left", "Right", "Page_Down", "Page_Up"):
            return None
        w = self.focus_get()
        cls = w.winfo_class() if w is not None else ""
        # A text field owns all four arrows: they move the caret, not the page.
        if cls in ("Entry", "TEntry", "Text", "Spinbox", "TSpinbox",
                   "Combobox", "TCombobox"):
            return None
        # Everything else -- the file list, a graph canvas, and the pick and
        # result tables -- pages the graphs.  The tables are the point: while
        # picking, the pick table holds focus, and its own Up/Down binding ate
        # the key before this toplevel one ever saw it, so the arrows felt dead
        # exactly where the operator was working.
        if keys in ("Down", "Right", "Page_Down"):
            self._page_by(1)
        else:
            self._page_by(-1)
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

    def toggle_tour(self):
        self.touring = not self.touring
        if self.touring:
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
        n = self._graph_count()
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
        self._refresh_channel_buttons()
        self.redraw()

    def _ensure_doc(self, path: Path) -> TraceDocument:
        key = str(path.resolve())
        if key not in self.docs:
            prev_status = self.status_var.get()
            self.status_var.set(f"Loading {path.name}…")
            self.update_idletasks()
            try:
                order, src = self._order_for(path)
                doc = load_trace(path, base_order=order)
                declared = declared_dye_order(getattr(doc, "footer_fields", None))
                if src == "default" and declared and declared != order:
                    # The file says how its channels are dyed; believe it
                    # over the general setting (a run-specific choice still
                    # wins -- see Comments ▸ Dye order of this run…).
                    order, src = declared, "file"
                    doc = load_trace(path, base_order=order)
                elif src == "default" and declared:
                    src = "file"
                doc.dye_order = order
                doc.dye_order_source = src
                self.docs[key] = doc
                self._restore_alignment(path)
                if self.status_var.get().startswith("Loading "):
                    self.status_var.set(prev_status)
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
                        run_basecall(doc, self._settings_for(doc, self.settings))
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

    def clear_basecalls(self):
        """Undo a base call: drop sequence/peaks/qualities so a run that was
        mistakenly call-ed (e.g. genotyping/fragment data) shows clean traces
        again.  Only the in-memory documents are touched; data files stay
        as they are."""
        if not self.selected:
            messagebox.showinfo("Clear base calls",
                                "Select one or more wells in the list first.")
            return
        cleared = 0
        for path in self.selected:
            doc = self.docs.get(str(path.resolve()))
            if doc is not None and doc.sequence:
                doc.sequence = ""
                doc.peak_positions = []
                doc.qualities = []
                doc.settings_used = None
                cleared += 1
        self.status_var.set(
            f"Cleared base calls for {cleared} well(s)" if cleared
            else "No base calls to clear on the selected wells")
        self.redraw()
        self._show_sequence()

    def toggle_genotyping_picking(self):
        """Main-window manual peak picking: the viewer stays exactly as it is
        (zoom bars, µA, Channels), but clicks on a subplot pick peaks for that
        well and the sequence pane shows the picked-peaks table instead.
        Calling it again (or Exit peak picking) returns to the viewer."""
        if self.genotyping_active and not self.area_mode:
            self.exit_genotyping_picking()
            return
        self.enter_genotyping_picking()

    def enter_genotyping_picking(self):
        if not self.selected:
            messagebox.showinfo("Genotyping",
                                "Select the wells to genotype in the list first.")
            return
        if self.genotyping_active:
            self.exit_genotyping_picking()     # switching between the two modes
        self.area_mode = False
        self._mode_pick.set(True)
        self._mode_area.set(False)
        self._build_pick_table()
        self._pick_cid = self.canvas.mpl_connect("button_press_event",
                                                 self._on_gen_pick)
        self._gen_motion_cid = self.canvas.mpl_connect(
            "motion_notify_event", self._on_gen_motion)
        self.genotyping_active = True
        self._below_refresh()          # the pick table takes the slot
        self.status_var.set("Manual peak picking — click each peak you want "
                            "(Ctrl+drag zooms), page on with ↓ / ↑ and pick the "
                            "next wells. R+click removes a pick; Esc leaves.")
        self.redraw()

    def toggle_area_picking(self):
        """Drag-to-measure, in the same window and the same mode plumbing as
        manual peak picking.  Holding the left button and dragging between two
        points measures the area above the straight line between them, so there
        is no second window to keep in sync with this one."""
        if self.genotyping_active and self.area_mode:
            self.exit_genotyping_picking()
            return
        self.enter_area_picking()

    def enter_area_picking(self):
        if not self.selected:
            messagebox.showinfo("Area measure",
                                "Select the wells to measure in the list first.")
            return
        if self.genotyping_active:
            self.exit_genotyping_picking()
        self.area_mode = True
        self._mode_area.set(True)
        self._mode_pick.set(False)
        self._build_pick_table(area=True)
        self._pick_cid = self.canvas.mpl_connect("button_press_event",
                                                 self._on_gen_pick)
        self._gen_motion_cid = self.canvas.mpl_connect(
            "motion_notify_event", self._on_gen_motion)
        self._area_release_cid = self.canvas.mpl_connect(
            "button_release_event", self._on_area_release)
        self.genotyping_active = True
        self._below_refresh()          # the area table takes the slot
        self.status_var.set("Drag to measure — hold the left button and drag "
                            "between the two points the baseline should pass "
                            "through, then release. Drag right-to-left works "
                            "too; page the batch as usual.")
        self.redraw()

    def exit_genotyping_picking(self):
        if not self.genotyping_active:
            return
        for attr in ("_pick_cid", "_gen_motion_cid", "_area_release_cid"):
            cid = getattr(self, attr, None)
            if cid is not None:
                try:
                    self.canvas.mpl_disconnect(cid)
                except Exception:
                    pass
                setattr(self, attr, None)
        self._gen_cursors = []
        self._area_cursors = []
        self._area_drag = None
        if self.pick_table is not None:
            try:
                self.pick_table.destroy()
            except Exception:
                pass
            self.pick_table = None
            self.pick_tree = None
        self.genotyping_active = False
        self.area_mode = False
        self._mode_pick.set(False)
        self._mode_area.set(False)
        self._below_refresh()          # hand the slot back to the table/sequence
        self.status_var.set("Back to the trace viewer.")
        self.redraw()

    def _build_pick_table(self, area=False):
        if self.pick_table is not None:
            # rebuilds used to leave the previous table packed below the new
            # one, so every mark-all or semi-automatic pass stacked a second
            # copy of the table into the window
            try:
                self.pick_table.destroy()
            except Exception:
                pass
            self.pick_table = None
            self.pick_tree = None
        tblf = ttk.Frame(self.below)
        tblf.pack(fill=tk.X, padx=4, pady=2)
        if area:
            cols = ("run", "well", "ch", "start", "stop", "mid", "h V", "area")
            widths = {"run": 90, "well": 58, "ch": 34, "start": 48, "stop": 48,
                      "mid": 50, "h V": 62, "area": 64}
            left_anchor = ("run", "well")
            stretch = ("run", "well")
        else:
            cols = ("well", "scan", "duplex", "ch", "kind", "h V", "area",
                    "frac", "MF")
            widths = {"well": 60, "scan": 54, "duplex": 58, "ch": 34,
                      "kind": 62, "h V": 66, "area": 66, "frac": 50,
                      "MF": 52}
            left_anchor = ("well", "kind")
            stretch = ("scan", "kind")
        self.pick_tree = ttk.Treeview(tblf, columns=cols, show="headings",
                                      height=4)
        for c in cols:
            self.pick_tree.heading(c, text=c)
            self.pick_tree.column(c, width=widths[c],
                                  anchor="w" if c in left_anchor else "e",
                                  stretch=(c in stretch))
        vs = ttk.Scrollbar(tblf, orient=tk.VERTICAL, command=self.pick_tree.yview)
        self.pick_tree.configure(yscrollcommand=vs.set)
        self.pick_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vs.pack(side=tk.RIGHT, fill=tk.Y)
        # The pick table takes focus when the operator clicks a peak; give it
        # the same arrow paging as the file list so the arrows keep working.
        for seq in ("<Up>", "<Down>", "<Left>", "<Right>",
                    "<Page_Up>", "<Page_Down>"):
            self.pick_tree.bind(seq, self._page_key)
        self.pick_table = tblf

    def _ensure_area_picker(self, path):
        """One DragAreaPicker per displayed well, built on first use."""
        from genotyping import DragAreaPicker
        key = str(path.resolve())
        pk = self._area_pickers.get(key)
        if pk is None:
            doc = self._ensure_doc(path)
            col2chan = self._col_to_chan(self._doc_order(doc))
            pk = DragAreaPicker(doc, path, colors=self._theme_colors(),
                                base_order=self._doc_order(doc),
                                theme_mode=self._theme_mode(),
                                show=lambda c: self.chan_show[
                                    col2chan.get(c, c)].get(),
                                run_name=path.parent.name)
            self._area_pickers[key] = pk
        return pk

    def _run_orders(self):
        """Dye orders chosen per run (folder name -> order), from the prefs."""
        cache = getattr(self, "_run_orders_cache", None)
        if cache is None:
            try:
                cache = dict(mark_library.load_prefs().get("run_dye_orders", {}))
            except Exception:                              # noqa: BLE001
                cache = {}
            self._run_orders_cache = cache
        return cache

    def _order_for(self, path):
        """(dye order, where it came from) for the run *path* belongs to.

        "run": chosen for this run; "default": the general setting.  _ensure_doc
        upgrades "default" to the order the file declares, when it declares one.
        """
        run = Path(path).parent.name
        chosen = self._run_orders().get(run)
        if chosen:
            return chosen.upper(), "run"
        return (self.base_order_var.get() or "ACTG").upper(), "default"

    def _doc_order(self, doc):
        return (getattr(doc, "dye_order", None) or self.base_order_var.get()
                or "ACTG").upper()

    def _settings_for(self, doc, settings):
        """*settings* with the dye order of *doc*'s own run."""
        import dataclasses
        order = self._doc_order(doc)
        if (settings.base_order or "").upper() == order:
            return settings
        return dataclasses.replace(settings, base_order=order)

    def _display_order(self):
        """The order the channel bar is labelled with: the first well on
        screen, else the general setting."""
        for p in (self.selected or [])[: self._graph_count()]:
            doc = self.docs.get(str(Path(p).resolve()))
            if doc is not None:
                return self._doc_order(doc)
        return (self.base_order_var.get() or "ACTG").upper()

    def set_run_dye_order(self):
        """Comments ▸ Dye order of this run…: fix the order for one run."""
        paths = list(self.selected or [])
        if not paths:
            messagebox.showinfo("Dye order of this run", "Select a well of the run first.")
            return
        run = Path(paths[0]).parent.name
        doc = self.docs.get(str(Path(paths[0]).resolve()))
        declared = declared_dye_order(getattr(doc, "footer_fields", None)) if doc else ""
        win = tk.Toplevel(self)
        win.title("Dye order of this run")
        win.transient(self)
        body = ttk.Frame(win, padding=12)
        body.pack(fill=tk.BOTH, expand=True)
        ttk.Label(body, text=f"Run: {run}", font=("", 10, "bold")).pack(anchor=tk.W)
        ttk.Label(body, wraplength=360, justify=tk.LEFT, text=(
            (f"The file declares {declared}." if declared else
             "The file does not declare a dye order (normal for genotyping runs).")
            + f"\nGeneral setting: {(self.base_order_var.get() or 'ACTG').upper()}."
        )).pack(anchor=tk.W, pady=(4, 8))
        choices = ["Automatic (file, else general setting)", "ACTG", "TGCA", "GATC", "CTAG"]
        cur = self._run_orders().get(run)
        var = tk.StringVar(value=cur if cur in choices else choices[0])
        ttk.Combobox(body, textvariable=var, values=choices, state="readonly",
                     width=36).pack(anchor=tk.W)

        def save():
            v = var.get()
            orders = self._run_orders()
            if v == choices[0]:
                orders.pop(run, None)
            else:
                orders[run] = v
            try:
                prefs = mark_library.load_prefs()
                prefs["run_dye_orders"] = dict(orders)
                mark_library.save_prefs(prefs)
            except Exception:                              # noqa: BLE001
                pass
            # drop this run's traces (and its pickers, unless they hold picks)
            for key in [k for k in self.docs if Path(k).parent.name == run]:
                self.docs.pop(key, None)
            for store in (self._gen_pickers, self._area_pickers):
                for key in [k for k, pk in store.items()
                            if Path(k).parent.name == run and not getattr(pk, "records", None)]:
                    store.pop(key, None)
            win.destroy()
            self.redraw()
            self.status_var.set(f"{run}: dye order {v}")
        btns = ttk.Frame(body)
        btns.pack(fill=tk.X, pady=(10, 0))
        ttk.Button(btns, text="Cancel", command=win.destroy).pack(side=tk.RIGHT)
        ttk.Button(btns, text="Save", command=save).pack(side=tk.RIGHT, padx=6)
        self._grab_soon(win)

    def _gen_paths(self):
        n = self._graph_count()
        sel = [Path(p) for p in (self.selected or [])]
        return sel[:n]

    def _ensure_picker(self, path):
        from genotyping import PeakPicker
        key = str(path.resolve())
        pk = self._gen_pickers.get(key)
        if pk is None:
            doc = self._ensure_doc(path)
            # the picker is handed acgt columns, the checkboxes are indexed by
            # physical channel -- translate, or a non-ACTG run hides one trace
            # and picks from another
            col2chan = self._col_to_chan(self._doc_order(doc))
            pk = PeakPicker(doc, path, colors=self._theme_colors(),
                            base_order=self._doc_order(doc),
                            theme_mode=self._theme_mode(),
                            show=lambda c: self.chan_show[
                                col2chan.get(c, c)].get(),
                            include_sh=self._gen_sh.get(),
                            show_d2=self._gen_d2.get())
            self._gen_pickers[key] = pk
            self._library_apply(path, pk)
        return pk

    def _library_wells(self, run):
        """One run's saved entries: read from disk the first time they are
        wanted, kept in memory after that."""
        wells = self._library_cache.get(run)
        if wells is None:
            wells = mark_library.load_run(run)
            self._library_cache[run] = wells
        return wells

    def _library_apply(self, path, pk):
        """Hand a picker back the marks this well was saved with, so a restart
        does not cost the operator the plate they already marked.  A library
        problem must never stop a well from opening."""
        try:
            entry = self._library_wells(path.parent.name).get(path.stem)
            if entry:
                mark_library.apply_picker(pk, entry, default_col=self._is_col(path))
        except Exception:
            pass

    def _is_col(self, path=None):
        """The ACGT column the internal standard lives on (in *path*'s run),
        or None if the settings are not readable yet."""
        try:
            import genotyping
            order = (self._order_for(path)[0] if path is not None
                     else self.base_order_var.get())
            doc = self.docs.get(str(Path(path).resolve())) if path is not None else None
            if doc is not None:
                order = self._doc_order(doc)
            return genotyping.acgt_index_for_channel(
                order, int(self._auto_is_channel.get()))
        except Exception:                              # noqa: BLE001
            return None

    def _library_store(self, *pairs):
        """Save marked wells, one file write per run.

        Only wells that carry a standard are stored: peaks picked but not yet
        marked a standard stay in the session, because the library is a record
        of answers, not of half-finished clicks.
        """
        dirty = {}
        items = []
        # Accept either (path, pk, ...) flat or ((path, pk), ...) tuples
        if len(pairs) == 2 and not isinstance(pairs[0], tuple):
            items = [(pairs[0], pairs[1])]
        else:
            for p in pairs:
                if isinstance(p, tuple) and len(p) >= 2:
                    items.append((p[0], p[1]))
                # else skip malformed
        for path, pk in items:
            if not pk.std:
                continue
            run = Path(path).parent.name
            wells = self._library_wells(run)
            entry = dict(wells.get(Path(path).stem) or {})
            entry.update(mark_library.entry_from_picker(pk))
            wells[Path(path).stem] = entry
            dirty[run] = wells
        for run, wells in dirty.items():
            try:
                mark_library.save_run(run, wells)
            except Exception as exc:                      # noqa: BLE001
                # A lost write means the operator's marks vanish on restart
                # and nothing on screen says so, so it goes to the log.
                print(f"library: could not save {run}: {exc}", file=sys.stderr)

    def _library_store_results(self, rows, label_source="engine"):
        """Save what auto-genotype measured and called for each well -- the
        labelled half of the machine-learning export.  ``label_source`` says
        whose call it is, so an operator's correction is never confused with
        the engine's guess."""
        dirty = {}
        for row in rows:
            run, well = row.get("run"), row.get("well")
            if not run or not well:
                continue
            wells = self._library_wells(run)
            entry = dict(wells.get(well) or {})
            entry.update(mark_library.entry_from_row(row, label_source))
            wells[well] = entry
            dirty[run] = wells
        for run, wells in dirty.items():
            try:
                mark_library.save_run(run, wells)
            except Exception as exc:                      # noqa: BLE001
                # A lost write means the operator's marks vanish on restart
                # and nothing on screen says so, so it goes to the log.
                print(f"library: could not save {run}: {exc}", file=sys.stderr)

    def _gen_axes_hit(self, event):
        """Resolve a mouse event to (index, trace axes, scan, volts).

        matplotlib reports ``event.inaxes`` as the last axes sharing the hit
        point, which is the uA twin when the current overlay is on -- so the
        event's own ``ydata`` would be microamps, not volts.  When display
        coords are available, match against the stacked trace axes and invert
        that axes' own transform instead.  Synthetic events (tests) without
        display coords fall back to the ``inaxes`` identity."""
        ex = getattr(event, "x", None)
        ey = getattr(event, "y", None)
        x = getattr(event, "xdata", None)
        y = getattr(event, "ydata", None)
        ax = None
        if ex is not None and ey is not None:
            for i, candidate in enumerate(self._plot_axes):
                try:
                    inside = candidate.get_window_extent().contains(ex, ey)
                except Exception:
                    inside = False
                if inside:
                    ax = candidate
                    try:
                        x, y = ax.transData.inverted().transform((ex, ey))
                    except Exception:
                        return None
                    break
            if ax is None:
                return None
        else:
            try:
                i = self._plot_axes.index(event.inaxes)
            except (ValueError, AttributeError):
                return None
            ax = self._plot_axes[i]
        if x is None or y is None:
            return None
        try:
            if not (np.isfinite(x) and np.isfinite(y)):
                return None
        except TypeError:
            return None
        return i, ax, x, y

    def _on_remove_key(self, _event=None):
        """R held: arm remove-on-click, unless a text field has focus."""
        w = self.focus_get()
        cls = w.winfo_class() if w is not None else ""
        if cls in ("Entry", "TEntry", "Text", "Spinbox", "TSpinbox",
                   "Combobox", "TCombobox"):
            return None
        self._remove_mod = True
        if getattr(self, "genotyping_active", False) and not self.area_mode:
            self.status_var.set("Remove: click a picked peak to delete it (with "
                                "its +A tag). Release R to go back to picking.")
        return None

    def _on_remove_release(self, _event=None):
        self._remove_mod = False
        return None

    def _on_gen_pick(self, event):
        if not self.genotyping_active:
            return
        if event.inaxes is None or event.xdata is None:
            return
        if getattr(event, "button", 1) != 1:
            return
        ctrl, shift = _event_mods(event)
        if ctrl or shift:
            return              # Ctrl+drag zooms, Shift+drag slides the graph
        hit = self._gen_axes_hit(event)
        if hit is None:
            return
        i, _ax, xs, vs = hit
        paths = self._gen_paths()
        if not (0 <= i < len(paths)):
            return
        path = paths[i]
        if self.area_mode:
            # arm a span on this subplot; _on_area_release commits it
            self._area_drag = {"i": i, "start": int(round(xs)),
                               "stop": int(round(xs))}
            return
        try:
            pk = self._ensure_picker(path)
        except Exception as e:
            self.status_var.set(f"Could not load {path.name}: {e}")
            return
        if self._remove_mod:
            gone = pk.remove_at(xs, col=self._gen_lock_col(), vol=vs)
            if gone is None:
                self.status_var.set(f"No picked peak near that scan on "
                                    f"{path.name} — click closer to a mark.")
                return
            self._gen_active_path = path
            self._library_store(path, pk)
            self.status_var.set(
                f"{path.parent.name}/{path.name} · well {pk.doc.well}: "
                f"removed peak scan {gone['scan']} · {gone['base']}. "
                "Pick the correct peak, then Mark IS in shown wells.")
            self.redraw()
            self._sync_pick_table()
            return
        rec = pk.pick(xs, vol=vs, only_col=self._gen_lock_col())
        if rec is None:
            if getattr(pk, "_reject", None) == "area":
                self.status_var.set(f"That area is already picked on {path.name} "
                                    "— undo it first to pick it again.")
            else:
                self.status_var.set(f"No peak found near that scan on {path.name} — "
                                    "try again closer to a hump.")
            return
        self._gen_active_path = path
        self.status_var.set(f"{path.parent.name}/{path.name} · well {pk.doc.well}: "
                            f"peak scan {rec['scan']} · {rec['base']} "
                            f"height {rec['height']:.3f} V")
        self.redraw()
        self._sync_pick_table()

    def _make_area_markers(self, ax):
        """Live preview artists for the span being dragged on one subplot.

        Persistent artists updated with set_data/set_xdata rather than a full
        figure rebuild, because this runs on every motion event and a redraw of
        eight stacked traces would stutter badly.
        """
        line = ax.plot([], [], color="#B00020", lw=1.2, ls="--", alpha=0.9,
                       zorder=23, visible=False)[0]
        left = ax.axvline(0, color="#B00020", lw=0.7, ls=":", alpha=0.7,
                          zorder=23, visible=False)
        right = ax.axvline(0, color="#B00020", lw=0.7, ls=":", alpha=0.7,
                           zorder=23, visible=False)
        ann = ax.annotate("", xy=(0.5, 0.90), xycoords="axes fraction",
                          ha="center", va="top", fontsize=7, color="#111",
                          zorder=25, visible=False,
                          bbox=dict(boxstyle="round,pad=0.3", fc="#FFF9C4",
                                    ec="#999", lw=0.6, alpha=0.9))
        return (ax, line, left, right, ann)

    def _area_preview(self, event):
        """Follow the pointer with the dashed baseline and a live area readout."""
        from genotyping import DRAG_MIN_SPAN, region_area
        drag = self._area_drag
        if drag is None or not self._area_cursors:
            return
        if getattr(event, "xdata", None) is None:
            return
        drag["stop"] = int(round(event.xdata))
        i = drag["i"]
        if not (0 <= i < len(self._area_cursors)):
            return
        ax, line, left, right, ann = self._area_cursors[i]
        a, b = sorted((drag["start"], drag["stop"]))
        for art in self._area_cursors:
            for k in art[1:]:
                k.set_visible(False)
        if b - a < DRAG_MIN_SPAN:
            self.canvas.draw_idle()
            return
        paths = self._gen_paths()
        if not (0 <= i < len(paths)):
            return
        try:
            pk = self._ensure_area_picker(paths[i])
        except Exception:
            return
        col = pk.best_col(a, b)
        if col is None:
            return
        y = np.asarray(pk.doc.acgt[:, col], dtype=float)
        a = int(np.clip(a, 0, y.size - 1))
        b = int(np.clip(b, a, y.size - 1))
        xs = np.arange(a, b + 1)
        bl = np.linspace(float(y[a]), float(y[b]), b - a + 1)
        line.set_data(xs, bl)
        line.set_visible(True)
        left.set_xdata([a, a])
        right.set_xdata([b, b])
        left.set_visible(True)
        right.set_visible(True)
        res = region_area(y, a, b)
        if res is not None:
            ann.set_text(f"mid {res['midpoint']} · {res['area']:.1f} V·scan")
            ann.set_visible(True)
        self.canvas.draw_idle()

    def _on_area_release(self, event):
        """Commit the dragged span as one measurement."""
        from genotyping import DRAG_MIN_SPAN
        drag = self._area_drag
        if not self.area_mode or drag is None:
            return
        self._area_drag = None
        for art in self._area_cursors:
            for k in art[1:]:
                k.set_visible(False)
        if getattr(event, "button", 1) != 1:
            self.canvas.draw_idle()
            return
        # The release position is the reliable end of the span: a brisk drag
        # can finish before any motion event arrives.
        x = getattr(event, "xdata", None)
        stop = int(round(x)) if x is not None else drag["stop"]
        i = drag["i"]
        paths = self._gen_paths()
        if not (0 <= i < len(paths)):
            self.canvas.draw_idle()
            return
        path = paths[i]
        a, b = sorted((drag["start"], stop))
        if b - a < DRAG_MIN_SPAN:
            self.status_var.set(f"Drag at least {DRAG_MIN_SPAN} scans to "
                                f"measure a peak on {path.name}.")
            self.canvas.draw_idle()
            return
        try:
            pk = self._ensure_area_picker(path)
        except Exception as e:
            self.status_var.set(f"Could not load {path.name}: {e}")
            return
        rec = pk.add(a, b)
        if rec is None:
            # Two quite different reasons nothing came back, and "nothing rises
            # above the baseline" is simply false when every channel is off.
            if not any(v.get() for v in self.chan_show):
                self.status_var.set(
                    "No channel is switched on, so there is nothing to measure "
                    f"on {path.name} — tick at least one in the Channels row.")
            else:
                self.status_var.set(f"Nothing rises above the baseline between "
                                    f"those scans on {path.name} — try a wider "
                                    f"or better-placed span.")
        else:
            self._gen_active_path = path
            self.status_var.set(
                f"{path.parent.name}/{path.name} · well {rec['well']}: "
                f"{rec['base']} mid {rec['midpoint']} (peak {rec['peak_scan']}) "
                f"height {rec['height_V']:.3f} V, area "
                f"{rec['area_Vscan']:.1f} V·scan")
        self.redraw()
        self._sync_pick_table()

    def _make_gen_cursor(self, ax):
        """Readout + crosshair for one stacked subplot: shows the scan under the
        pointer, its voltage, and the nearest real peak (picked or not) so the
        user can tell which hump they are aiming at before clicking."""
        ann = ax.annotate(
            "", xy=(0.996, 0.985), xycoords="axes fraction", ha="right",
            va="top", fontsize=7, color="#111", zorder=25, visible=False,
            bbox=dict(boxstyle="round,pad=0.3", fc="#FFF9C4", ec="#999",
                      lw=0.6, alpha=0.9))
        vline = ax.axvline(0, color="#B00020", lw=0.7, ls="--", alpha=0.65,
                           zorder=24, visible=False)
        hline = ax.axhline(0, color="#B00020", lw=0.7, ls=":", alpha=0.55,
                           zorder=24, visible=False)
        return (ax, ann, vline, hline)

    def _hide_gen_cursor(self):
        for (_ax, ann, vline, hline) in self._gen_cursors:
            ann.set_visible(False)
            vline.set_visible(False)
            hline.set_visible(False)

    def _gen_nearest_peak(self, pk, scan, vol):
        """The dominant local max within one base on any visible channel, as a
        short label.  Ranked by height (what the user is aiming at), not by
        scan distance, so a flat-channel wobble cannot win over the real peak.
        """
        if vol is None or scan is None:
            return None
        best = None
        rad = max(6, int(round(pk._spacing())))
        col2chan = self._col_to_chan(pk.base_order)
        for col in sorted(pk.col_color):
            try:
                if not self.chan_show[col2chan.get(col, col)].get():
                    continue
            except Exception:
                pass
            for (s, y) in pk._window_peaks(col, scan - rad, scan + rad):
                key = (-y, abs(s - scan))
                if best is None or key < best[0]:
                    best = (key, col, s, y)
        if best is None:
            return None
        (_k, col, s, y) = best
        picked = any(r["col"] == col and abs(r["scan"] - s) <= 2
                     for r in pk.records)
        return (f"peak  Ch{col2chan.get(col, col) + 1} "
                f"{CHANNEL_ORDER[col]}  scan {s} · {y:.0f} V"
                + ("  [picked]" if picked else ""))

    def _on_gen_motion(self, event):
        if not self.genotyping_active:
            return
        if self.area_mode:
            # area mode keeps its own preview artists, so the crosshair guard
            # below must not apply to it
            if self._area_cursors:
                self._area_preview(event)
            return
        if not self._gen_cursors:
            return
        hit = self._gen_axes_hit(event)
        if hit is None:
            if any(ann.get_visible() for (_a, ann, _v, _h) in self._gen_cursors):
                self._hide_gen_cursor()
                self.canvas.draw_idle()
            return
        i, ax, xs, vs = hit
        paths = self._gen_paths()
        if not (0 <= i < len(paths)):
            return
        try:
            pk = self._ensure_picker(paths[i])
        except Exception:
            self._hide_gen_cursor()
            return
        scan = int(round(xs))
        vol = float(vs)
        text = f"cursor  scan {scan} · {vol:.0f} V"
        near = self._gen_nearest_peak(pk, scan, vol)
        if near:
            text += "\n" + near
        for (cax, ann, vline, hline) in self._gen_cursors:
            if cax is ax:
                ann.set_text(text)
                ann.set_visible(True)
                vline.set_xdata([scan, scan])
                vline.set_visible(True)
                hline.set_ydata([vol, vol])
                hline.set_visible(True)
            else:
                ann.set_visible(False)
                vline.set_visible(False)
                hline.set_visible(False)
        self.canvas.draw_idle()

    def _redraw_genotyping(self):
        """Stack the current batch and draw each well's raw channels plus the
        picked-peak overlays, so zoom / µA / Reset view keep working while the
        user picks peaks on the main window"""
        self.fig.clear()
        self._plot_axes = []
        self._axes_keys = []
        self._axes_titles = []
        self._gen_cursors = []
        self._area_cursors = []
        self._full_xlim = None
        self._full_ylim = None
        self._wrap_active = False
        paths = self._gen_paths()
        self._update_page_label(paths)
        if not paths:
            self.canvas.draw_idle()
            self._sync_pick_table()
            return
        gx0 = gx1 = gy0 = gy1 = None
        for i, path in enumerate(paths):
            ax = self.fig.add_subplot(len(paths), 1, i + 1)
            try:
                pk = (self._ensure_area_picker(path) if self.area_mode
                      else self._ensure_picker(path))
                doc = pk.doc
            except Exception as e:
                ax.text(0.5, 0.5, f"Load error: {e}",
                        ha="center", transform=ax.transAxes)
                continue
            n = doc.acgt.shape[0]
            x = np.arange(n)
            plotted = False
            # draw by PHYSICAL channel (run dye order) so pick mode shows the
            # same legend, colours, stacking and show/hide as the viewer
            for ci, base, col, color in self._channel_lines(
                    self._settings_for(doc, self._settings_from_ui())):
                if not self.chan_show[ci].get():
                    continue
                seg = doc.acgt[:, col]
                ax.plot(x, seg, color=color, lw=0.7,
                        label=f"Ch{ci + 1} {base}")
                m0 = float(np.nanmin(seg))
                m1 = float(np.nanmax(seg))
                gy0 = m0 if gy0 is None else min(gy0, m0)
                gy1 = m1 if gy1 is None else max(gy1, m1)
                plotted = True
            if n:
                gx0 = 0 if gx0 is None else 0
                gx1 = n if gx1 is None else max(gx1, n)
            self._plot_axes.append(ax)
            self._axes_keys.append(self._akey(path))
            if not plotted:
                ax.text(0.5, 0.5, "(all channels hidden)", ha="center",
                        va="center", transform=ax.transAxes, color="#888")
            if self.show_current.get() and doc.current_ua is not None:
                axc = ax.twinx()
                cu = doc.current_ua
                axc.plot(np.arange(n), cu, color=CURRENT_COLOR, lw=0.8,
                         alpha=0.8, linestyle=":", label="I (µA)")
                axc.set_ylabel("µA", color=CURRENT_COLOR, fontsize=6)
                axc.tick_params(axis="y", labelcolor=CURRENT_COLOR, labelsize=6)
            if self.area_mode:
                ap = self._area_pickers.get(str(path.resolve()))
                if ap is not None:
                    ap.plot_overlay(ax)
                self._area_cursors.append(self._make_area_markers(ax))
            else:
                pk.plot_overlay(ax)
            self._draw_ctce_marks(ax, path)          # CTCE species, hump, block baselines
            ax.tick_params(labelsize=7)
            base_title = self._graph_title(path, doc.path.parent.name)
            self._axes_titles.append((ax.text(
                0.004, 0.995, base_title, transform=ax.transAxes, ha="left",
                va="top", fontsize=8, color="#222", zorder=6), base_title))
            if not self.area_mode:
                self._gen_cursors.append(self._make_gen_cursor(ax))
            self._style_x_axis(ax, i == len(paths) - 1)
        if gx1 is not None and gy0 is not None and gx1 > gx0:
            pad = 0.02 * (gy1 - gy0) or 1.0
            self._full_xlim = (float(gx0), float(gx1))
            self._full_ylim = (float(gy0 - pad), float(gy1 + pad))
        self._data_full = self._full_xlim
        self._update_ref_full()
        self._refresh_align_titles()
        self.fig.tight_layout(rect=(0.045, 0, 1, 1), h_pad=0.25, pad=0.3)
        self._add_volt_label()
        self._apply_zoom()
        self.canvas.draw_idle()
        self._sync_pick_table()

    def _sync_pick_table(self):
        if getattr(self, "pick_tree", None) is None:
            return
        self.pick_tree.delete(*self.pick_tree.get_children())
        if self.area_mode:
            for path in self._gen_paths():
                pk = self._area_pickers.get(str(path.resolve()))
                if pk is None:
                    continue
                for r in sorted(pk.records,
                                key=lambda r: (int(r["start_scan"]),
                                               int(r["channel"]))):
                    self.pick_tree.insert("", tk.END, values=(
                        r["run"] or "-", r["well"], r["base"],
                        r["start_scan"], r["end_scan"], r["midpoint"],
                        f"{r['height_V']:.3f}", f"{r['area_Vscan']:.1f}"))
            return
        for path in self._gen_paths():
            pk = self._gen_pickers.get(str(path.resolve()))
            if pk is None:
                continue
            for r in sorted(pk.records,
                            key=lambda r: (int(r["scan"]), int(r["col"]),
                                           int(r.get("gid", 0)))):
                fr = pk.clust_frac(r) if r["kind"] == "main" else None
                ma = pk.mass_action(r) if r["kind"] == "main" else None
                self.pick_tree.insert("", tk.END, values=(
                    pk.doc.well, r["scan"], pk.duplex_of(r), r["base"],
                    r["kind"], f"{r['height']:.3f}", f"{r['area']:.1f}",
                    f"{fr:.3f}" if fr else "",
                    f"{ma['mf']:.3f}" if ma else ""))

    def _gen_picker_active(self):
        """The picker to act on: the well of the last click, else the first
        displayed well that already has picks, else the first displayed well."""
        if self._gen_active_path is not None:
            pk = self._gen_pickers.get(str(self._gen_active_path.resolve()))
            if pk is not None:
                return pk
        for path in self._gen_paths():
            pk = self._gen_pickers.get(str(path.resolve()))
            if pk is not None and pk.records:
                return pk
        for path in self._gen_paths():
            pk = self._gen_pickers.get(str(path.resolve()))
            if pk is not None:
                return pk
        return None

    def _sync_gen_opts(self):
        for pk in self._gen_pickers.values():
            pk.include_sh = bool(self._gen_sh.get())
            pk.show_d2 = bool(self._gen_d2.get())
        self.redraw()

    def _gen_undo(self):
        if not self._require_picking("Undo last pick"):
            return
        if self.area_mode:
            rec = None
            if self._gen_active_path is not None:
                pk = self._area_pickers.get(str(self._gen_active_path.resolve()))
                if pk is not None:
                    rec = pk.undo_last()
            if rec is None:
                for path in self._gen_paths():
                    pk = self._area_pickers.get(str(path.resolve()))
                    if pk is not None and pk.records:
                        rec = pk.undo_last()
                        self._gen_active_path = path
                        break
            if rec is None:
                self.status_var.set("Nothing to undo.")
                return
            self.status_var.set(f"Removed {rec['base']} at mid "
                                f"{rec['midpoint']}.")
            self.redraw()
            self._sync_pick_table()
            return
        pk = self._gen_picker_active()
        if pk is None or not pk.undo_last():
            self.status_var.set("Nothing to undo.")
            return
        self.status_var.set("Removed last picked peak.")
        self.redraw()
        self._sync_pick_table()

    def _gen_clear(self):
        if not self._require_picking("Clear picks"):
            return
        store = self._area_pickers if self.area_mode else self._gen_pickers
        if store:
            store.clear()
            self._gen_active_path = None
            self.status_var.set("Measurements cleared." if self.area_mode
                                else "Picks cleared.")
        self.redraw()
        self._sync_pick_table()

    def _gen_clear_shown(self):
        """Remove the picks / measurements of the wells on screen only, so a
        well can be picked or dragged again from scratch (the menu's Clear
        removes them everywhere)."""
        if not self._require_picking("Clear shown wells"):
            return
        store = self._area_pickers if self.area_mode else self._gen_pickers
        n = 0
        for path in self._gen_paths():
            pk = store.pop(str(Path(path).resolve()), None)
            n += len(getattr(pk, "records", []) or []) if pk is not None else 0
        self._gen_active_path = None
        what = "measurement(s)" if self.area_mode else "pick(s)"
        self.status_var.set(f"Cleared {n} {what} in the shown wells.")
        self.redraw()
        self._sync_pick_table()

    def _gen_mark_std_batch(self):
        """Place the IS in every displayed well from the plate's learned shape.

        96 wells is 384 clicks; nobody does that twice.  This uses the geometry
        fitted by "Semi automatic standard…", so the wells are placed the way
        the operator's own marks say rather than by whatever spacing detection
        happens to guess on each trace -- the guess is wrong often enough
        (it locks onto noise in a third of the ABCC2 wells) that leaning on it
        for all 96 is not the same thing as semi-automatic.

        Picked-but-unmarked wells are marked straight away; only a displayed
        well with nothing picked and no shape learned refuses, and it says so
        rather than guessing.  Every mark this places is tagged ``std_auto``,
        so it shows up in the table as detection rather than as something the
        operator picked.
        """
        import genotyping
        if not self._require_picking("Mark IS in all shown wells"):
            return
        paths = self._gen_paths()
        if not paths:
            self.status_var.set("Select the wells to mark first.")
            return

        # The operator's own clicks are the authority: any displayed well that
        # carries picked mains but no standard yet is marked right now, so
        # "pick four peaks on five wells, then Mark IS in all shown wells"
        # works as written instead of refusing until a model exists.
        picked = []
        for path in paths:
            try:
                pk = self._ensure_picker(path)
            except Exception:
                continue
            mains = [r for r in pk.records if r["kind"] == "main"]
            if not mains:
                continue
            if pk.std:
                # Adding a forgotten IS peak to an already-marked well has to
                # work.  If the operator has picked a main the mark does not
                # cover, re-mark from the picks so the new peak joins the
                # standard; an untouched mark (its mains are exactly the ones
                # already in ``std``) is left alone.
                std_scans = [int(x) for x, _ in pk.std]
                if not any(min(abs(m["scan"] - s) for s in std_scans) > 2
                           for m in mains):
                    continue
            picked.append(pk)
        if picked:
            marked, failed, pairs = [], [], []
            for pk in picked:
                try:
                    pk.mark_std()
                    marked.append(pk.doc.well)
                    pairs.append((Path(pk.path), pk))
                except ValueError as e:
                    failed.append((pk.doc.well, str(e)))
            self._library_store(*pairs)
            self._build_pick_table()
            self.redraw()
            msg = (f"Marked {len(marked)} of {len(picked)} picked well(s) "
                   "as the internal standard.")
            if failed:
                msg += ("\nCould not mark: "
                        + ", ".join(f"{w} ({why})" for w, why in failed[:12]))
            msg += ("\n\nThese are saved in the library. Run Auto genotyping ▸ "
                    "Auto-genotype sample window…\nto genotype the plate — "
                    "it learns the standard from these marks itself.")
            messagebox.showinfo("Mark IS in all shown wells", msg)
            return

        channel = simpledialog.askinteger(
            "Standard channel",
            "Which physical channel carries the internal standard?\n\n"
            f"(dye order {self._display_order()}: "
            + ", ".join(f"Ch{c + 1}={b}" for c, b in
                        enumerate(self._display_order()[:4])) + ")",
            initialvalue=int(self._auto_is_channel.get()),
            minvalue=1, maxvalue=4, parent=self)
        if channel is None:
            return
        model = self._model_for(paths[0])
        if model is None:
            messagebox.showinfo(
                "Mark IS in all shown wells",
                "Mark a few wells by hand first (or let Auto-genotype learn "
                "from marks\nit already has), then use this. Nothing is "
                "placed without a shape to\nplace it with, because the "
                "spacing is what varies between fragments.")
            return
        marked, failed, skipped, pairs = [], [], [], []
        self._gen_std_batch_progress = (0, len(paths))
        for path in paths:
            try:
                pk = self._ensure_picker(path)
            except Exception as e:
                failed.append((path.stem, str(e)))
                continue
            if pk.std:
                skipped.append(path.stem)
                continue
            try:
                col = genotyping.acgt_index_for_channel(
                    self._doc_order(pk.doc), channel)
                self._model_mark(pk, self._model_for(path), col)
                marked.append(path.stem)
                pairs.append((path, pk))
            except ValueError as e:
                failed.append((path.stem, str(e)))
            self._gen_std_batch_progress = (len(marked) + len(failed),
                                            len(paths))
            self.update_idletasks()
        self._library_store(*pairs)
        self._build_pick_table()
        self.redraw()
        msg = (f"Marked {len(marked)} of {len(paths)} wells from this plate's "
               f"learned standard.")
        if skipped:
            msg += f"  {len(skipped)} already had an IS set (left alone)."
        if failed:
            msg += f"  {len(failed)} had no quartet found — click those by hand."
            messagebox.showwarning(
                "Mark IS in all shown wells",
                msg + "\n\nNo standard found in:\n"
                + ", ".join(w for w, _ in failed[:24])
                + ("\n…" if len(failed) > 24 else "")
                + "\n\nTurn on Manual peak picking… and click the four IS peaks "
                  "on the standard channel, then Mark peaks as standard…")
        else:
            messagebox.showinfo("Mark IS in all shown wells", msg)

    def report_marked_wells(self):
        """List the wells carrying an IS mark, and say what the model knows.

        Marks live in the pick tables (and the mark library on disk), so this
        is the one place to count them before running Auto over the plate: Auto will learn the
        standard from these marks when none has been learned yet, so what is
        marked here is what the whole plate is fitted to.
        """
        by_run = {}
        for pk in self._gen_pickers.values():
            if not pk.std:
                continue
            run = Path(pk.path).parent.name
            how = pk.std_auto or "hand-marked"
            by_run.setdefault(run, []).append(
                f"{pk.doc.well} ({len(pk.std)} peak"
                f"{'s' if len(pk.std) != 1 else ''}, {how})")
        total = sum(len(v) for v in by_run.values())
        self.status_var.set(f"{total} well(s) marked as the internal standard.")
        if not total:
            messagebox.showinfo(
                "Marked IS wells",
                "No wells are marked as the internal standard in this "
                "session.\n\nTurn on Manual genotyping ▸ Manual peak picking…, click "
                "the four IS peaks\nin a well, then Manual genotyping ▸ Mark peaks as "
                "standard…\nDo that in at least two wells, then Auto-genotype "
                "learns the\nstandard from them and finds it in the rest.")
            return
        lines = []
        for run in sorted(by_run):
            marks = by_run[run]
            lines.append(f"\n{run} — {len(marks)} marked:")
            lines.append("  " + ", ".join(marks))
        model = (f"Learned: {self._auto_is_model_note}"
                 if self._auto_is_model is not None else
                 "Not learned yet — Auto-genotype will fit it from these marks.")
        messagebox.showinfo(
            "Marked IS wells",
            f"{total} well{'' if total == 1 else 's'} marked as the internal "
            f"standard in this session:\n" + "".join(lines) +
            f"\n\nPlate standard: {model}\n")

    def library_show(self):
        """What the library on disk holds: which runs are marked, how much of
        each has a call."""
        try:
            rows = mark_library.summary()
        except Exception:
            rows = []
        head = f"Library: {mark_library.LIBRARY_DIR}"
        if not rows:
            messagebox.showinfo(
                "Standard library",
                head + "\n\n(empty)\n\nMarked internal standards and "
                "auto-genotype calls are saved\nthere as you make them, so a "
                "restart never costs a marked plate.")
            return
        lines = [head, ""]
        for run, marked, results, updated in rows[:20]:
            bits = f"{marked} marked"
            if results:
                bits += f", {results} with a call"
            lines.append(f"  {run}\n      {bits}   {str(updated)[:16]}")
        if len(rows) > 20:
            lines.append(f"  … and {len(rows) - 20} more runs")
        messagebox.showinfo("Standard library", "\n".join(lines))

    def library_export_ml(self):
        """Write the library's labelled rows out as training data.

        One row per genotyped well: what was measured (peak heights, SNRs,
        the standard's scans) plus the call and whose call it is, so the
        operator's corrections can be used as labels and the engine's own
        calls filtered out.
        """
        wells = {}
        try:
            for data in mark_library.runs():
                wells[str(data.get("run", "?"))] = data.get("wells") or {}
        except Exception:
            pass
        for run, entries in getattr(self, "_library_cache", {}).items():
            wells.setdefault(run, {}).update(entries)
        rows = mark_library.ml_rows(wells)
        if not rows:
            messagebox.showinfo(
                "Export library (ML)",
                "Nothing to export yet.\n\nRun Auto-genotype over a plate "
                "first — every call it makes,\nalong with your corrections to "
                "it, becomes one training row.")
            return
        try:
            mark_library.LIBRARY_DIR.mkdir(parents=True, exist_ok=True)
        except OSError:
            pass
        try:
            path = filedialog.asksaveasfilename(
                title="Export library (ML)", defaultextension=".csv",
                initialdir=str(mark_library.LIBRARY_DIR),
                initialfile="library.csv",
                filetypes=[("CSV", "*.csv"), ("JSON Lines", "*.jsonl")],
                parent=self)
        except tk.TclError:
            path = ""
        if not path:
            return
        try:
            n = mark_library.write_ml(path, rows)
        except OSError as e:
            messagebox.showerror("Export library (ML)", str(e))
            return
        self.status_var.set(f"Exported {n} library rows to {path}")
        messagebox.showinfo("Export library (ML)", f"Wrote {n} rows to\n{path}")

    def _gen_lock_col(self):
        """The locked acgt column for picking, or None to let the picker guess."""
        v = (self._gen_lock_var.get() or "any").strip()
        if not v.isdigit():
            return None
        return self._col_to_chan().get(int(v) - 1)

    def _require_picking(self, what):
        """Guard for the manual picking commands.

        These all used to ``return`` silently when picking was off, so the menu
        item looked dead with no explanation -- worst right after a restart,
        where picking starts off and every manual command silently no-ops.
        """
        if self.genotyping_active:
            return True
        messagebox.showinfo(
            what,
            "Turn on Manual genotyping ▸ Manual peak picking… first.\n\n"
            "Then click each of the four IS peaks on the standard channel,\n"
            "then Manual genotyping ▸ Mark peaks as standard…")
        return False

    def _gen_mark_std(self):
        """Mark the standard in every well the operator picked peaks in.

        Not just the last well clicked: picking four peaks on four wells and
        pressing this once has to mark four wells, otherwise the operator
        gets one mark and three silent wells and has no way to tell that
        anything was left out.  A well that already carries a mark is left
        alone unless it is the one being worked on, so re-marking is still
        possible without disturbing the rest of the plate.
        """
        if not self._require_picking("Mark peaks as standard"):
            return
        active = self._gen_picker_active()
        targets = []
        for path in self._gen_paths():
            pk = self._gen_pickers.get(str(path.resolve()))
            if pk is None:
                continue
            mains = [r for r in pk.records if r["kind"] == "main"]
            if not mains:
                continue
            if pk.std and pk is not active:
                # An operator's mark is never overwritten -- unless they have
                # picked a main it does not cover, which is how a forgotten IS
                # peak gets added back to an already-marked well.
                std_scans = [int(x) for x, _ in pk.std]
                if not any(min(abs(m["scan"] - s) for s in std_scans) > 2
                           for m in mains):
                    continue
            targets.append(pk)
        if not targets:
            self.status_var.set("Select wells first.")
            return
        length = simpledialog.askstring(
            "Internal standard",
            "Fragment length (bp) of the standard — the same for the four\n"
            "CTC-CE duplex peaks (HOM1/HOM2/HET1/HET2).\n"
            "Leave empty if unknown.",
            parent=self)
        if length is None:
            return
        length_bp = None
        if length.strip():
            try:
                length_bp = float(length.strip())
            except ValueError:
                messagebox.showerror("Internal standard",
                                     f"'{length}' is not a number.")
                return
        done, failed, pairs = [], [], []
        for pk in targets:
            try:
                pk.mark_std(length_bp)
            except ValueError as e:
                failed.append(f"{pk.doc.well}: {e}")
                continue
            done.append(pk.doc.well)
            self._gen_active_path = Path(pk.path)
            pairs.append((Path(pk.path), pk))
        self._library_store(*pairs)
        self._sync_pick_table()
        if len(done) == 1 and not failed:
            self.status_var.set(f"Well {done[0]}: Standard set.")
        elif done:
            self.status_var.set(
                f"Standard set in {len(done)} well(s): {', '.join(done)}.")
        else:
            self.status_var.set("No standard set.")
        if failed:
            messagebox.showwarning("Internal standard",
                                   "\n".join(failed[:12]))

    def _gen_mark_duplex(self):
        """Tag the active well's picked mains as the duplex species of one
        allelic position, so it reports the CTCE mass-action MF."""
        if not self._require_picking("Tag duplex species for MF"):
            return
        pk = self._gen_picker_active()
        if pk is None:
            self.status_var.set("Select wells first.")
            return
        try:
            msg = pk.mark_duplex()
        except ValueError as e:
            messagebox.showwarning("Duplex species", str(e))
            return
        self._gen_active_path = Path(pk.path)
        ma = next((pk.mass_action(m) for m in pk.labelled_species()), None)
        if ma:
            msg += f"   MF {ma['mf']:.3f}"
        self.status_var.set(f"Well {pk.doc.well}: {msg}")
        self._sync_pick_table()
        self.redraw()

    # ------------------------------------------------------- auto-genotyping
    def _auto_targets(self):
        """The wells to score: the selection, or the whole sample window.

        Selecting files one at a time to genotype a plate was pure friction --
        the sample list already shows exactly the run being worked on, so with
        nothing selected the whole window is taken instead.
        """
        picked = [Path(p) for p in (self.selected or []) if Path(p).is_file()]
        if picked:
            return picked
        return [Path(p) for p in getattr(self, "files", [])
                if Path(p).is_file()]

    def _model_for(self, path):
        """The learned internal standard of the run this file belongs to."""
        return (self._auto_models.get(Path(path).parent.name)
                or self._auto_is_model)

    def auto_genotype_channels(self):
        """Ask which physical channels carry the standard and the sample.

        The engine needs to know which trace holds the internal-standard
        quartet and which holds the sample's four duplexes.  That is a property
        of the assay on the plate, not of the dye order, so it is asked for
        explicitly rather than assumed.  The pairing is checked here, before
        any well is scored: if the two land on the same acgt column the engine
        would compare a channel against itself and every well would come back
        "no standard quartet" for no visible reason.
        """
        win = tk.Toplevel(self)
        win.title("Auto-genotyping — channel roles")
        win.transient(self)
        win.resizable(False, False)
        body = ttk.Frame(win, padding=12)
        body.pack(fill="both", expand=True)
        ttk.Label(
            body,
            text="Which physical channel carries what on this run?\n"
                 "The dye order only says which base sits on Ch1–Ch4; it "
                 "does not\nsay which channel holds the standard.",
            justify="left").grid(row=0, column=0, columnspan=2, sticky="w",
                                 pady=(0, 10))
        for row, (var, label) in enumerate(
                ((self._auto_is_channel, "Internal-standard channel"),
                 (self._auto_sample_channel, "Sample channel")), start=1):
            text = f"{label} ({self._channel_label(var.get())}):"
            ttk.Label(body, text=text).grid(row=row, column=0, sticky="w",
                                             pady=3)
            box = ttk.Combobox(body, textvariable=var, state="readonly",
                               values=[1, 2, 3, 4], width=6)
            box.grid(row=row, column=1, sticky="w", padx=(8, 0), pady=3)
            var.trace_add(
                "write",
                lambda *_a, var=var, row=row: self._auto_label_text(body, row, var))
        self._auto_hint = ttk.Label(body, text="", foreground="#8A2A0A")
        self._auto_hint.grid(row=3, column=0, columnspan=2, sticky="w",
                             pady=(8, 0))
        btns = ttk.Frame(body)
        btns.grid(row=4, column=0, columnspan=2, sticky="e", pady=(12, 0))
        ttk.Button(btns, text="Cancel",
                   command=win.destroy).pack(side="right", padx=(6, 0))
        ttk.Button(btns, text="Save", command=self._auto_save_channels).pack(
            side="right")

        def check(*_a):
            is_c = int(self._auto_is_channel.get())
            sa_c = int(self._auto_sample_channel.get())
            if is_c == sa_c:
                self._auto_hint.config(
                    text=f"Pick two different channels — Ch{is_c} cannot be "
                         "both the standard and the sample.")
            else:
                self._auto_hint.config(text="")

        self._auto_is_channel.trace_add("write", check)
        self._auto_sample_channel.trace_add("write", check)
        check()

    def auto_genotype_min_frac(self):
        """Ask how big the minor allele has to be before a het counts.

        A bump beside a strong peak is called a het by any threshold that has
        never seen this operator's instrument, and the caller had no way to
        say "not on my plate" -- the ``ai`` flag only said the balance was
        odd.  So the floor is asked for: 0 keeps the engine's own behaviour,
        anything above it stops a small second band being read as a second
        allele.  It is remembered with the marks, because a plate scored under
        one floor should not quietly be rescored under another.
        """
        current = float(self._auto_min_het_frac.get())
        value = simpledialog.askfloat(
            "Heterozygote floor",
            "Smallest minor-allele fraction still called a het:\n\n"
            "0 keeps the engine's own behaviour -- every het is called and\n"
            "an off-balance one is flagged \"ai\".  0.10, say, asks that the\n"
            "smaller allele carry at least 10% of the position; below that\n"
            "the well falls through to the homozygote tests and is flagged\n"
            "\"het-floor\" so the table says why it is not a het.\n\n"
            "Applies from the next Auto-genotype on.",
            parent=self, initialvalue=current, minvalue=0.0, maxvalue=0.5)
        if value is None:                                 # cancelled
            return
        value = float(value)
        self._auto_min_het_frac.set(value)
        prefs = mark_library.load_prefs()
        prefs["min_het_frac"] = value
        try:
            mark_library.save_prefs(prefs)
        except Exception as exc:                          # noqa: BLE001
            print(f"library: could not save prefs: {exc}", file=sys.stderr)
        self.status_var.set(
            "Het floor: " + ("engine default (0)" if value <= 0
                             else f"minor allele {value:.2f}")
            + " — applies from the next Auto-genotype")

    def _auto_label_text(self, body, row, var):
        """Keep each channel row's label showing the base that channel holds."""
        labels = ("Internal-standard channel", "Sample channel")
        name = labels[row - 1]
        for w in body.winfo_children():
            if isinstance(w, ttk.Label) and w.cget("text").startswith(name):
                w.config(text=f"{name} ({self._channel_label(var.get())}):")
                break

    def _channel_label(self, channel):
        """'Ch3 (T)' for the given physical channel, per the run's dye order."""
        from genotyping import acgt_index_for_channel
        try:
            base = "ACGT"[acgt_index_for_channel(
                self._display_order(), channel)]
        except Exception:
            base = "?"
        return f"Ch{channel} ({base})"

    def _auto_save_channels(self):
        is_c = int(self._auto_is_channel.get())
        sa_c = int(self._auto_sample_channel.get())
        if is_c == sa_c:
            messagebox.showerror(
                "Auto-genotyping — channel roles",
                f"Ch{is_c} is set as both the internal standard and the "
                "sample.\n\nPick two different channels: the engine locates "
                "the standard quartet\non one channel and measures the "
                "sample's duplexes on another.")
            return
        for w in self.winfo_children():
            if isinstance(w, tk.Toplevel) and w.title().startswith(
                    "Auto-genotyping"):
                w.destroy()
        self.status_var.set(
            f"Auto-genotyping: internal standard Ch{is_c}, sample Ch{sa_c} "
            f"({self.base_order_var.get()}).")

    def semi_auto_is(self):
        """Learn the plate's internal standard from a few marked wells, then
        place it in the rest.

        This is the semi-automatic mode.  There is deliberately no automatic
        fallback to any built-in geometry: the internal standard's spacing is a
        property of the fragment and the run, and one assay's numbers are wrong
        for another.  So the operator marks it where they can see it, and this
        propagates that knowledge to the wells they did not mark -- reporting
        how each well was matched, so the ones found only by position can be
        spotted and corrected rather than trusted blindly.

        Every well the model places is written back into the pick table as an
        ordinary IS mark, so it can be reviewed, moved or cleared like any
        hand mark.  Nothing here overwrites a mark the operator made.
        """
        import genotyping
        targets = self._auto_targets()
        if not targets:
            messagebox.showinfo(
                "Semi automatic internal standard",
                "Select the wells of one run in the list first.\n\nThe "
                "standard is learned per plate, so pick one run's worth of "
                "wells (or a subset of it).")
            return
        is_c = int(self._auto_is_channel.get())
        sa_c = int(self._auto_sample_channel.get())
        try:
            run_order = self._doc_order(self._ensure_doc(Path(targets[0])))
            is_col = genotyping.acgt_index_for_channel(run_order, is_c)
            samp_col = genotyping.acgt_index_for_channel(run_order, sa_c)
        except ValueError as e:
            messagebox.showerror("Semi automatic internal standard", str(e))
            return

        runs = {p.parent.name for p in targets}
        if len(runs) > 1:
            messagebox.showerror(
                "Semi automatic internal standard",
                f"The selection spans {len(runs)} runs:\n"
                f"{', '.join(sorted(runs)[:4])}"
                f"{'…' if len(runs) > 4 else ''}\n\n"
                "Each run has its own capillary drift and needs its own marks, "
                "so\nlearn one run at a time.")
            return

        # Seeds = wells the operator already marked, by hand or by IS.csv.
        seeds, traces, load_errors = {}, {}, []
        for path in targets:
            try:
                doc = self._ensure_doc(path)
            except Exception as e:
                load_errors.append(f"{path.stem}: {e}")
                continue
            traces[path.stem] = np.asarray(doc.acgt, dtype=float)
            try:
                pk = self._ensure_picker(path)
            except Exception:
                continue
            if pk.std:
                seeds[path.stem] = [int(pair[0]) for pair in pk.std]
        if load_errors:
            self.status_var.set(
                f"{len(load_errors)} of {len(targets)} wells could not be "
                f"read; skipped.")
        if len(seeds) < 2:
            if seeds:
                where = "Only " + ", ".join(sorted(seeds))
            else:
                where = "No marked wells found"
            messagebox.showinfo(
                "Semi automatic internal standard",
                f"{where}.\n\n"
                "Mark the internal standard in at least two wells first.\n\n"
                "On the standard channel, pick the standard's main peaks and "
                "hit\nManual genotyping ▸ Mark peaks as standard… for each well. "
                "Two is\nthe minimum and three to five is better: the more "
                "you mark, the\ntighter the learned shape and the fewer wells "
                "need checking.")
            self.status_var.set(
                f"Semi automatic: {len(seeds)} marked well"
                f"{'s' if len(seeds) != 1 else ''}, need at least 2.")
            return

        self.status_var.set(
            f"Learning the standard on {Path(targets[0]).parent.name} from "
            f"{len(seeds)} marked well{'' if len(seeds) == 1 else 's'}…")
        self.update_idletasks()
        try:
            model, found = genotyping.semi_auto_plate_model(
                traces, seeds, is_col, samp_col=samp_col)
        except ValueError as e:
            messagebox.showerror("Semi automatic internal standard", str(e))
            return

        # Write the learned standards back into the pickers, as ordinary marks.
        names4 = ["HOM1", "HOM2", "HET1", "HET2"]
        names3 = ["HOM1", "HOM2", "HET"]
        added, pairs = 0, []
        for well, (scans, how, label, _score) in found.items():
            path = next((p for p in targets if p.stem == well), None)
            if path is None:
                continue
            pk = self._ensure_picker(path)
            if pk.std:
                continue            # an operator mark is never overwritten
            names = names4 if len(scans) >= 4 else names3
            pk.std = [(x, names[i]) for i, x in enumerate(scans[:len(names)])]
            pk.length_bp = None
            pk.std_auto = how
            for x in scans[:len(names)]:
                pk.add_std_mark(x, is_col)
            pairs.append((path, pk))
            added += 1

        run = Path(targets[0]).parent.name
        self._library_store(*pairs)
        self._remember_model(run, model)
        self._auto_is_model = model
        self._auto_is_model_note = model.describe_plate()
        self._build_pick_table()
        self.redraw()

        counts = {}
        for _s, how, label, _sc in found.values():
            counts[label] = counts.get(label, 0) + 1
        summary = ", ".join(f"{v} {k}" for k, v in sorted(counts.items()))
        msg = (f"Learned the standard on {Path(targets[0]).parent.name}:\n\n"
               f"  {model.describe_plate()}\n\n"
               f"Placed it in {len(found)} of {len(traces)} wells "
               f"({summary}).\n"
               f"{added} of them are now marked in the pick table.\n\n")
        if counts.get("positional"):
            msg += (f"{counts['positional']} well(s) matched on position "
                    "only --\nequimolar peaks near where this plate's "
                    "standards always are,\nwithout a shape to confirm them. "
                    "Check those by eye; they carry a\nstd-positional flag in "
                    "the auto-genotype table.\n\n")
        if counts.get("fair") or counts.get("weak"):
            msg += ("A few wells are weak; they carry std-weak and are worth "
                    "a look.\n\n")
        msg += "Now run Auto genotyping ▸ Auto-genotype sample window…."
        messagebox.showinfo("Semi automatic internal standard", msg)

    def _model_mark(self, pk, model, col):
        """Place the plate's learned standard on one picker, as a tagged mark.

        Writes ``std`` and ``std_auto``, so the mark behaves like any other:
        reviewable, movable and clearable in the pick table, but distinguishable
        from something the operator picked.  Raises ValueError with a reason
        when the well has no standard at that position.
        """
        import genotyping
        y = np.asarray(pk.doc.acgt[:, col], dtype=float)
        found = model.find(y, cut=genotyping.DEFAULT_IS_CUT)
        if found is None:
            raise ValueError("no standard at this well's usual position")
        scans, _heights, how = found
        names = (["HOM1", "HOM2", "HET1", "HET2"] if len(scans) >= 4
                 else ["HOM1", "HOM2", "HET"][:len(scans)])
        pk.std = [(x, names[i]) for i, x in enumerate(scans[:len(names)])]
        pk.length_bp = None
        pk.std_auto = how
        for x in scans[:len(names)]:
            pk.add_std_mark(x, col)

    def _remember_model(self, run, model):
        """Register a fitted standard under its run and write it to the library.

        The in-memory copy is what the current session uses; the library copy is
        what lets the fit seed a later plate and survive a restart.  The write
        is best-effort -- a model that cannot be saved costs a re-fit, never a
        plate, so :mod:`mark_library` swallows the failure itself.
        """
        self._auto_models[run] = model
        try:
            mark_library.save_model(run, model.to_dict())
        except Exception as exc:                          # noqa: BLE001
            print(f"library: could not save model {run}: {exc}", file=sys.stderr)

    def forget_plate_model(self):
        """Drop the learned plate model and its marks, back to hand-marking.

        The learned IS marks stay in the pick table until they are cleared
        there; this only forgets the model that propagates to other wells, so
        the operator can teach it again from a different set of seeds.
        """
        self._auto_is_model = None
        self._auto_is_model_note = ""
        self._auto_models = {}
        self._auto_saved_model = None
        self._auto_saved_run = ""
        self.status_var.set(
            "Forgot the learned plate model. The IS marks already written to "
            "the pick table\nremain until you clear them there.")

    def auto_use_saved_model(self):
        """Load a standard fitted on an earlier plate, to seed this one.

        Every fitted standard is written to the library, so a shape learned
        once can be handed to a later plate of the same assay without
        re-marking the wells it took to fit it.  It is used only for a run that
        has no marks of its own to fit from -- a plate the operator has already
        marked always wins -- and is applied by the next Auto-genotype.
        "Forget learned standard" drops it again.
        """
        try:
            stored = mark_library.saved_models()
        except Exception:                                 # noqa: BLE001
            stored = []
        import genotyping
        picks = []
        for run, updated, data in reversed(stored):        # newest first
            try:
                model = genotyping.PlateISModel.from_dict(data)
            except Exception:                              # noqa: BLE001
                model = None
            if model is not None:
                picks.append((run, str(updated)[:16], model))
        if not picks:
            messagebox.showinfo(
                "Use a saved standard",
                "No standard is stored yet.\n\nA plate's standard is saved "
                "the first time it is\nlearned (Semi automatic, or Auto "
                "genotype on a marked\nplate) and can then seed a later "
                "plate here.")
            return
        lines = [f"{i}) {run}   ({updated})\n      {model.describe_plate()}"
                 for i, (run, updated, model) in enumerate(picks, start=1)]
        choice = simpledialog.askinteger(
            "Use a saved standard",
            "Seed this plate with a standard fitted earlier.\nA run with marks "
            "of its own still fits them first.\n\n"
            + "\n".join(lines) + "\n\nNumber:",
            parent=self, minvalue=1, maxvalue=len(picks))
        if not choice:
            return
        run, _updated, model = picks[choice - 1]
        self._auto_saved_model = model
        self._auto_saved_run = run
        self._auto_is_model = model
        self._auto_is_model_note = f"seeded from {run}: {model.describe_plate()}"
        self.status_var.set(
            f"Saved standard from {run} loaded — the next Auto-genotype on a "
            f"run with no marks uses it ({model.describe_plate()}).")

    def _auto_learn_from_marks(self, targets, is_c, sa_c, base):
        """Fit the plate model from the marks on this plate, when there are any.

        Marking a few wells and then hitting Auto is the workflow people
        actually use, so the model is fitted here too rather than only by the
        separate Semi automatic item -- otherwise those marks silently protect
        only the marked wells and the rest of the plate falls back to the
        no-standard path, which is how a plate comes back as a handful of calls
        and ninety-some no-calls.

        The geometry still comes from nowhere but the operator's marks: this
        fits exactly what Semi automatic would have fitted, from the same
        seeds, and writes the placements into the pick table the same way.
        Marks restored from the library count like hand marks -- that is the
        point of saving them.  Nothing is invented when there is nothing (or
        almost nothing) to learn from: it just returns and the plate is
        genotyped as marked.  ``targets`` must be one run's wells; the model
        is fitted per run and registered under it.
        """
        import genotyping
        runs = {p.parent.name for p in targets}
        if len(runs) != 1:
            return
        run = next(iter(runs))
        try:
            is_col = genotyping.acgt_index_for_channel(base, is_c)
            samp_col = genotyping.acgt_index_for_channel(base, sa_c)
        except ValueError:
            return
        seeds, traces = {}, {}
        for path in targets:
            try:
                doc = self._ensure_doc(path)
                pk = self._ensure_picker(path)
            except Exception:
                continue
            traces[path.stem] = np.asarray(doc.acgt, dtype=float)
            if pk.std:
                seeds[path.stem] = [int(pair[0]) for pair in pk.std]
        if len(seeds) < 2:
            return
        self.status_var.set(
            f"Learning the standard from {len(seeds)} marked well"
            f"{'s' if len(seeds) != 1 else ''}…")
        self.update_idletasks()
        try:
            model, found = genotyping.semi_auto_plate_model(
                traces, seeds, is_col, samp_col=samp_col)
        except ValueError:
            return
        added, pairs = 0, []
        for well, (scans, how, label, _score) in found.items():
            path = next((p for p in targets if p.stem == well), None)
            if path is None:
                continue
            pk = self._ensure_picker(path)
            if pk.std:
                continue            # an operator mark is never overwritten
            names = ["HOM1", "HOM2", "HET1", "HET2"] if len(scans) >= 4 \
                else ["HOM1", "HOM2", "HET"]
            pk.std = [(x, names[i]) for i, x in enumerate(scans[:len(names)])]
            pk.length_bp = None
            pk.std_auto = how
            for x in scans[:len(names)]:
                pk.add_std_mark(x, is_col)
            pairs.append((path, pk))
            added += 1
        self._library_store(*pairs)
        self._remember_model(run, model)
        self._auto_is_model = model
        self._auto_is_model_note = model.describe_plate()
        self._build_pick_table()
        self.redraw()
        self.status_var.set(
            f"Learned the standard from {len(seeds)} marked well(s): "
            f"{model.describe_plate()} — found in {len(found)} well(s), "
            f"{added} newly marked.")
        return model

    def auto_genotype_wells(self):
        """Score the sample window with no clicking, and show the results.

        Only the selection is taken when there is one; otherwise every file in
        the sample list, because selecting 96 files one at a time was pure
        friction.  Progress goes to the status bar as it goes and the rows land
        in their own table, separate from the manual pick/area table, so the
        two never overwrite each other.
        """
        import genotyping
        targets = self._auto_targets()
        if not targets:
            messagebox.showinfo(
                "Auto-genotype sample window",
                "The sample window is empty — add a data folder or run first.")
            return
        is_c = int(self._auto_is_channel.get())
        sa_c = int(self._auto_sample_channel.get())
        if is_c == sa_c:
            messagebox.showerror(
                "Auto-genotype sample window",
                f"Ch{is_c} is set as both the internal standard and the "
                "sample.\n\nUse Auto genotyping ▸ Channel "
                "roles to\npoint them at two different channels.")
            return
        # One standard shape per run: fit each run's from its own marks (marks
        # restored from the library count), so a window holding several runs
        # gets all of them genotyped instead of one run's shape applied
        # everywhere.  Marking a few wells and pressing Auto is the way this
        # gets used, so the fitting happens here and not in a separate step.
        groups = {}
        for p in targets:
            groups.setdefault(p.parent.name, []).append(p)
        models = {}
        for run, paths in groups.items():
            # each run in its own dye order (declared by its files, or chosen)
            run_base = self._doc_order(self._ensure_doc(paths[0]))
            models[run] = (self._auto_learn_from_marks(paths, is_c, sa_c, run_base)
                           or self._auto_models.get(run)
                           or self._auto_saved_model)
        self._auto_rows = []
        self.status_var.set(f"Auto-genotyping {len(targets)} wells…")
        self.update_idletasks()
        for i, path in enumerate(targets, start=1):
            try:
                doc = self._ensure_doc(path)
            except Exception as e:            # unreadable file, not fatal
                self._auto_rows.append(
                    {"run": path.parent.name, "well": path.stem,
                     "call": "no-call", "frac": 0.0, "flags": "load-error",
                     "hom1": 0.0, "hom2": 0.0, "het1": 0.0, "het2": 0.0,
                     "snr1": 0.0, "snr2": 0.0, "snr3": 0.0, "snr4": 0.0,
                     "is_channel": is_c, "sample_channel": sa_c,
                     "std_scans": "", "reason": str(e)})
            else:
                # Hand over the operator's own STD picks for this well when they
                # exist.  find_is_quartet accepts any equimolar set, so on a
                # weak-IS plate it can prefer four equal noise ripples to the
                # real standard; the operator's marks are authoritative.
                #
                # Marks the semi-automatic pass wrote are NOT treated as manual.
                # They came out of the plate model, so re-detecting them through
                # it is what produces the source and confidence columns -- and
                # feeding them back as "manual" would report every learned well
                # as an operator mark and hide exactly the wells worth checking.
                manual_std = None
                try:
                    pk = self._ensure_picker(path)
                    if pk.std and not pk.std_auto:
                        # pk.std is [(scan, name), ...]; the engine wants scans.
                        manual_std = [int(pair[0]) for pair in pk.std]
                except Exception:
                    manual_std = None
                self._auto_rows.append(genotyping.auto_genotype(
                    doc, is_channel=is_c, sample_channel=sa_c,
                    base_order=self._doc_order(doc), cut=genotyping.DEFAULT_IS_CUT,
                    run_name=path.parent.name,
                    std_scans_manual=manual_std,
                    is_model=models.get(path.parent.name),
                    min_het_frac=float(self._auto_min_het_frac.get())))
            if i % 8 == 0 or i == len(targets):
                self.status_var.set(f"Auto-genotyping… {i}/{len(targets)}")
                self.update_idletasks()
        # What was measured and called goes to the library as it stands, so
        # the machine-learning export has a labelled row for every well the
        # plate produced -- not only the ones with marks.
        self._library_store_results(self._auto_rows)
        self._build_auto_table()
        calls = {}
        for r in self._auto_rows:
            calls[r["call"]] = calls.get(r["call"], 0) + 1
        summary = ", ".join(f"{v} {k}" for k, v in sorted(calls.items()))
        # Without a learned model, wells with no hand mark cannot be scored at
        # all -- there is no built-in geometry to fall back on, deliberately.
        # Say so up front rather than returning a table of no-calls the user has
        # to work out the meaning of.
        unplaced = sum(1 for r in self._auto_rows
                       if r.get("std_source") in ("none", ""))
        note = ""
        if unplaced and self._auto_is_model is None:
            note = (f"  {unplaced} well(s) had no hand-marked standard: run "
                    "Auto genotyping ▸\n  Semi automatic "
                    "standard first.")
        elif unplaced:
            note = f"  {unplaced} well(s) had no standard placed."
        if self._auto_is_model_note:
            note = f"  {self._auto_is_model_note}" + note
        self.status_var.set(
            f"Auto-genotyped {len(self._auto_rows)} wells "
            f"(standard Ch{is_c}, sample Ch{sa_c}): {summary}.{note}")
        self.redraw()

    def _build_auto_table(self):
        """The genotyping table under the graphs.

        Every row the plate has, with the call, the standard it was read
        against and the peaks the call came from -- the three things an
        operator needs to judge a well, in one place under the traces that
        show them.

        It must not replace the pick table.  Manual STD picks are the ground
        truth for this assay -- auto-genotyping defers to them -- and the pick
        table is what ``_gen_save`` writes, so manual picking takes this slot
        over while it is on and hands it back on exit, rather than the two
        tables fighting for the same space.

        Built once and refilled by ``_sync_geno_table``; the rows come from
        this session's auto-genotype overlaid on the library, so a restart
        still shows the plate that was scored before it.
        """
        if self._auto_tree is not None:
            try:
                self._auto_tree.master.destroy()
            except tk.TclError:
                pass
            self._auto_tree = None
            self._auto_table = None
        cols = ("well", "call", "frac", "hom1", "hom2", "het1", "het2",
                "snr1", "snr2", "snr3", "snr4", "ispeaks", "called", "conf",
                "flags", "reason")
        self._auto_cols = cols              # the editors address columns by
        #                                name, so a new column cannot shift one
        tblf = ttk.Frame(self.below)
        bar = ttk.Frame(tblf)
        bar.pack(fill=tk.X, padx=2, pady=(2, 0))
        ttk.Label(bar, text="double-click a call or the IS peaks to correct "
                            "them").pack(side=tk.LEFT)
        ttk.Button(bar, text="Save table\u2026",
                   command=self.auto_genotype_save).pack(side=tk.RIGHT)
        # The pick table has to exist and be visible for manual picks to stay
        # saveable, so build it if auto-genotyping is the first thing run.
        if self.pick_tree is None and self.genotyping_active and not self.area_mode:
            self._build_pick_table()
        tree = ttk.Treeview(tblf, columns=cols, show="headings", height=8)
        widths = {"well": 56, "call": 74, "frac": 48, "hom1": 56, "hom2": 56,
                  "het1": 56, "het2": 56, "snr1": 46, "snr2": 46, "snr3": 46,
                  "snr4": 46, "ispeaks": 104, "called": 104, "conf": 76,
                  "flags": 90, "reason": 200}
        left = ("well", "call", "ispeaks", "called", "conf", "flags", "reason")
        headings = {"ispeaks": "IS peaks", "called": "called peaks"}
        for c in cols:
            tree.heading(c, text=headings.get(c, c))
            tree.column(c, width=widths[c],
                        anchor="w" if c in left else "e",
                        stretch=(c in ("reason", "flags", "ispeaks", "called")))
        # The operator's own standard marks are the ground truth this assay
        # defers to; tint those rows so they stand apart from learned ones.
        tree.tag_configure("manual", background="#EAF4FF")
        # A flagged well is the one worth a second look.  As plain text in a
        # column of mostly-blank cells it was easy to miss, so the row itself
        # turns red.
        tree.tag_configure("flagged", foreground="#B00020")
        # Double-click opens the editor for call or IS
        tree.bind("<Double-1>", self._auto_edit_call)
        # Selecting a line navigates: the table is a second view of the same
        # wells, so it moves the graphs rather than only highlighting itself.
        tree.bind("<<TreeviewSelect>>", self._geno_row_select)
        # Same arrow paging as the pick table, so the graph keeps moving under
        # the cursor whichever table the operator is reading from.
        for seq in ("<Up>", "<Down>", "<Left>", "<Right>",
                    "<Page_Up>", "<Page_Down>"):
            tree.bind(seq, self._page_key)
        vs = ttk.Scrollbar(tblf, orient=tk.VERTICAL, command=tree.yview)
        tree.configure(yscrollcommand=vs.set)
        tree.pack(side=tk.LEFT, fill=tk.X, expand=True)
        vs.pack(side=tk.RIGHT, fill=tk.Y)
        self._auto_tree = tree
        self._auto_table = tblf
        # Which pane owns the slot underneath the graphs is decided in one
        # place, so building the table -- after a run, or on a restart's first
        # load -- can never leave it sitting there unattached while a manual
        # mode has already claimed the space.
        self._below_refresh()
        self._sync_geno_table()

    def _geno_row_source(self):
        """The table's rows: this session's results over the library's.

        The library half is what makes a restart still show the plate -- the
        scored calls live on disk -- and this session's own rows win where
        they overlap, because they are what the engine just measured.  Only
        the runs in play are read, so an unrelated plate in the library does
        not appear under whichever traces are open.
        """
        rows, index = [], {}
        runs = {r.get("run") for r in self._auto_rows if r.get("run")}
        runs |= {p.parent.name for p in (self.files or [])}
        runs |= {p.parent.name for p in (self.selected or [])}
        for run in sorted(x for x in runs if x):
            wells = self._library_wells(run) or {}
            for well in sorted(wells):
                res = (wells.get(well) or {}).get("result")
                if not res:
                    continue
                row = dict(res)
                row["run"], row["well"] = run, well
                rows.append(row)
                index[(run, well)] = row
        for r in self._auto_rows:
            key = (r.get("run"), r.get("well"))
            old = index.get(key)
            if old is None:
                rows.append(r)
            else:
                for j, existing in enumerate(rows):
                    if existing is old:
                        rows[j] = r
                        break
            index[key] = r
        return rows

    def _sync_geno_table(self):
        """Refill the table -- after a run, an edit, or a restart's load."""
        tree = getattr(self, "_auto_tree", None)
        if tree is None:
            return
        for iid in tree.get_children():
            tree.delete(iid)
        self._geno_rows = self._geno_row_source()
        self._auto_iids = {}
        for i, r in enumerate(self._geno_rows):
            # IS peak column: where the standard was found, so the operator can
            # eyeball it against the traces -- ``std_scans`` is set for both
            # detection and manual marks, ``std_scans_manual`` backs it up.
            ispeaks = r.get("std_scans", "") or r.get("std_scans_manual", "")
            conf = " ".join(x for x in (
                r.get("std_source", ""),
                r.get("std_conf", ""),
                (f"{r.get('std_score', 0.0):.2f}" if r.get("std_score") else ""),
            ) if x)
            iid = str(i)
            tree.insert("", tk.END, iid=iid, values=(
                r.get("well", ""), r.get("call", ""),
                f"{r.get('frac', 0.0):.3f}",
                f"{r.get('hom1', 0.0):.0f}", f"{r.get('hom2', 0.0):.0f}",
                f"{r.get('het1', 0.0):.0f}", f"{r.get('het2', 0.0):.0f}",
                f"{r.get('snr1', 0.0):.0f}", f"{r.get('snr2', 0.0):.0f}",
                f"{r.get('snr3', 0.0):.0f}", f"{r.get('snr4', 0.0):.0f}",
                ispeaks, r.get("sample_scans", ""), conf,
                r.get("flags", ""), r.get("reason", "")))
            tags = []
            if r.get("std_source") == "manual":
                tags.append("manual")
            if (r.get("flags") or "").strip():
                tags.append("flagged")
            if tags:
                tree.item(iid, tags=tuple(tags))
            self._auto_iids[(r.get("run"), r.get("well"))] = iid
        self._plate_map_refresh()
        self._below_refresh()
        self._below_follow()

    def _below_target(self):
        """Which pane belongs under the graphs right now, or None."""
        if getattr(self, "genotyping_active", False):
            return "pick"
        if getattr(self, "_below_mode", "geno") == "sequence":
            return "seq"
        if getattr(self, "_auto_table", None) is None:
            return "seq"
        shown_runs = {Path(p).parent.name
                      for p in (self.selected or [])[: self._graph_count()]}
        if self.always_table.get() or any(
                r.get("run") in shown_runs for r in (getattr(self, "_geno_rows", None) or [])):
            return "geno"
        # Nothing genotyped for these runs: the table would be an empty strip.
        # Show the called sequence if a shown well has one, else nothing.
        for p in (self.selected or [])[: self._graph_count()]:
            doc = self.docs.get(str(Path(p).resolve()))
            if doc is not None and getattr(doc, "sequence", ""):
                return "seq"
        return None

    def _set_below_visible(self, show):
        """Add or take away the lower pane, keeping the height it was given."""
        vs = getattr(self, "vsplit", None)
        if vs is None or show == self._below_shown:
            return
        if show:
            self.vsplit.add(self.below, weight=1)
            self._below_shown = True

            def _restore():
                try:
                    total = self.vsplit.winfo_height()
                    want = self._below_sash or max(120, int(total * 0.24))
                    self.vsplit.sashpos(0, max(80, total - want))
                except tk.TclError:
                    pass
            self.after_idle(_restore)
        else:
            try:
                total = self.vsplit.winfo_height()
                self._below_sash = max(60, total - self.vsplit.sashpos(0))
            except tk.TclError:
                pass
            self.vsplit.forget(self.below)
            self._below_shown = False

    def _below_refresh(self):
        """Pack exactly one pane under the graphs, and title the slot.

        Manual picking and area measuring take the slot for their own table;
        otherwise it is the genotyping table, or the basecalled sequence when
        the toggle asks for it.  One place decides this so entering a mode,
        leaving it and swapping the toggle can never leave two panes stacked.
        When none of them has anything to say, the slot is removed and the
        graphs get the room.
        """
        target = self._below_target()
        upd = getattr(self, "_update_mode_banner", None)
        if upd is not None:
            upd()
        self._set_below_visible(target is not None)
        pick = getattr(self, "pick_table", None)
        geno = getattr(self, "_auto_table", None)
        seq = getattr(self, "seq_text", None)
        want = {"pick": pick, "geno": geno, "seq": seq}.get(target)
        if want is None or getattr(self, "_below_packed", None) is not want:
            for w in (pick, geno, seq):
                if w is None or w is want:
                    continue
                try:
                    w.pack_forget()
                except tk.TclError:
                    pass
            if want is not None:
                try:
                    want.pack(fill=tk.BOTH, expand=True, padx=4, pady=2)
                except tk.TclError:
                    want = None
            self._below_packed = want
        title = {"pick": ("Measured areas" if getattr(self, "area_mode", False)
                          else "Picked peaks"),
                 "seq": "Called sequence", "geno": "Genotyping"}.get(target, "")
        hdr = getattr(self, "seq_hdr", None)
        if hdr is not None:
            try:
                hdr.config(text=title)
            except tk.TclError:
                pass
        btn = getattr(self, "seq_btn", None)
        if btn is not None:
            try:
                btn.config(text=("Show sequence" if getattr(
                    self, "_below_mode", "geno") == "geno"
                    else "Show genotyping"))
            except tk.TclError:
                pass

    def _toggle_below(self):
        """Swap the pane under the graphs between the genotyping table and the
        basecalled sequence.  Manual picking owns the slot while it is on."""
        if getattr(self, "genotyping_active", False):
            return
        self._below_mode = ("sequence"
                            if getattr(self, "_below_mode", "geno") == "geno"
                            else "geno")
        self._below_refresh()
        if self._below_mode == "sequence":
            self._show_sequence()
        else:
            self._below_follow()

    def _below_follow(self):
        """Move the table to the wells the graphs are showing.

        Paging through the traces should not leave the operator hunting for
        matching rows: the table scrolls to the batch on screen and marks it,
        so the two always describe the same wells.
        """
        tree = getattr(self, "_auto_tree", None)
        iids = getattr(self, "_auto_iids", None)
        if tree is None or not iids:
            return
        if getattr(self, "_below_mode", "geno") != "geno":
            return
        if getattr(self, "genotyping_active", False):
            return                       # the pick table owns this slot
        paths = list((getattr(self, "selected", None) or
                      [])[: self._graph_count()])
        shown = [iids[(p.parent.name, p.stem)]
                 for p in paths
                 if (p.parent.name, p.stem) in iids]
        if not shown or tuple(shown) == getattr(self, "_followed", None):
            return
        self._followed = tuple(shown)
        # Remember what this marking is, so the selection handler can tell the
        # table's own batch mark from an operator clicking a line: the two look
        # identical to <<TreeviewSelect>>, and mistaking one for the other
        # would drag the graphs to whichever row the batch happened to start
        # at every time the page moved.
        self._geno_programmatic = tuple(shown)
        self._geno_selecting = True
        try:
            tree.selection_set(shown)
            first = tree.index(shown[0]) if shown else 0
            total = max(1, len(tree.get_children()))
            tree.yview_moveto(float(first) / total)
            # keep selection
            tree.selection_set(shown)
        except tk.TclError:
            pass
        finally:
            self._geno_selecting = False

    def _geno_row_select(self, _evt=None):
        """Clicking a line in the table brings that well's trace up, first.

        The table and the graphs describe the same wells, so picking a line
        should show the trace it is about instead of making the operator hunt
        for it in the file list; with several traces on screen it goes to the
        top, which is the slot the eye is already on.  The rest of the window
        stays as it was -- the point is to look at one well more closely, not
        to turn the page.
        """
        if getattr(self, "_geno_selecting", False):
            return                      # our own marking, mid-flight
        if getattr(self, "genotyping_active", False):
            return                      # manual picking owns this slot
        tree = getattr(self, "_auto_tree", None)
        if tree is None:
            return
        sel = tuple(tree.selection())
        if not sel or sel == getattr(self, "_geno_programmatic", None):
            return                      # a batch mark, not an operator's click
        try:
            row = self._geno_row(int(sel[0]))
        except (TypeError, ValueError):
            return
        if row is None:
            return
        well = row.get("well", "")
        path = self._path_for_run_well(row.get("run"), well)
        if path is None:
            self.status_var.set(
                f"{well or 'That well'}'s trace is not open in this window.")
            return
        n = self._graph_count()
        current = [Path(p) for p in (getattr(self, "selected", None) or [])]
        if path in current:
            window = [path] + [p for p in current if p != path][:n - 1]
        else:
            files = [Path(p) for p in (getattr(self, "files", None) or [])]
            if path in files:
                k = files.index(path)
                window = [files[(k + j) % len(files)] for j in range(n)]
            else:
                window = [path] + current[:n - 1]
        if window == current:
            return                      # already at the top of the window
        self.selected = window
        # The file list has to agree with the graphs, or the next page starts
        # from where the old selection was and the two drift apart.
        idxs = []
        for p in window:
            for i, q in enumerate(getattr(self, "files", None) or []):
                if Path(q) == p:
                    idxs.append(i)
                    break
        if idxs:
            try:
                self.file_list.selection_clear(0, tk.END)
                for i in idxs:
                    self.file_list.selection_set(i)
                self.file_list.see(idxs[0])
            except tk.TclError:
                pass
            self.cycle_pos = idxs[0]
        self.status_var.set(f"{well}: top of the graphs "
                            f"({len(window)} shown)")
        self.redraw()

    def _geno_row(self, i):
        """The table's row *i*: the dict the table is showing.

        The table is this session's rows laid over the library's, so its index
        is not an index into ``_auto_rows`` -- editing through the session list
        reached the wrong well as soon as any library-only row was on screen.
        """
        rows = getattr(self, "_geno_rows", None) or []
        return rows[i] if 0 <= i < len(rows) else None

    def _adopt_row(self, row):
        """Take *row* into this session, so edits and Save can reach it.

        A row only the library has is a copy of what is on disk; changing it
        would change nothing anyone can reload.  Adopting swaps it in under
        its own (run, well), so the next save carries the operator's version
        instead of the engine's.
        """
        key = (row.get("run"), row.get("well"))
        rows = getattr(self, "_auto_rows", None)
        if rows is None:
            rows = self._auto_rows = []
        for j, existing in enumerate(rows):
            if (existing.get("run"), existing.get("well")) == key:
                rows[j] = row
                return row
        rows.append(row)
        return row

    def _path_for_run_well(self, run, well):
        """The open trace file for (run, well), or None when it is not loaded."""
        for p in list(self.selected or []) + list(self.files or []):
            p = Path(p)
            if p.parent.name == run and p.stem == well:
                return p
        return None

    def _geno_col(self, name):
        """Treeview column id for a named table column (``"#2"`` is call)."""
        cols = getattr(self, "_auto_cols", None) or ()
        if name not in cols:
            return ""
        return f"#{cols.index(name) + 1}"

    def _std_scans_of(self, row):
        """A row's standard as four ints, empty when it has none."""
        val = row.get("std_scans") or row.get("std_scans_manual") or ""
        parts = (val.replace(",", "/").split("/") if isinstance(val, str)
                 else list(val))
        try:
            return [int(float(p)) for p in parts if str(p).strip()][:4]
        except (TypeError, ValueError):
            return []

    def _auto_edit_call(self, event):
        """Double-click a table cell to correct it.

        The engine's call is a suggestion scored from one plate's noise and one
        well's standard; the operator's eye and the clinical context override
        it.  The override edits the row in place -- the same dict the Save
        button writes -- so a corrected call survives export instead of living
        only on screen.  A small dialog (not an inline widget) is used because
        the combobox popdown steals focus and would clobber an inline editor.

        Which column was hit picks the editor: the call, or the standard's four
        scans (which re-read the well against them).
        """
        if getattr(self, "_auto_tree", None) is None or \
                getattr(self, "_auto_edit_win", None) is not None:
            return
        rowid = self._auto_tree.identify_row(event.y)
        if not rowid:
            return
        i = int(rowid)
        if self._geno_row(i) is None:
            return
        colid = self._auto_tree.identify_column(event.x)
        if colid == self._geno_col("call"):
            self._auto_edit_call_at(i, event)
        elif colid == self._geno_col("ispeaks"):
            self._auto_edit_is_at(i, event)

    def _auto_edit_call_at(self, i, event=None):
        """Build the genotype-edit dialog for table row *i*.

        Positioned over the clicked cell when the table is really on screen;
        near the pointer otherwise (headless runs have no cell geometry)."""
        row = self._geno_row(i)
        if row is None or getattr(self, "_auto_edit_win", None) is not None:
            return
        from scorer import CALLS
        x = y = 0
        if event is not None:
            try:
                bbox = self._auto_tree.bbox(str(i), self._geno_col("call"))
                if bbox:
                    # bbox() is tree-relative; geometry() takes root-relative
                    # coordinates, so without the offset the dialog lands at
                    # the top-left corner of the screen instead of on the cell.
                    x = self._auto_tree.winfo_rootx() + bbox[0]
                    y = self._auto_tree.winfo_rooty() + bbox[1]
            except tk.TclError:
                pass
            if not x:
                x, y = int(event.x_root), int(event.y_root)
        win = tk.Toplevel(self)
        self._auto_edit_win = win
        win.title(f"Genotype \u2014 {row.get('well', '')}")
        win.transient(self)
        # Closing with the window-manager X has to release the editor slot,
        # or the next double-click is refused and the table looks dead.
        win.protocol("WM_DELETE_WINDOW", self._auto_edit_close)
        ttk.Label(win, text="Genotype:").pack(side=tk.LEFT, padx=(8, 4), pady=6)
        var = tk.StringVar(value=row.get("call", ""))
        om = ttk.Combobox(win, textvariable=var, values=list(CALLS),
                          state="readonly", width=12)
        om.pack(side=tk.LEFT, padx=(0, 8), pady=6)
        win._apply = lambda: self._auto_edit_apply(i, var)
        win._var = var
        btn = ttk.Button(win, text="Apply", command=win._apply)
        btn.pack(side=tk.LEFT, padx=2, pady=6)
        btn = ttk.Button(win, text="Cancel", command=self._auto_edit_close)
        btn.pack(side=tk.LEFT, padx=2, pady=6)
        # Enter commits the highlighted value and picking one applies at once,
        # the same single-gesture rhythm as a spreadsheet dropdown.
        om.bind("<Return>", lambda _e: win._apply())
        om.bind("<<ComboboxSelected>>", lambda _e: win._apply())
        win.bind("<Escape>", lambda _e: self._auto_edit_close())
        win.geometry(f"+{max(0, int(x))}+{max(0, int(y))}")
        om.focus_set()
        self._grab_when_viewable(win)

    def _auto_edit_apply(self, i, var):
        """Commit the operator's override into the row that Save writes.

        The library records whose call it is: an operator's answer is a label
        worth training on, the engine's guess only a prediction.
        """
        row = self._geno_row(i)
        if row is None:
            self._auto_edit_close()
            return
        row["call"] = var.get()
        self._adopt_row(row)
        self._library_store_results([row], label_source="edited")
        col = self._geno_col("call")
        if getattr(self, "_auto_tree", None) is not None and col:
            self._auto_tree.set(str(i), col, row["call"])
        self._auto_edit_close()

    def _auto_edit_close(self):
        win = getattr(self, "_auto_edit_win", None)
        if win is None:
            return
        self._auto_edit_win = None
        try:
            win.grab_release()
        except tk.TclError:
            pass
        try:
            win.destroy()
        except tk.TclError:
            pass

    def _auto_edit_is_at(self, i, event=None):
        """Build the standard-editor dialog for table row *i*.

        The standard is four scan numbers -- one per band -- which is exactly
        how it reads on the trace, so the operator can correct it straight from
        the peaks they are looking at instead of re-picking all four by hand.
        """
        row = self._geno_row(i)
        if row is None or getattr(self, "_auto_edit_win", None) is not None:
            return
        scans = self._std_scans_of(row)
        x = y = 0
        if event is not None:
            try:
                bbox = self._auto_tree.bbox(str(i), self._geno_col("ispeaks"))
                if bbox:
                    # bbox() is tree-relative; geometry() is root-relative,
                    # so the cell's own position needs the window offset or
                    # the dialog opens at the screen's top-left corner.
                    x = self._auto_tree.winfo_rootx() + bbox[0]
                    y = self._auto_tree.winfo_rooty() + bbox[1]
            except tk.TclError:
                pass
            if not x:
                x, y = int(event.x_root), int(event.y_root)
        win = tk.Toplevel(self)
        self._auto_edit_win = win
        win.title(f"Standard \u2014 {row.get('well', '')}")
        win.transient(self)
        win.protocol("WM_DELETE_WINDOW", self._auto_edit_close)
        ttk.Label(win, text="IS scans:").pack(side=tk.LEFT, padx=(8, 4),
                                              pady=6)
        vars_ = []
        first = None
        for n in range(4):
            var = tk.StringVar(value=str(scans[n]) if n < len(scans) else "")
            entry = ttk.Entry(win, textvariable=var, width=7)
            entry.pack(side=tk.LEFT, padx=2, pady=6)
            if first is None:
                first = entry
            vars_.append(var)
        win._is_vars = vars_
        win._apply = lambda: self._auto_edit_is_apply(i, vars_)
        btn = ttk.Button(win, text="Apply", command=win._apply)
        btn.pack(side=tk.LEFT, padx=2, pady=6)
        btn = ttk.Button(win, text="Cancel", command=self._auto_edit_close)
        btn.pack(side=tk.LEFT, padx=2, pady=6)
        win.bind("<Return>", lambda _e: win._apply())
        win.bind("<Escape>", lambda _e: self._auto_edit_close())
        if event is not None and not x:
            x, y = int(event.x_root), int(event.y_root)
        win.geometry(f"+{max(0, int(x))}+{max(0, int(y))}")
        if first is not None:
            first.focus_set()
        self._grab_when_viewable(win)

    def _grab_soon(self, win, tries=60):
        """``grab_set`` for any popup, retried until the window is mapped.

        Same XWayland "window not viewable" story as ``_grab_when_viewable``,
        without tying the popup to the table editor's slot."""
        try:
            if not win.winfo_exists():
                return
            win.grab_set()
        except tk.TclError:
            if tries > 0:
                try:
                    win.after(50, lambda: self._grab_soon(win, tries - 1))
                except tk.TclError:
                    pass

    def _grab_when_viewable(self, win, tries=60):
        """Make a freshly built popup modal once the window manager maps it.

        ``grab_set()`` on a brand-new ``Toplevel`` raises
        ``TclError: grab failed: window not viewable`` under XWayland because
        the window is not on screen yet.  That error used to escape *before*
        the dialog's widgets were built, leaving an empty grey box with no
        controls; and because the editor slot stayed claimed, every later
        double-click was refused.  Retrying from the event loop keeps the
        original modal behaviour without depending on map timing, and a popup
        that never becomes viewable still works, just non-modally.
        """
        if getattr(self, "_auto_edit_win", None) is not win:
            return
        try:
            if not win.winfo_exists():
                return
        except tk.TclError:
            return
        try:
            win.grab_set()
        except tk.TclError:
            if tries > 0:
                try:
                    win.after(50, lambda: self._grab_when_viewable(win,
                                                                   tries - 1))
                except tk.TclError:
                    pass

    def _auto_edit_is_apply(self, i, vars_):
        """Re-read the well against the four scans the operator typed.

        The scans become the well's standard mark -- the picker carries them
        so the next Auto-genotype defers to them, the library keeps them so a
        restart does not lose them -- and the well is scored again against
        exactly those peaks, which is what makes the row's call and features
        describe the standard now on screen rather than the one before.
        """
        row = self._geno_row(i)
        if row is None:
            self._auto_edit_close()
            return
        try:
            scans = [int(round(float(v.get()))) for v in vars_]
        except (TypeError, ValueError):
            scans = []
        if len(scans) != 4 or any(s < 0 for s in scans) or \
                scans != sorted(scans) or len(set(scans)) != 4:
            messagebox.showerror(
                "Standard",
                "A standard is four different, increasing scan numbers,\n"
                "one per band: hom1, het1, het2, hom2.")
            return
        run, well = row.get("run"), row.get("well")
        path = self._path_for_run_well(run, well)
        if path is None:
            messagebox.showerror(
                "Standard",
                f"{well or 'This well'}'s trace is not open in this window.\n"
                "Open the run first, then correct the standard.")
            return
        try:
            doc = self._ensure_doc(path)
            pk = self._ensure_picker(path)
        except Exception as exc:                       # noqa: BLE001
            messagebox.showerror("Standard", f"{well}: {exc}")
            return
        try:
            pk.std = [(s, "") for s in scans]
            pk.std_auto = False                        # the operator's mark
            self._library_store(path, pk)
        except Exception as exc:                       # noqa: BLE001
            print(f"library: could not save {run}/{well}: {exc}", file=sys.stderr)
        import genotyping
        models = getattr(self, "_auto_models", None) or {}
        try:
            new = genotyping.auto_genotype(
                doc, is_channel=int(self._auto_is_channel.get()),
                sample_channel=int(self._auto_sample_channel.get()),
                base_order=self._doc_order(doc),
                cut=genotyping.DEFAULT_IS_CUT, run_name=run,
                std_scans_manual=list(scans),
                is_model=models.get(run) or getattr(self, "_auto_is_model",
                                                    None),
                min_het_frac=float(self._auto_min_het_frac.get()))
        except Exception as exc:                       # noqa: BLE001
            messagebox.showerror("Standard", f"{well}: {exc}")
            return
        self._adopt_row(new)
        self._library_store_results([new], label_source="edited")
        if getattr(self, "_auto_tree", None) is not None:
            self._sync_geno_table()
        self.status_var.set(
            f"{well}: standard set to " + "/".join(str(s) for s in scans))
        self._auto_edit_close()

    def auto_genotype_save(self):
        """Write the auto-genotype rows out through the shared save_table, so
        they land in the same CSV/XLSX/JSON shapes as the manual pick table."""
        if not self._auto_rows:
            messagebox.showinfo(
                "Save auto-genotype table",
                "Run Auto genotyping ▸ Auto-genotype selected "
                "wells\nfirst.")
            return
        from genotyping import save_table
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
        try:
            save_table(Path(path), self._auto_rows)
        except Exception as e:
            messagebox.showerror("Save auto-genotype table", str(e))
            return
        self.status_var.set(f"Saved {len(self._auto_rows)} auto-genotype "
                            f"rows to {Path(path).name}.")

    def _gen_save(self):
        from genotyping import save_table
        rows = []
        if self.area_mode:
            # run folder first, then well, then left-to-right along the trace
            for key, pk in sorted(self._area_pickers.items(),
                                  key=lambda kv: (kv[1].run_name, kv[1].well,
                                                  str(kv[1].path))):
                if pk.records:
                    rows.extend(pk.export_rows())
        else:
            # group by sample first (run folder, then well name) so one
            # sample's peaks are never interleaved with another's, then scan
            # order inside
            for key, pk in sorted(self._gen_pickers.items(),
                                  key=lambda kv: str(kv[1].path)):
                if pk.records:
                    rows.extend(pk.export_rows())
        if not rows:
            messagebox.showinfo(
                "Save peaks table" if not self.area_mode else "Save areas",
                "Pick some peaks first, then save." if not self.area_mode
                else "Measure something first, then save.")
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
        try:
            save_table(path, rows)
        except Exception as e:
            messagebox.showerror("Save peaks table",
                                 f"Could not write file:\n{e}")
            return
        wells = sorted({r["well"] for r in rows})
        unit = "area" if self.area_mode else "peak"
        self.status_var.set(f"Saved {len(rows)} {unit} rows ({len(wells)} well"
                            f"{'s' if len(wells) != 1 else ''}) to {path}")

    def size_selected_fragments(self):
        """Size the selected wells against a size-standard ladder.

        Sizing is a *different* analysis from peak picking -- no internal
        standard, no prior sample knowledge -- so it works on any selection
        and writes a fresh sized-peak table (CSV/JSON) or an HTML report.  The
        channels are physical (1-4 in the run's dye order) and are asked for
        explicitly because the ladder/sample pair is a property of the kit.
        """
        if not self.selected:
            messagebox.showinfo("Fragment sizing",
                                "Select one or more wells first.")
            return
        from fragment_sizing import list_ladders, load_ladder, size_trace

        # the last ladder and channels used: a lab sizes plate after plate
        # with the same kit, so retyping them each time is only a chance to
        # get one wrong
        try:
            last = dict(mark_library.load_prefs().get("fragment_sizing", {}))
        except Exception:                                  # noqa: BLE001
            last = {}
        names = ", ".join(list_ladders())
        spec = simpledialog.askstring(
            "Fragment-length sizing",
            "Ladder: a built-in name,\n"
            f"  {names},\n"
            "a comma-separated length list, or a path to a ladder JSON file.\n"
            "Run `python scorer.py ladders` for the MegaBACE ET standards.",
            initialvalue=str(last.get("ladder") or "genescan500_rox"),
            parent=self)
        if spec is None:
            return
        try:
            ladder = load_ladder(spec.strip())
        except (ValueError, OSError) as e:                 # noqa: BLE001
            messagebox.showerror("Fragment sizing", f"Ladder:\n{e}")
            return
        ladder_channel = simpledialog.askinteger(
            "Fragment sizing", "Ladder channel (1-4, physical):",
            initialvalue=int(last.get("ladder_channel") or 4),
            minvalue=1, maxvalue=4, parent=self)
        if ladder_channel is None:
            return
        sample_channel = simpledialog.askinteger(
            "Fragment sizing", "Sample channel (1-4, physical):",
            initialvalue=int(last.get("sample_channel") or 2),
            minvalue=1, maxvalue=4, parent=self)
        if sample_channel is None:
            return
        try:
            prefs = mark_library.load_prefs()
            prefs["fragment_sizing"] = {"ladder": spec.strip(),
                                        "ladder_channel": ladder_channel,
                                        "sample_channel": sample_channel}
            mark_library.save_prefs(prefs)
        except Exception:                                  # noqa: BLE001
            pass
        out = filedialog.asksaveasfilename(
            parent=self, defaultextension=".csv",
            filetypes=[("CSV (sized peaks)", "*.csv"), ("JSON", "*.json"),
                       ("HTML report", "*.html")])
        if not out:
            return

        from scorer import sizing_html, write_rows

        results, failures = [], []
        for path in list(self.selected):
            try:
                doc = self._ensure_doc(Path(path))
                results.append(size_trace(
                    doc, ladder, ladder_channel=ladder_channel,
                    sample_channel=sample_channel,
                    base_order=self._doc_order(doc)))
            except Exception as e:                          # noqa: BLE001
                failures.append(f"{Path(path).name}: {e}")
        if not results:
            messagebox.showerror(
                "Fragment sizing",
                "No wells could be sized." +
                ("\n" + "\n".join(failures) if failures else ""))
            return
        try:
            if out.lower().endswith((".html", ".htm")):
                Path(out).write_text(sizing_html(results), encoding="utf-8")
            else:
                write_rows(out, [r for res in results for r in res.rows])
        except Exception as e:                              # noqa: BLE001
            messagebox.showerror("Fragment sizing",
                                 f"Could not write file:\n{e}")
            return

        n_peaks = sum(len(res.rows) for res in results)
        n_warn = sum(len(res.warnings) for res in results)
        warned = [res.well or res.file for res in results if res.warnings]
        if len(warned) > 12:
            warned = warned[:12] + [f"… {len(warned) - 12} more"]
        self.status_var.set(
            f"Sized {n_peaks} peak(s) in {len(results)} well(s) → {out}")
        messagebox.showinfo(
            "Fragment sizing",
            f"Sized {n_peaks} peak(s) across {len(results)} well(s) "
            f"against {ladder.name}.\n" +
            (f"{n_warn} warning(s) in well(s) {', '.join(warned)} — their "
             "rows have ladder_ok = False; open the report or run the "
             "sizing CLI for details.\n" if n_warn else "") +
            (f"{len(failures)} well(s) skipped.\n" if failures else "") +
            f"Written to {out}")

    def _col_to_chan(self, order=None):
        """acgt matrix column -> physical channel index of that trace.

        The matrix is always stored A,C,G,T, but the *instrument* channel a
        dye sits on follows the run's base order, so physical Ch1 is not
        necessarily the A column.  Everything that keys off a channel number
        (the show/hide checkboxes, the labels) must go through this map."""
        order = (order or self.base_order_var.get() or "ACTG").upper()
        out: dict = {}
        for ci, base in enumerate(order[:4]):
            if base in CHANNEL_ORDER:
                out[CHANNEL_ORDER.index(base)] = ci
        return out

    def _on_base_order_changed(self, *_args):
        """Dye order changed: the channel bar, the legends and the peak colours
        all follow it.  Pickers and traces baked in the old order, so drop them
        -- but never silently throw away peaks the user has already picked."""
        self._refresh_channel_labels()
        if any(pk.records for pk in self._gen_pickers.values()):
            self.status_var.set(
                f"Dye order is now {self.base_order_var.get()} — the already "
                "picked peaks keep the colours they were picked with. Clear "
                "picks to rebuild them in the new order.")
            return
        self._gen_pickers = {}
        # Every cached trace mapped its columns with the old order, and the
        # cache is keyed by path alone, so reload them.  Area pickers hold the
        # same docs; keep them (and their docs) only when they carry regions.
        if not any(getattr(pk, "records", None)
                   for pk in self._area_pickers.values()):
            self._area_pickers = {}
            self.docs.clear()
        self.redraw()

    def channel_bar_order(self, order=None):
        """The channel checkboxes as (channel, label), left to right.

        Widget k drives chan_show[k], and chan_show is indexed by *physical
        channel*, so widget k is already Ch k+1 and the widgets are already in
        the right order.  The only thing that was ever wrong is the label: it
        was keyed off the matrix column the widget was born with rather than
        off the channel, which is what put "Ch4" to the left of "Ch3".  So this
        relabels in place and never re-packs -- the checkboxes share their
        parent with side=RIGHT siblings, and re-packing them would move the
        whole group to the end of the toolbar, away from its "Channels:"
        label.  Free of Tk so it can be tested.
        """
        order = (order or self._display_order() or "ACTG").upper()[:4]
        return [(ci, f"  Ch{ci + 1} {order[ci]}  ")
                for ci in range(min(len(order), len(CHANNEL_ORDER)))]

    def _refresh_channel_labels(self):
        """Relabel the channel checkboxes to match the run's dye order."""
        cbs = getattr(self, "_chan_cbs", [])
        for ci, text in self.channel_bar_order():
            if ci < len(cbs):
                try:
                    cbs[ci].config(text=text)
                except Exception:
                    pass

    def _chan_changed(self):
        """Channel on/off toggles refresh the viewer (and the picked-peak
        overlays when main-window picking is active)."""
        self.redraw()

    def _show_sequence(self):
        self.seq_text.delete("1.0", tk.END)
        for path in self.selected[: self._graph_count()]:
            key = str(path.resolve())
            doc = self.docs.get(key)
            if not doc:
                continue
            if key not in self._comment_cache:
                self._comment_cache[key] = self._read_comment(path)
            comment = self._comment_cache.get(key, "")
            if not doc.sequence:
                head = f">{doc.well} ({doc.source}) — use Basecall ▸ Basecall selected"
                if comment:
                    head += f"\n; {comment}"
                self.seq_text.insert(tk.END, head + "\n\n")
                continue
            head = f">{doc.well} {self._well_stats(doc)}  [{doc.source}]"
            if comment:
                head += f"\n; {comment}"
            self.seq_text.insert(
                tk.END, f"{head}\n{doc.sequence}\n\n")
        self._below_follow()       # the table follows the batch either way

    # ------------------------------------------------------- comments / info
    def _comment_file(self, path: Path) -> Path:
        """Comments are kept next to the data as '<file>.comment.txt' — they
        travel with the run and never touch the binary .rsd/.scf file."""
        return Path(str(path) + ".comment.txt")

    def _read_comment(self, path: Path) -> str:
        cf = self._comment_file(path)
        try:
            if cf.exists():
                return cf.read_text(encoding="utf-8", errors="replace").strip()
        except OSError:
            pass
        return ""

    def edit_comments(self):
        """View / write the run comment for the currently selected file."""
        if not self.selected:
            messagebox.showinfo("Run comments", "Select a well/file first.")
            return
        path = self.selected[0]
        key = str(path.resolve())
        if key not in self._comment_cache:
            self._comment_cache[key] = self._read_comment(path)
        comment = self._comment_cache.get(key, "")

        win = tk.Toplevel(self)
        win.title("Run comments")
        win.geometry("520x300")
        win.transient(self)
        self._grab_soon(win)
        tk.Label(win, text=f"{path.name}: comment for this run", anchor=tk.W).pack(
            fill=tk.X, padx=10, pady=(8, 2))
        box = tk.Text(win, wrap=tk.WORD, undo=True)
        box.insert("1.0", comment)
        box.pack(fill=tk.BOTH, expand=True, padx=10, pady=4)
        note = (
            "Stored beside the data file:\n"
            f"{self._comment_file(path).name}\n"
            "You can delete the file to remove the comment.")
        tk.Label(win, text=note, fg="#555", justify=tk.LEFT, anchor=tk.W).pack(
            fill=tk.X, padx=10, pady=(0, 2))
        bar = ttk.Frame(win)
        bar.pack(fill=tk.X, padx=10, pady=6)

        def save():
            text = box.get("1.0", "end-1c").strip()
            self._comment_cache[key] = text
            try:
                if text:
                    self._comment_file(path).write_text(text, encoding="utf-8")
                else:
                    try:
                        self._comment_file(path).unlink()
                    except OSError:
                        pass
                self.status_var.set(f"Comment saved for {path.name}")
            except OSError as e:
                messagebox.showerror("Run comments", f"Could not write comment:\n{e}")
            win.destroy()
            self._show_sequence()

        ttk.Button(bar, text="Save", command=save).pack(side=tk.LEFT)
        ttk.Button(bar, text="Discard", command=win.destroy).pack(side=tk.LEFT, padx=6)

    def show_run_info(self):
        """Read-only rundown of the selected well's run parameters."""
        if not self.selected:
            messagebox.showinfo("Run info", "Select a well/file first.")
            return
        path = self.selected[0]
        key = str(path.resolve())
        try:
            doc = self._ensure_doc(path)
        except Exception as e:
            messagebox.showerror("Run info", f"{path.name}:\n{e}")
            return
        comment = self._comment_cache.get(key) or self._read_comment(path)

        add = lambda k, v: chunks.append(f"{k:<20} {v}")
        chunks = ["Run info", "   " + str(path), ""]
        add("Well", doc.well)
        add("Source", doc.source)
        add("Scans", f"{doc.n_scans}  ({doc.n_scans / SCAN_RATE_HZ / 60:.1f} min "
                     f"at {SCAN_RATE_HZ:.2f} Hz)")
        # The instrument's own run fields, decoded from the RSD footer. These
        # are the record of what was actually run: plate, sample, machine,
        # method and base caller. Anything already in `meta` is redundant now.
        fld = getattr(doc, "footer_fields", None) or {}
        if fld:
            chunks.append("")
            add("Run", fld.get("R CODE") or fld.get("PLATE ID") or "-")
            for key in ("PLATE ID", "SAMPLE NAME", "WELL ID", "MACHINE ID",
                        "APPLICATION", "BASE CALLER", "NAME", "COMMENT",
                        "LASER MODE", "BEAMSPLITTER A", "BEAMSPLITTER B",
                        "FILTER"):
                if fld.get(key):
                    add(key, fld[key])
            for key, label in (("CHANNEL DYES", "Dyes (Ch1–Ch4)"),
                               ("CHANNEL BASES", "Bases (Ch1–Ch4)")):
                vals = fld.get(key) or []
                if any(vals):
                    add(label, " · ".join(v or "–" for v in vals))
        src = {"run": "chosen for this run (Comments ▸ Dye order of this run…)",
               "file": "declared by the file",
               "default": "general setting (⚙ Basecall settings…)"}.get(
                   getattr(doc, "dye_order_source", "default"), "")
        add("Dye order (channels)", f"{self._doc_order(doc)}  — {src}")
        # Run voltage / time / temperature / PMT come out of the .rsd itself;
        # a text export is merged in only for the extras it adds.
        rp = getattr(doc, "run_params", None) or {}
        rp_rows = run_params_summary(rp)
        chunks.append("")
        if rp_rows:
            chunks.append("Run conditions (read from the .rsd)")
            for k, v in rp_rows:
                if "\n" in v:
                    first, *rest = v.splitlines()
                    add(k, first)
                    for extra in rest:
                        chunks.append(f"{'':<20} {extra}")
                else:
                    add(k, v)
            # The export is also the authoritative per-run dye order; flag a
            # disagreement with the current setting rather than hiding it.
            exp_order = ""
            for key, val in rp.items():
                if key.strip().lower() == "base order":
                    exp_order = val.strip().upper()
                    break
            cur = self._doc_order(doc)
            if exp_order and exp_order != cur:
                add("Dye order note", f"run says {exp_order}, using {cur}")
            nl = ""
            for key, val in rp.items():
                if key.strip().lower() == "number of lines" and val.strip().isdigit():
                    nl = val.strip()
                    break
            if nl and int(nl) != doc.n_scans:
                add("Scan count note",
                    f"export says {nl} lines, trace has {doc.n_scans}")
        else:
            chunks.append("Run conditions")
            add("Voltage / time / temp", "not recoverable from this file")
        add("Column layout", doc.base_order)
        add("Basecaller", self.basecaller.get())
        if doc.sequence:
            n = len(doc.sequence)
            qa = np.asarray(doc.qualities[:n], dtype=float)
            sp = np.median(np.diff(np.asarray(doc.peak_positions[:n], float)))
            add("Sequence length", str(n))
            add("Quality mean", f"{qa.mean():.1f}" if qa.size else "-")
            add("Quality min", f"{qa.min():.0f}" if qa.size else "-")
            add("N count", str(doc.sequence.count("N")))
            add("Peak spacing", f"{sp:.2f} scans" if np.isfinite(sp) else "-")
        else:
            add("Sequence", "not base-called yet")
        tr = doc.acgt
        if tr.size:
            # Report per *channel*, not per matrix column: doc.base_order is
            # the A,C,G,T matrix layout, so pairing it with tr[:, c] would
            # label the columns as channels and disagree with the channel bar.
            order = self._doc_order(doc)[:4]
            col2chan = self._col_to_chan(order)
            per = ", ".join(
                f"Ch{ci + 1} {order[ci]}={float(np.nanmax(tr[:, col])):.3f} V"
                for col, ci in sorted(col2chan.items(), key=lambda kv: kv[1]))
            add("Signal max (V)", per)
        cu = doc.current_ua
        if cu is not None and cu.size:
            add("Current (µA)", f"mean {np.nanmean(cu):.2f}, max {np.nanmax(cu):.2f}")
        if comment:
            chunks.append("")
            chunks.append("Comment:")
            chunks.append(comment)

        win = tk.Toplevel(self)
        win.title("Run info")
        win.geometry("680x420")
        win.transient(self)
        self._grab_soon(win)
        txt = scrolledtext.ScrolledText(
            win, wrap=tk.WORD, font=("DejaVu Sans Mono", 9))
        txt.insert("1.0", "\n".join(chunks))
        txt.config(state=tk.DISABLED)
        txt.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)
        ttk.Button(win, text="Close", command=win.destroy).pack(anchor=tk.E, padx=10, pady=6)

    # ------------------------------------------------------------------ plot
    def _theme_colors(self) -> dict:
        return dict(TRACE_THEMES.get(self.theme.get(), "Classic"))

    def _theme_mode(self) -> str:
        """'base' themes paint each trace by the colour of its base letter
        (MegaBACE SA sequencing colours are fixed per base); 'channel' themes
        paint Channel1..4 by position (genotyping dye sets)."""
        return TRACE_THEME_MODE.get(self.theme.get(), "base")

    def _channel_lines(self, settings):
        """Physical channel order Ch1..Ch4 -> (channel_idx, base letter,
        acgt column, colour), so traces render exactly as MegaBACE SA shows
        them: channels in run order, each coloured by the base it carries."""
        order = (settings.base_order or "ACTG").upper()
        colors = self._theme_colors()
        mode = self._theme_mode()
        lines = []
        for ci, base in enumerate(order[:4]):
            try:
                col = CHANNEL_ORDER.index(base)
            except ValueError:
                continue
            color = colors[CHANNEL_ORDER[ci]] if mode == "channel" else colors.get(base, "#000000")
            lines.append((ci, base, col, color))
        return lines

    def _refresh_channel_buttons(self):
        """Label/colour the Channel checkbuttons in the run's physical channel
        order (Ch1 = first tick, always shown first), using the colour of the
        base each channel carries (Classic: Ch1 green, Ch2 blue, Ch3 red, Ch4
        black, whatever the chemistry)."""
        cbs = getattr(self, "_chan_cbs", None)
        if not cbs:
            return
        order = (self._display_order() or "ACTG").upper()
        colors = self._theme_colors()
        mode = self._theme_mode()
        for ci, base in enumerate(order[:4]):
            if ci >= len(cbs):
                break
            color = colors[CHANNEL_ORDER[ci]] if mode == "channel" else colors.get(base, "#000000")
            cbs[ci].config(text=f"  Ch{ci + 1} {base}  ", fg=color)

    def _well_stats(self, doc) -> str:
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
            dset = self._settings_for(doc, settings)
            tr = display_trace(doc, dset)
            n = tr.shape[0]
            gx1 = max(gx1, n)
            segs = np.array_split(np.arange(n), rows)
            for r in range(rows):
                ax = axes[k * rows + r]
                sel = segs[r]
                for ci, base, col, color in self._channel_lines(dset):
                    if self.chan_show[ci].get():
                        y = tr[sel, col]
                        ax.plot(sel, y, color=color, lw=0.6,
                                label=f"Ch{ci + 1} {base}" if r == 0 else None)
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
                if self.show_letters.get():
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
            if draw_let:
                ax.text(p, ytop + 0.015 * ytop, base_, color=colors[base_],
                        ha="center", va="bottom", fontsize=6.5, zorder=5, clip_on=True)
        return med_sp

    def _draw_esd_peaks(self, ax, doc, s0, s1, channel_lines):
        """ESD-style view — the MegaBACE basecaller's processed picture: the
        four dye traces with a clear mark on every peak.  When the run is
        base-called the marks sit at the CALLER's peaks, coloured like the
        trace that carries them; a heterozygous position then shows as two
        rings at one scan and a mutation as a ring in an unexpected channel.
        Nothing called yet?  Falls back to simple per-channel peak detection
        so the pattern still reads like an ESD."""
        tr = doc.acgt
        if tr is None or tr.size == 0 or s1 <= s0:
            return
        col_color = {col: color for ci, base, col, color in channel_lines}
        colors = self._theme_colors()
        if doc.sequence and doc.peak_positions:
            pos = np.asarray(doc.peak_positions, dtype=float)
            for pi, base_ in enumerate(doc.sequence):
                if pi >= len(pos):
                    break
                if base_ not in CHANNEL_ORDER:
                    continue
                p = int(pos[pi])
                if not (s0 <= p < s1) or p >= tr.shape[0]:
                    continue
                col = CHANNEL_ORDER.index(base_)
                color = col_color.get(col, colors.get(base_, "#444"))
                ax.scatter([p], [tr[p, col]], s=18, marker="o",
                           facecolors="none", edgecolors=color, linewidths=0.8,
                           zorder=5, rasterized=True, clip_on=True)
            return
        seg = tr[s0:s1]
        for ci, base, col, color in channel_lines:
            y = seg[:, col]
            if y.size < 3:
                continue
            y0 = float(np.nanmin(y))
            y1 = float(np.nanmax(y))
            if not np.isfinite(y1 - y0) or y1 - y0 <= 1e-9:
                continue
            rng = y1 - y0
            # strict local maxima lifted off the noise floor (>5% of range)
            cand = np.flatnonzero(
                (y[1:-1] >= y[:-2]) & (y[1:-1] > y[2:])
                & (y[1:-1] - y0 >= 0.05 * rng))
            cand += 1
            keep = []
            for p in cand:
                if not keep or p - keep[-1] >= 3:
                    keep.append(p)
            if not keep:
                continue
            kp = np.asarray(keep, dtype=int)
            peak_v = y[kp]
            ax.scatter(s0 + kp, peak_v, s=14, marker="o",
                       facecolors="none", edgecolors=color, linewidths=0.7,
                       zorder=5, rasterized=True, clip_on=True)

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

    def _wheel_zoom(self, axis, direction, anchor=None):
        """Zoom one axis in (direction > 0) or out, keeping *anchor* -- a
        fraction of the full range, e.g. the point under the mouse -- where it
        is on screen.  Without an anchor the view zooms about its centre."""
        v = self._view_x if axis == "x" else self._view_y
        lo, hi = v
        width = max(1e-9, hi - lo)
        new_w = min(1.0, max(width / 1000.0, width * (0.8 if direction > 0 else 1.25)))
        if anchor is None:
            anchor = 0.5 * (lo + hi)
        rel = 0.5 if width <= 0 else min(1.0, max(0.0, (anchor - lo) / width))
        lo2 = anchor - rel * new_w
        lo2 = min(max(0.0, lo2), 1.0 - new_w)
        v[:] = [lo2, lo2 + new_w]
        self._apply_zoom()

    def _graph_count(self) -> int:
        """The Graphs box, clamped to 1..MAX_GRAPHS (and safe while the box
        is being edited and holds no number)."""
        try:
            n = int(self.n_graphs.get())
        except (tk.TclError, ValueError):
            n = 1
        return max(1, min(MAX_GRAPHS, n))

    def _save_view_prefs(self, *_):
        try:
            prefs = mark_library.load_prefs()
            prefs["view"] = {"y_scale": self.y_scale.get(),
                             "fit_y": bool(self.fit_y.get()),
                             "overlay_runs": bool(self.overlay_runs.get()),
                             "always_table": bool(self.always_table.get())}
            mark_library.save_prefs(prefs)
        except Exception:                                  # noqa: BLE001
            pass

    def _view_option_changed(self, redraw=False):
        self._save_view_prefs()
        if redraw:
            self.redraw()
        else:
            self._apply_zoom()

    @staticmethod
    def _trace_yrange(ax, xlo, xhi):
        """Min/max of the trace lines (labelled "Ch…") drawn on *ax* within
        [xlo, xhi], or None when there are none.  X data is the scan axis,
        ascending, so the window is a slice, not a scan of every point."""
        lo = hi = None
        for line in ax.get_lines():
            if not str(line.get_label()).startswith("Ch"):
                continue
            xd = np.asarray(line.get_xdata(), dtype=float)
            yd = np.asarray(line.get_ydata(), dtype=float)
            if xd.size == 0:
                continue
            a = int(np.searchsorted(xd, xlo, side="left"))
            b = int(np.searchsorted(xd, xhi, side="right"))
            seg = yd[a:b]
            seg = seg[np.isfinite(seg)] if seg.size else seg
            if seg.size == 0:
                continue
            m0, m1 = float(seg.min()), float(seg.max())
            lo = m0 if lo is None else min(lo, m0)
            hi = m1 if hi is None else max(hi, m1)
        return None if lo is None else (lo, hi)

    def _y_references(self, xlo, xhi):
        """Per-graph Y range the vertical view fraction applies to.

        Shared: every graph uses one range, so heights compare across wells.
        Each graph: every graph uses its own, so a weak well is not flattened
        by a strong one.  With "Fit Y to view" the range is taken over the
        scans on screen, so zooming into a region also scales to its peaks.
        """
        full = self._full_ylim
        if self._wrap_active:
            return [full] * len(self._plot_axes)
        if not self.fit_y.get():
            xlo, xhi = self._full_xlim
        # each graph over its OWN scans for that reference window (alignment)
        ranges = [self._trace_yrange(ax, self._from_ref(i, xlo), self._from_ref(i, xhi))
                  for i, ax in enumerate(self._plot_axes)]
        out = []
        if self.y_scale.get() == "each":
            for r in ranges:
                out.append(full if r is None else r)
        else:
            got = [r for r in ranges if r is not None]
            shared = (min(r[0] for r in got), max(r[1] for r in got)) if got else None
            out = [shared or full] * len(self._plot_axes)
        padded = []
        for r in out:
            lo, hi = r
            pad = 0.04 * (hi - lo) or 1.0
            padded.append((lo - pad * 0.5, hi + pad))
        return padded

    def _zoom_to_scans(self, a, b):
        """Show scans a..b (on the full X range of the graphs)."""
        if self._full_xlim is None:
            return
        x0, x1 = self._full_xlim
        a, b = sorted((float(a), float(b)))
        a, b = max(x0, a), min(x1, b)
        if b - a < 2:
            return
        span = x1 - x0
        self._view_x[:] = [(a - x0) / span, (b - x0) / span]
        self._apply_zoom()

    def _scan_entry_go(self, *_):
        try:
            a = float(self.scan_from.get())
            b = float(self.scan_to.get())
        except (tk.TclError, ValueError):
            self.status_var.set("Scan range: type two scan numbers, e.g. 1800 and 2700")
            return "break"
        k = len(self._plot_axes) - 1               # typed in the bottom graph's scans
        self._zoom_to_scans(self._to_ref(k, a), self._to_ref(k, b))
        self.canvas.get_tk_widget().focus_set()
        return "break"

    # --- mouse on the graphs: wheel zooms at the cursor, drag zooms a box ---
    def _on_canvas_resize(self, _event=None):
        """Margins were fitted for the old size; refit once resizing pauses,
        so axis labels are not clipped after the table pane is dragged."""
        job = getattr(self, "_relayout_job", None)
        if job is not None:
            try:
                self.after_cancel(job)
            except tk.TclError:
                pass
        self._relayout_job = self.after(150, self._relayout)

    def _relayout(self):
        self._relayout_job = None
        if not self._plot_axes:
            return
        try:
            if self._wrap_active:
                self.fig.tight_layout(rect=(0.045, 0, 1, 1), h_pad=0.3, w_pad=0.2, pad=0.4)
            else:
                self.fig.tight_layout(rect=(0.045, 0, 1, 1), h_pad=0.25, pad=0.3)
            self.canvas.draw_idle()
        except Exception:                                  # noqa: BLE001
            pass

    def _plot_hit(self, event):
        """(graph axes, x, y) under the mouse, in that graph's data units.

        The µA overlay is a twin axes drawn on top, so matplotlib reports the
        twin as ``inaxes``; resolve by position to the graph underneath."""
        for ax in self._plot_axes:
            try:
                if ax.get_window_extent().contains(event.x, event.y):
                    xd, yd = ax.transData.inverted().transform((event.x, event.y))
                    return ax, float(xd), float(yd)
            except Exception:                              # noqa: BLE001
                continue
        return None

    def _on_plot_scroll(self, event):
        """Wheel: zoom X at the cursor.  Ctrl+wheel: zoom Y.  Shift+wheel: pan."""
        hit = self._plot_hit(event) if self._full_xlim is not None else None
        if hit is None:
            return
        hax, hx, hy = hit
        up = getattr(event, "step", 0) > 0 or getattr(event, "button", "") == "up"
        ctrl, shift = _event_mods(event)
        x0, x1 = self._full_xlim
        if shift:
            lo, hi = self._view_x
            w = hi - lo
            step = 0.15 * w * (-1 if up else 1)
            lo2 = min(max(0.0, lo + step), 1.0 - w)
            self._view_x[:] = [lo2, lo2 + w]
            self._apply_zoom()
        elif ctrl:
            ylo, yhi = hax.get_ylim()
            frac_in_view = (hy - ylo) / max(1e-9, yhi - ylo)
            vlo, vhi = self._view_y
            self._wheel_zoom("y", 1 if up else -1,
                             anchor=vlo + frac_in_view * (vhi - vlo))
        else:
            hr = self._to_ref(self._ax_index(hax), hx)
            self._wheel_zoom("x", 1 if up else -1,
                             anchor=(hr - x0) / max(1e-9, x1 - x0))

    def _box_zoom_allowed(self, event):
        """Plain drag zooms in the viewer; while picking or measuring, the
        left button belongs to that tool, so box-zoom needs Ctrl."""
        ctrl, _shift = _event_mods(event)
        if getattr(self, "genotyping_active", False):
            return ctrl
        return True

    def _on_plot_press(self, event):
        hit = self._plot_hit(event) if self._full_xlim is not None else None
        if hit is None:
            return
        hax, hx, hy = hit
        ctrl, shift = _event_mods(event)
        idx = self._ax_index(hax)
        if event.button == 2:                           # middle drag: pan
            self._plot_drag = {"kind": "pan", "ax": hax,
                               "px": event.x, "view": list(self._view_x)}
        elif event.button == 1 and shift and not ctrl:  # slide this graph
            key = self._axes_keys[idx] if idx is not None and idx < len(self._axes_keys) else None
            if key is None:
                self.status_var.set("This view cannot be aligned (wrap rows).")
                return
            self._plot_drag = {"kind": "align", "ax": hax, "key": key, "idx": idx,
                               "px": event.x,
                               "base": {k: self._align.get(k, (1.0, 0.0))
                                        for k in self._group_keys(idx)},
                               "window": self._ref_window()}
        elif event.button == 1 and self._box_zoom_allowed(event):
            if getattr(event, "dblclick", False):
                self._reset_zoom()
                return
            from matplotlib.patches import Rectangle
            ax = hax
            rect = Rectangle((hx, hy), 0, 0, fill=True,
                             facecolor="#2F6FB5", alpha=0.12,
                             edgecolor="#2F6FB5", lw=0.8, zorder=20)
            ax.add_patch(rect)
            self._plot_drag = {"kind": "box", "ax": ax, "x": hx,
                               "y": hy, "px": event.x, "py": event.y,
                               "rect": rect}

    def _on_plot_motion(self, event):
        d = self._plot_drag
        if d is None:
            hit = self._plot_hit(event) if self._plot_axes else None
            self._hover_idx = self._ax_index(hit[0]) if hit else None
            return
        if d["kind"] == "align":
            lo, hi = d["window"]
            bbox = d["ax"].get_window_extent()
            shift_ref = (event.x - d["px"]) / max(1.0, bbox.width) * (hi - lo)
            self._shift_group(d["idx"], shift_ref, base=d["base"])
            self._update_ref_full(keep=d["window"])
            self._refresh_align_titles()
            self._apply_zoom()
            return
        if d["kind"] == "pan":
            x0, x1 = self._full_xlim
            ax = d["ax"]
            bbox = ax.get_window_extent()
            lo, hi = d["view"]
            w = hi - lo
            dfrac = (event.x - d["px"]) / max(1.0, bbox.width) * w
            lo2 = min(max(0.0, lo - dfrac), 1.0 - w)
            self._view_x[:] = [lo2, lo2 + w]
            self._apply_zoom()
            return
        ax = d["ax"]
        ex, ey = ax.transData.inverted().transform((event.x, event.y))
        xlo, xhi = ax.get_xlim()
        ylo, yhi = ax.get_ylim()
        ex = min(max(ex, xlo), xhi)
        ey = min(max(ey, ylo), yhi)
        x_only = abs(event.y - d["py"]) < 15
        r = d["rect"]
        r.set_x(min(d["x"], ex))
        r.set_width(abs(ex - d["x"]))
        if x_only:
            r.set_y(ylo)
            r.set_height(yhi - ylo)
        else:
            r.set_y(min(d["y"], ey))
            r.set_height(abs(ey - d["y"]))
        self.canvas.draw_idle()

    def _on_plot_release(self, event):
        d, self._plot_drag = self._plot_drag, None
        if d is not None and d["kind"] == "align":
            self._persist_alignment(list(d["base"]))
            a, b = self._align.get(d["key"], (1.0, 0.0))
            self.status_var.set(f"{Path(d['key']).parent.name}/{Path(d['key']).stem}: "
                                f"shifted {b:+.0f} scans  (Alt+Home resets this graph)")
            return
        if d is None or d["kind"] != "box":
            return
        r = d["rect"]
        xa, w = r.get_x(), r.get_width()
        ya, h = r.get_y(), r.get_height()
        try:
            r.remove()
        except (ValueError, AttributeError):
            pass
        if abs(event.x - d["px"]) < 6 or w <= 0:
            self.canvas.draw_idle()        # a click, not a drag
            return
        ylo, yhi = d["ax"].get_ylim()
        box_y = h > 0 and (ya > ylo + 1e-9 or ya + h < yhi - 1e-9)
        x0, x1 = self._full_xlim
        k_ax = self._ax_index(d["ax"])
        ra, rb = sorted((self._to_ref(k_ax, xa), self._to_ref(k_ax, xa + w)))
        a, b = max(x0, ra), min(x1, rb)
        if b - a < 2:
            self.canvas.draw_idle()
            return
        self._view_x[:] = [(a - x0) / (x1 - x0), (b - x0) / (x1 - x0)]
        if box_y:
            # The Y reference may itself follow the new X window ("Fit Y to
            # view"), so take the drawn band against the range it will have.
            try:
                k = self._plot_axes.index(d["ax"])
            except ValueError:
                k = 0
            ref_lo, ref_hi = self._y_references(a, b)[k]
            span = max(1e-9, ref_hi - ref_lo)
            self._view_y[:] = [min(max(0.0, (ya - ref_lo) / span), 1.0),
                               min(1.0, max(0.0, (ya + h - ref_lo) / span))]
            if self._view_y[1] - self._view_y[0] < 1e-6:
                self._view_y[:] = [0.0, 1.0]
        self._apply_zoom()

    # ---------------------------------------------------------- alignment
    @staticmethod
    def _akey(path):
        return str(Path(path).resolve())

    def _ax_index(self, ax):
        try:
            return self._plot_axes.index(ax)
        except ValueError:
            return None

    def _ab(self, i):
        keys = getattr(self, "_axes_keys", None) or []
        key = keys[i] if i is not None and 0 <= i < len(keys) else None
        return self._align.get(key, (1.0, 0.0)) if key else (1.0, 0.0)

    def _to_ref(self, i, s):
        a, b = self._ab(i)
        return a * s + b

    def _from_ref(self, i, r):
        a, b = self._ab(i)
        return (r - b) / a

    def _ref_window(self):
        x0, x1 = self._full_xlim
        f0, f1 = self._view_x
        return x0 + (x1 - x0) * f0, x0 + (x1 - x0) * f1

    def _update_ref_full(self, keep=None):
        """The shared X range is the union of every graph's aligned range;
        with *keep* (a reference window) the visible scans stay put."""
        if self._data_full is None:
            return
        g0, g1 = self._data_full
        ends = [self._to_ref(i, s) for i in range(len(self._plot_axes)) for s in (g0, g1)]
        self._full_xlim = (min(ends), max(ends)) if ends else (g0, g1)
        if keep is not None:
            x0, x1 = self._full_xlim
            lo, hi = keep
            span = max(1e-9, x1 - x0)
            f0 = min(max(0.0, (lo - x0) / span), 1.0)
            f1 = min(max(f0 + 1e-6, (hi - x0) / span), 1.0)
            self._view_x[:] = [f0, f1]

    def _set_alignment(self, key, ab, keep=None, redraw=False):
        a, b = ab
        if abs(a - 1.0) < 1e-9 and abs(b) < 1e-9:
            self._align.pop(key, None)
        else:
            self._align[key] = (float(a), float(b))
        if redraw:
            self.redraw()
            return
        self._update_ref_full(keep=keep if keep is not None else self._ref_window())
        self._refresh_align_titles()
        self._apply_zoom()

    def _restore_alignment(self, path):
        """Give a freshly loaded trace back the alignment saved for its well
        (library/<run>.json, key "align": [a, b]).  A session change wins."""
        key = self._akey(path)
        if key in self._align:
            return
        try:
            entry = (self._library_wells(Path(path).parent.name) or {}).get(Path(path).stem) or {}
            a, b = entry.get("align") or (1.0, 0.0)
            a, b = float(a), float(b)
        except Exception:                                  # noqa: BLE001
            return
        if a > 0 and (abs(a - 1.0) > 1e-9 or abs(b) > 1e-9):
            self._align[key] = (a, b)

    def _persist_alignment(self, keys):
        """Write the alignments of *keys* to the library -- batched, so a run
        of Alt-nudges or an auto-align is one write per run, not one per step."""
        pending = getattr(self, "_align_pending", None)
        if pending is None:
            pending = self._align_pending = set()
        pending.update(k for k in keys if k)
        job = getattr(self, "_align_save_job", None)
        if job is not None:
            try:
                self.after_cancel(job)
            except tk.TclError:
                pass
        try:
            self._align_save_job = self.after(500, self._flush_alignment)
        except tk.TclError:
            self._flush_alignment()

    def _flush_alignment(self):
        self._align_save_job = None
        pending, self._align_pending = getattr(self, "_align_pending", set()) or set(), set()
        by_run = {}
        for key in pending:
            by_run.setdefault(Path(key).parent.name, []).append(key)
        for run, keys in by_run.items():
            try:
                wells = self._library_wells(run)
                for key in keys:
                    well = Path(key).stem
                    entry = dict(wells.get(well) or {})
                    ab = self._align.get(key)
                    if ab is None:
                        entry.pop("align", None)
                    else:
                        entry["align"] = [round(ab[0], 8), round(ab[1], 4)]
                    if entry:
                        wells[well] = entry
                    else:
                        wells.pop(well, None)
                mark_library.save_run(run, wells)
            except Exception as exc:                      # noqa: BLE001
                print(f"library: could not save alignment for {run}: {exc}", file=sys.stderr)

    def _refresh_align_titles(self):
        for i, (artist, base) in enumerate(getattr(self, "_axes_titles", []) or []):
            a, b = self._ab(i)
            suffix = ""
            if abs(a - 1.0) > 1e-6 or abs(b) > 0.05:
                suffix = f"   ⇄ {b:+.0f} scans"
                if abs(a - 1.0) > 1e-6:
                    suffix += f" ×{a:.4f}"
            try:
                artist.set_text(base + suffix)
            except Exception:                              # noqa: BLE001
                pass

    def _group_keys(self, i):
        """Every file drawn in graph *i*: with Overlay runs, all its runs (so
        moving the graph moves them together); otherwise just its own file."""
        keys = self._axes_keys
        if i is None or i >= len(keys) or keys[i] is None:
            return []
        if self.overlay_runs.get() and not getattr(self, "genotyping_active", False):
            groups = self._view_groups(self._graph_count())
            if i < len(groups):
                return [self._akey(p) for p in groups[i]]
        return [keys[i]]

    def _shift_group(self, i, delta, base=None):
        """Move graph *i* (and its overlaid runs) by *delta* reference scans,
        from *base* ({key: (a, b)}, default: the current alignments)."""
        for k in self._group_keys(i):
            a, b = (base or {}).get(k, self._align.get(k, (1.0, 0.0)))
            if abs(a - 1.0) < 1e-9 and abs(b + delta) < 1e-9:
                self._align.pop(k, None)
            else:
                self._align[k] = (a, b + delta)

    def _nudge(self, scans):
        """Alt+←/→ (Alt+Shift: ×10): move the graph under the mouse."""
        i = self._hover_idx
        keys = self._axes_keys
        if i is None or i >= len(keys) or keys[i] is None or self._full_xlim is None:
            self.status_var.set("Point at a graph, then Alt+← / Alt+→ to move it.")
            return "break"
        win = self._ref_window()
        self._shift_group(i, scans)
        self._update_ref_full(keep=win)
        self._refresh_align_titles()
        self._apply_zoom()
        self._persist_alignment(self._group_keys(i))
        return "break"

    def _reset_alignment(self, which="hover"):
        if which == "all":
            keys = {k for i in range(len(self._plot_axes)) for k in self._group_keys(i)}
            for k in list(self._align):
                if k in keys:
                    self._align.pop(k)
            self._persist_alignment(keys)
            self.redraw()
            self.status_var.set("Alignment reset for the graphs on screen.")
            return "break"
        i = self._hover_idx
        if i is not None and i < len(self._axes_keys) and self._axes_keys[i]:
            keys = self._group_keys(i)
            for k in keys:
                self._align.pop(k, None)
            self._persist_alignment(keys)
            self.redraw()
        return "break"

    def _envelope(self, i, grid, own_runs_only=True, run=None):
        """Summed trace signal of graph *i* sampled at reference positions."""
        ax = self._plot_axes[i]
        a, b = self._ab(i)
        own = (grid - b) / a                      # this graph's own scans
        total = np.zeros_like(grid, dtype=float)
        for line in ax.get_lines():
            lab = str(line.get_label())
            if not lab.startswith("Ch"):
                continue
            if run is None and own_runs_only and "(" in lab:
                continue                         # an overlaid run, not this graph's own
            if run is not None and f"({run})" not in lab:
                continue
            xd = np.asarray(line.get_xdata(), float)
            yd = np.asarray(line.get_ydata(), float)
            if xd.size < 2:
                continue
            total += np.interp(own, xd, yd, left=np.nan, right=np.nan)
        return total

    @staticmethod
    def _best_lag(e0, e1, max_lag):
        """Samples to move *e1* right so its features sit on *e0*'s."""
        ok = np.isfinite(e0) & np.isfinite(e1)
        if ok.sum() < 10:
            return None
        a = np.where(ok, e0, 0.0)
        b = np.where(ok, e1, 0.0)
        a = (a - a[ok].mean()) / (a[ok].std() or 1.0)
        b = (b - b[ok].mean()) / (b[ok].std() or 1.0)
        a[~ok] = 0.0
        b[~ok] = 0.0
        n = a.size
        c = np.correlate(a, b, mode="full")
        lags = np.arange(-(n - 1), n)
        keep = np.abs(lags) <= max_lag
        return int(lags[keep][np.argmax(c[keep])])

    def align_to_top(self):
        """Line up what is in view by cross-correlation of the summed traces.

        With Overlay runs: first each overlaid run onto its graph's first run,
        then every graph (with all its runs) onto the top graph; without it,
        every graph onto the top one."""
        if len(self._plot_axes) < 1 or self._full_xlim is None:
            return
        lo, hi = self._ref_window()
        n = int(min(6000, max(50, hi - lo)))
        grid = np.linspace(lo, hi, n)
        step = (hi - lo) / max(1, n - 1)
        max_lag = int(0.3 * n)
        moved = []
        overlay = self.overlay_runs.get() and not getattr(self, "genotyping_active", False)
        if overlay:
            # 1) runs within each graph onto the graph's first run
            for i, group in enumerate(self._view_groups(self._graph_count())):
                if i >= len(self._plot_axes):
                    break
                own = self._envelope(i, grid)
                a0 = self._ab(i)[0]
                for other in group[1:]:
                    lag = self._best_lag(own, self._envelope(i, grid, run=other.parent.name), max_lag)
                    if lag:
                        k = self._akey(other)
                        a, b = self._align.get(k, (1.0, 0.0))
                        self._align[k] = (a, b + lag * step * a0)
                        moved.append(f"{other.parent.name[-5:]}/{other.stem} {lag * step:+.0f}")
        # 2) every graph (all its runs) onto the top one
        e0 = self._envelope(0, grid)
        for i in range(1, len(self._plot_axes)):
            if self._axes_keys[i] is None:
                continue
            lag = self._best_lag(e0, self._envelope(i, grid), max_lag)
            if lag:
                self._shift_group(i, lag * step)
                moved.append(f"{Path(self._axes_keys[i]).stem} {lag * step:+.0f}")
        keys = [k for i in range(len(self._plot_axes)) for k in self._group_keys(i)]
        self._persist_alignment(keys)
        self.redraw()
        self.status_var.set("Aligned on the scans in view: " + (", ".join(moved) or "already aligned"))

    def _is_scans(self, key):
        pk = self._gen_pickers.get(key) if hasattr(self, "_gen_pickers") else None
        if pk is not None and getattr(pk, "std", None):
            return sorted(float(s) for s, _n in pk.std)
        path = Path(key)
        try:
            entry = (self._library_wells(path.parent.name) or {}).get(path.stem) or {}
        except Exception:                                  # noqa: BLE001
            entry = {}
        pairs = entry.get("std") or []
        return sorted(float(p[0]) for p in pairs if p)

    def align_on_standard(self):
        """Line the graphs up on their marked internal standard: the bands of
        each well are mapped onto the top marked well's (offset, plus stretch
        when three or more bands agree on it)."""
        keys = [k for k in self._axes_keys if k]
        marks = {k: self._is_scans(k) for k in keys}
        ref = next((k for k in keys if len(marks[k]) >= 2), None)
        if ref is None:
            messagebox.showinfo("Align on internal standard",
                                "None of the graphs on screen has a marked internal standard.")
            return
        r = np.asarray(marks[ref])
        a0, b0 = self._align.get(ref, (1.0, 0.0))
        done, skipped = [], []
        for k in keys:
            if k == ref:
                continue
            s = np.asarray(marks[k])
            n = min(len(s), len(r))
            if n < 2 or len(s) != len(r):
                skipped.append(Path(k).stem)
                continue
            m, c = np.polyfit(s, r, 1) if n >= 3 else (1.0, float(np.median(r - s)))
            if not (0.9 < m < 1.1):
                m, c = 1.0, float(np.median(r - s))
            self._align[k] = (a0 * m, a0 * c + b0)
            done.append(Path(k).stem)
        self._persist_alignment(keys)
        self.redraw()
        msg = f"Aligned on the internal standard of {Path(ref).stem}: {', '.join(done) or 'nothing to move'}"
        if skipped:
            msg += f" · skipped (no matching IS marks): {', '.join(skipped)}"
        self.status_var.set(msg)

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
        """Push the X/Y view onto every visible graph + both bars."""
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
        fx0, fx1 = self._view_x
        fy0, fy1 = self._view_y
        w = (x1 - x0) * max(1e-9, fx1 - fx0)
        xlo = x0 + (x1 - x0) * fx0
        refs = self._y_references(xlo, xlo + w)
        for i, (ax, (y0, y1)) in enumerate(zip(self._plot_axes, refs)):
            h = (y1 - y0) * max(1e-9, fy1 - fy0)
            ax.set_xlim(self._from_ref(i, xlo), self._from_ref(i, xlo + w))
            ax.set_ylim(y0 + (y1 - y0) * fy0, y0 + (y1 - y0) * fy0 + h)
        if hasattr(self, "xbar"):
            self.xbar.set(fx0, fx1)
            self.ybar.set(fy0, fy1)
        self._show_scan_range(xlo, xlo + w)
        self.canvas.draw_idle()

    def _show_scan_range(self, a, b):
        """Mirror the visible scans in the range boxes (unless being typed in)."""
        if not hasattr(self, "scan_from"):
            return
        try:
            focus = self.focus_get()
        except (KeyError, tk.TclError):
            focus = None
        if focus in (getattr(self, "_scan_from_entry", None),
                     getattr(self, "_scan_to_entry", None)):
            return
        k = len(self._plot_axes) - 1
        self.scan_from.set(str(int(round(self._from_ref(k, a)))))
        self.scan_to.set(str(int(round(self._from_ref(k, b)))))

    def _fmt_time(self, scan, _pos=None):
        """Scan number -> minutes:seconds label (whole seconds)."""
        total = int(round(float(scan) / SCAN_RATE_HZ))
        return f"{total // 60}:{total % 60:02d}"

    def _style_x_axis(self, ax, is_bottom, label=True):
        """Show x tick labels on the bottom pane only; label scans or time."""
        ax.tick_params(axis="x", labelbottom=bool(is_bottom))
        if not is_bottom:
            return
        if self.x_time.get():
            ax.xaxis.set_major_formatter(FuncFormatter(self._fmt_time))
        if label:
            ax.set_xlabel("Time (min)" if self.x_time.get() else "Scan")

    def _add_volt_label(self):
        """One shared Y-axis label, centred on the whole stack with some air."""
        self.fig.supylabel("Volt", fontsize=8, x=0.014, color="#333")

    def redraw(self):
        """Repaint the graphs, then move the table to match them.

        Paging, list selection, the tour -- everything that changes which
        wells are on screen arrives here, so the table follows the figures
        rather than being scrolled separately and left behind.  One funnel is
        what makes the two read as one view instead of two lists.
        """
        self._redraw_plots()
        self._plate_map_refresh()
        order = self._display_order()
        if order != getattr(self, "_labelled_order", None):
            # the channel bar names each channel's base for the wells on
            # screen, which may be a run with its own dye order
            self._labelled_order = order
            self._refresh_channel_labels()
            self._refresh_channel_buttons()
        self._below_refresh()
        self._below_follow()

    def _redraw_plots(self):
        if self.genotyping_active:
            self._redraw_genotyping()
            return
        self.fig.clear()
        self._plot_axes = []
        self._axes_keys = []
        self._axes_titles = []
        self._data_full = None
        self._full_xlim = None
        self._full_ylim = None
        colors = self._theme_colors()
        n = self._graph_count()
        groups = self._view_groups(n)
        paths = [p for g in groups for p in g]
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
            self._wrap_active = False
            self._update_page_label([])
            self.canvas.draw_idle()
            self.seq_text.delete("1.0", tk.END)
            return

        settings = self._settings_from_ui()
        wrap = settings.view_mode == "wrap"

        self._wrap_active = bool(wrap)
        self._update_page_label(paths)
        if wrap:
            self._draw_wrap(paths, settings, colors)
            # wrap rows are slices of one well: nothing to align
            self._axes_keys = [None] * len(self._plot_axes)
            self._data_full = self._full_xlim
            self.fig.tight_layout(rect=(0.045, 0, 1, 1), h_pad=0.3, w_pad=0.2, pad=0.4)
            self._add_volt_label()
            self._apply_zoom()
            self.canvas.draw_idle()
            self._show_sequence()
            return

        gx0 = gx1 = gy0 = gy1 = None
        n_ax = len(groups)
        for i, group in enumerate(groups):
            ax = self.fig.add_subplot(n_ax, 1, i + 1)
            path = group[0]
            try:
                doc = self._ensure_doc(path)
            except Exception as e:
                ax.text(0.5, 0.5, f"Load error: {e}", ha="center", transform=ax.transAxes)
                continue
            for k, other in enumerate(group[1:], start=1):
                gy = self._overlay_run(ax, other, settings, RUN_STYLES[k % len(RUN_STYLES)],
                                       primary=self._akey(path))
                if gy is not None:
                    gy0 = gy[0] if gy0 is None else min(gy0, gy[0])
                    gy1 = gy[1] if gy1 is None else max(gy1, gy[1])
                    gx1 = gy[2] if gx1 is None else max(gx1, gy[2])
                    gx0 = 0 if gx0 is None else min(gx0, 0)
            dset = self._settings_for(doc, settings)
            tr = display_trace(doc, dset)
            s0 = max(0, settings.signal_start)
            s1 = settings.signal_end if settings.signal_end > 0 else tr.shape[0]
            s1 = min(tr.shape[0], s1)
            x = np.arange(s0, s1)
            plotted = False
            for ci, base, col, color in self._channel_lines(dset):
                if not self.chan_show[ci].get():
                    continue
                seg = tr[s0:s1, col]
                ax.plot(x, seg, color=color, lw=0.7, label=f"Ch{ci + 1} {base}")
                if seg.size:
                    m0 = float(np.nanmin(seg)); m1 = float(np.nanmax(seg))
                    gy0 = m0 if gy0 is None else min(gy0, m0)
                    gy1 = m1 if gy1 is None else max(gy1, m1)
                plotted = True
            if s1 > s0:
                gx0 = s0 if gx0 is None else min(gx0, s0)
                gx1 = s1 if gx1 is None else max(gx1, s1)
            self._plot_axes.append(ax)
            self._axes_keys.append(self._akey(path))
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

            ax.tick_params(labelsize=7)
            if len(group) > 1:
                # one plate per graph: name it once, then each run's own part
                plate = self._plate_of(path)
                names = "   ".join(
                    f"{'━' if k == 0 else {'--': '╍', ':': '┈', '-.': '┄'}.get(RUN_STYLES[k % len(RUN_STYLES)], '━')} "
                    f"{g.parent.name[len(plate):].lstrip('_- ') or g.parent.name}"
                    for k, g in enumerate(group))
                label = self._graph_title(path, f"{plate}    {names}")
            else:
                label = self._graph_title(path, doc.path.parent.name)
            self._axes_titles.append((ax.text(
                0.004, 0.995, label, transform=ax.transAxes, ha="left", va="top",
                fontsize=8, color="#222", zorder=6), label))
            self._draw_ctce_marks(ax, path)
            self._style_x_axis(ax, i == n_ax - 1)

            if settings.view_mode == "esd":
                self._draw_esd_peaks(
                    ax, doc, s0, s1, self._channel_lines(dset))
            elif doc.sequence and settings.view_mode == "called" and doc.peak_positions:
                if self.show_letters.get():
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
                    and settings.view_mode not in ("raw", "esd")):
                self._draw_quality(ax, doc, s0, s1)
        if gx0 is not None and gy0 is not None and gx1 > gx0:
            pad = 0.02 * (gy1 - gy0) or 1.0
            self._full_xlim = (float(gx0), float(gx1))
            self._full_ylim = (float(gy0 - pad), float(gy1 + pad))
        self._data_full = self._full_xlim
        self._update_ref_full()
        self._refresh_align_titles()
        self.fig.tight_layout(rect=(0.045, 0, 1, 1), h_pad=0.25, pad=0.3)
        self._add_volt_label()
        self._apply_zoom()
        self.canvas.draw_idle()
        self._show_sequence()

    def _well_label(self, path):
        m = WELL_RE.search(Path(path).stem)
        if not m:
            return Path(path).stem
        w = m.group(1).upper()
        return f"{w[0]}{int(w[1:]):02d}"

    @staticmethod
    def _plate_of(path):
        """The plate a run folder belongs to: its name without a trailing
        run number, so ..._270910Run01 and ..._270910Run02 are one plate."""
        return re.sub(r"[_\- ]?run\s*\d+$", "", Path(path).parent.name, flags=re.I)

    def _graph_title(self, path, rest):
        """Graph corner label: the well in bold, then *rest* (run name or the
        overlay's run legend).  Only a clean well id goes through mathtext."""
        well = self._well_label(path)
        if re.fullmatch(r"[A-Za-z0-9]+", well or ""):
            return f"$\\bf{{{well}}}$   {str(rest).replace('$', '')}"
        return f"{well}   {str(rest).replace('$', '')}"

    def _view_groups(self, n):
        """What the viewer draws: n graphs, one file each -- or, with
        "Overlay runs", n wells, each holding every selected run of it."""
        sel = list(self.selected or [])
        if not (self.overlay_runs.get() and not getattr(self, "genotyping_active", False)):
            return [[p] for p in sel[:n]]
        groups, index = [], {}
        for p in sel:
            key = (self._plate_of(p), self._well_label(p))
            if key in index:
                groups[index[key]].append(p)
            elif len(groups) < n:
                index[key] = len(groups)
                groups.append([p])
        return groups

    def _overlay_run(self, ax, path, settings, style, primary=None):
        """Another run of the same well on *ax*: channel traces only, in the
        run's own line style, placed through its own alignment into the
        primary run's scan coordinates.  Returns (ymin, ymax, scans) or None."""
        try:
            doc = self._ensure_doc(path)
        except Exception:                                  # noqa: BLE001
            return None
        dset = self._settings_for(doc, settings)
        tr = display_trace(doc, dset)
        s0 = max(0, settings.signal_start)
        s1 = settings.signal_end if settings.signal_end > 0 else tr.shape[0]
        s1 = min(tr.shape[0], s1)
        x = np.arange(s0, s1)
        if primary is not None:
            a_k, b_k = self._align.get(self._akey(path), (1.0, 0.0))
            a_0, b_0 = self._align.get(primary, (1.0, 0.0))
            x = (a_k * x + b_k - b_0) / a_0
        lo = hi = None
        for ci, base, col, color in self._channel_lines(dset):
            if not self.chan_show[ci].get():
                continue
            seg = tr[s0:s1, col]
            ax.plot(x, seg, color=color, lw=0.7, linestyle=style, alpha=0.85,
                    label=f"Ch{ci + 1} {base} ({path.parent.name})")
            if seg.size:
                lo = float(np.nanmin(seg)) if lo is None else min(lo, float(np.nanmin(seg)))
                hi = float(np.nanmax(seg)) if hi is None else max(hi, float(np.nanmax(seg)))
        if primary is not None:
            a_k, b_k = self._align.get(self._akey(path), (1.0, 0.0))
            a_0, b_0 = self._align.get(primary, (1.0, 0.0))
            self._draw_ctce_marks(ax, path, xmap=lambda s: (a_k * s + b_k - b_0) / a_0,
                                  run=path.parent.name)
        else:
            self._draw_ctce_marks(ax, path, run=path.parent.name)
        return None if lo is None else (lo, hi, s1)

    # ------------------------------------------------------------------ I/O
    # ------------------------------------------------------------- CTCE
    CTCE_STYLE = {"melting": ("v", "#C62828", 9), "satellite": ("v", "#EF8A17", 6),
                  "static": ("o", "#8A8F98", 5)}

    def _draw_ctce_marks(self, ax, path, xmap=None, run=None):
        """Mark the CTCE classes of *path*'s peaks on *ax*: peaks of interest
        (red, with their template label), satellites (orange), static (grey).
        *xmap* places an overlaid run's scans into the graph's coordinates."""
        marks = (getattr(self, "_ctce_marks", None) or {}).get(self._akey(path))
        if not marks:
            return
        lines = [ln for ln in ax.get_lines() if str(ln.get_label()).startswith("Ch")
                 and ((run is not None and f"({run})" in str(ln.get_label()))
                      or (run is None and "(" not in str(ln.get_label())))]
        ov = (getattr(self, "_ctce_overlay", None) or {}).get(self._akey(path))
        if ov and xmap is None:
            self._draw_ctce_overlay(ax, ov, lines)
        for raw, role, label in marks:
            style = self.CTCE_STYLE.get(role)
            if style is None:
                continue
            x = xmap(raw) if xmap else raw
            y = max((float(np.interp(x, ln.get_xdata(), ln.get_ydata())) for ln in lines),
                    default=0.0)
            m, col, size = style
            ax.plot([x], [y], marker=m, color=col, markersize=size, linestyle="none",
                    zorder=8, markerfacecolor=col if role != "static" else "none",
                    label="_ctce")
            if label:
                ax.annotate(label, (x, y), xytext=(0, 7), textcoords="offset points",
                            ha="center", fontsize=7, color=col, zorder=9)

    def _draw_ctce_overlay(self, ax, ov, lines):
        """What the species measurement took away and kept, on the sample
        channel's own trace: the hump under the species (blue, not counted),
        and each block's straight baseline (red) -- its area is what lies
        between that line and the trace."""
        ch = ov.get("channel")
        ln = next((ln for ln in lines if str(ln.get_label()).startswith(f"Ch{ch} ")), None)
        blocks = ov.get("blocks") or []
        if ln is None or not blocks:
            return
        xd = np.asarray(ln.get_xdata(), float)
        yd = np.asarray(ln.get_ydata(), float)
        yq, hump = np.asarray(ov["yq"], float), np.asarray(ov["hump"], float)
        lo = max(0, int(blocks[0][1]) - 60)
        hi = min(len(yq) - 1, int(blocks[-1][2]) + 60)
        x = np.arange(lo, hi + 1)
        disp = np.interp(x, xd, yd)
        q = yq[lo:hi + 1]
        # the display may carry its own baseline / scale: fit display = k*yq + floor
        dd, dq = np.diff(disp), np.diff(q)
        ok = (q[1:] < 0.9 * q.max()) & (q[:-1] < 0.9 * q.max())   # not the clipped tops
        k = float((dd[ok] * dq[ok]).sum() / max((dq[ok] ** 2).sum(), 1e-9)) if ok.any() else 1.0
        if not 0.5 < k < 2.0:
            k = 1.0
        floor = disp - k * q
        good = q < 0.9 * q.max()
        if good.sum() > 5:                                 # bridge clipped stretches
            floor = np.interp(x, x[good], floor[good])
        from scipy.ndimage import median_filter
        floor = median_filter(floor, size=9)
        top = floor + k * hump[lo:hi + 1]
        ax.fill_between(x, floor, top, color="#1E88E5", alpha=0.28, lw=0, zorder=3,
                        label="_ctce_hump")
        ax.plot(x, top, color="#1E88E5", lw=0.9, ls="--", zorder=4, label="_ctce_hump")
        ax.annotate("hump (not counted)", (x[0], top[0]), xytext=(2, 4), textcoords="offset points",
                    fontsize=7, color="#1E88E5", zorder=9)
        for name, b0, b1, a, b in blocks:
            b0, b1 = int(b0), int(b1)
            xs = np.arange(b0, b1 + 1)
            fl = np.interp(xs, x, floor)
            base = fl + k * np.linspace(a, b, len(xs))
            tr = np.interp(xs, xd, yd)
            ax.fill_between(xs, base, np.maximum(tr, base), color="#C62828", alpha=0.12,
                            lw=0, zorder=3, label="_ctce_block")
            ax.plot([b0, b1], [base[0], base[-1]], color="#C62828", lw=1.1, zorder=5,
                    label="_ctce_block")
            ax.plot([b0, b0], [base[0], tr[0]], color="#C62828", lw=0.6, zorder=5, label="_ctce_block")
            ax.plot([b1, b1], [base[-1], tr[-1]], color="#C62828", lw=0.6, zorder=5, label="_ctce_block")

    def _ctce_plates(self):
        """{plate: [run folders]} for the runs in the file list."""
        import ctce
        return ctce.group_runs(sorted({Path(p).parent for p in (self.files or [])}))

    def _ctce_picks(self, run_folders):
        """Peaks picked in Manual peak picking on this plate's runs: the
        template seeds, {trace path: [scans]}."""
        folders = {str(Path(f).resolve()) for f in run_folders}
        picks = {}
        for pk in (getattr(self, "_gen_pickers", None) or {}).values():
            for rec in getattr(pk, "records", []) or []:
                if rec.get("kind", "main") != "main":
                    continue
                f = Path(rec.get("file", ""))
                if str(f.parent.resolve()) in folders:
                    picks.setdefault(str(f), []).append(float(rec["scan"]))
        return picks

    def show_ctce(self):
        """Auto genotyping ▸ CTCE analysis…: find the peaks of interest of a
        plate from how they move between runs at different temperatures (or
        by shape at one temperature, or from picked peaks), and quantify."""
        import ctce
        plates = self._ctce_plates()
        if not plates:
            messagebox.showinfo("CTCE analysis", "Add the run folders of a plate first.")
            return
        win = tk.Toplevel(self)
        win.title("CTCE analysis")
        win.geometry("1180x620")
        self._ctce_win = win
        top = ttk.Frame(win, padding=8)
        top.pack(fill=tk.X)

        def label_for(plate):
            temps = []
            for f in plates[plate]:
                try:
                    rp = load_trace(next(Path(f).glob("*.rsd"))).run_params or {}
                    temps.append(str(rp.get("Temperature", "?")))
                except Exception:                          # noqa: BLE001
                    temps.append("?")
            return f"{plate}   ({len(plates[plate])} runs: {' / '.join(temps)} °C)"
        names = {label_for(k): k for k in sorted(plates)}
        ttk.Label(top, text="Plate").grid(row=0, column=0, sticky="w")
        plate_var = tk.StringVar(value=next(iter(names)))
        ttk.Combobox(top, textvariable=plate_var, values=list(names), state="readonly",
                     width=70).grid(row=0, column=1, columnspan=5, sticky="w", padx=6)
        ttk.Label(top, text="Sample channel").grid(row=1, column=0, sticky="w", pady=(6, 0))
        ch_var = tk.StringVar(value="Auto (strongest)")
        ttk.Combobox(top, textvariable=ch_var, state="readonly", width=16,
                     values=["Auto (strongest)", "1", "2", "3", "4"]).grid(row=1, column=1, sticky="w", padx=6, pady=(6, 0))
        ttk.Label(top, text="Internal standard").grid(row=1, column=2, sticky="w", pady=(6, 0))
        is_var = tk.StringVar(value="None")
        ttk.Combobox(top, textvariable=is_var, state="readonly", width=8,
                     values=["None", "1", "2", "3", "4"]).grid(row=1, column=3, sticky="w", padx=6, pady=(6, 0))
        merge_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(top, text="Add stutter to its parent peak", variable=merge_var).grid(
            row=1, column=4, sticky="w", padx=(12, 0), pady=(6, 0))
        tpl_var = tk.BooleanVar(value=True)
        tpl_chk = ttk.Checkbutton(top, text="Use my picked peaks as template", variable=tpl_var)
        tpl_chk.grid(row=2, column=1, columnspan=3, sticky="w", padx=6, pady=(4, 0))
        species_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(top, text="Measure Hom2 / Hom3 / heteroduplex blocks (plate pattern)",
                        variable=species_var).grid(row=3, column=1, columnspan=4, sticky="w",
                                                   padx=6, pady=(4, 0))
        info = ttk.Label(top, foreground="#555")
        info.grid(row=2, column=4, columnspan=2, sticky="w", padx=(12, 0))

        def refresh_info(*_):
            n = len(self._ctce_picks(plates[names[plate_var.get()]]))
            info.config(text=(f"{n} picked well(s) on this plate" if n else
                              "no picks — Manual genotyping ▸ Manual peak picking to guide it"))
            tpl_chk.state(["!disabled"] if n else ["disabled"])
        plate_var.trace_add("write", refresh_info)
        refresh_info()

        cols = ("well", "temperatures", "method", "regimes", "measured_at", "species",
                "fractions", "fraction_sd", "minor_fraction", "shifts", "flags")
        heads = {"well": "well", "temperatures": "°C", "method": "method", "regimes": "peaks per run",
                 "measured_at": "measured", "species": "species", "fractions": "fractions",
                 "fraction_sd": "± SD", "minor_fraction": "minor", "shifts": "shift (scans)",
                 "flags": "flags"}
        widths = {"well": 50, "temperatures": 80, "method": 70, "regimes": 150, "measured_at": 70,
                  "species": 150, "fractions": 170, "fraction_sd": 120, "minor_fraction": 60,
                  "shifts": 90, "flags": 360}
        mid = ttk.Frame(win)
        mid.pack(fill=tk.BOTH, expand=True, padx=8)
        tree = ttk.Treeview(mid, columns=cols, show="headings", height=18)
        for c in cols:
            tree.heading(c, text=heads[c])
            tree.column(c, width=widths[c], stretch=(c == "flags"), anchor="w")
        vs = ttk.Scrollbar(mid, orient=tk.VERTICAL, command=tree.yview)
        tree.configure(yscrollcommand=vs.set)
        tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vs.pack(side=tk.RIGHT, fill=tk.Y)
        bot = ttk.Frame(win, padding=8)
        bot.pack(fill=tk.X)
        status = ttk.Label(bot, text="Choose the plate and press Analyse.")
        status.pack(side=tk.LEFT)
        state = {"results": [], "folders": []}

        def analyse():
            plate = names[plate_var.get()]
            folders = plates[plate]
            channel = None if ch_var.get().startswith("Auto") else int(ch_var.get())
            is_ch = None if is_var.get() == "None" else int(is_var.get())
            template = None
            picks = self._ctce_picks(folders) if tpl_var.get() else {}
            if picks:
                template = ctce.template_from_picks(picks, channel)
            wells = ctce.plate_wells(folders)
            tree.delete(*tree.get_children())
            results = []
            self._busy, self._cancel = True, False
            try:
                for k, w in enumerate(wells, 1):
                    if self._cancel:
                        status.config(text=f"Cancelled after {k - 1} wells.")
                        break
                    status.config(text=f"Analysing {w} ({k}/{len(wells)})…  Esc cancels")
                    win.update()
                    try:
                        results.append(ctce.analyse_well(
                            [Path(f) / f"{w}.rsd" for f in folders], channel, is_ch,
                            merge_var.get(), template))
                    except Exception as e:                 # noqa: BLE001
                        print(f"ctce: {w}: {e}", file=sys.stderr)
            finally:
                self._busy = False
            if species_var.get() and results:
                status.config(text="Learning the plate's species pattern…")
                win.update()
                ctce.quantify_species(results)
            ctce.flag_no_product(results)
            state["results"], state["folders"] = results, folders
            marks = {}
            for r in results:
                for run in r.runs:
                    marks[self._akey(run.path)] = [(p.raw_pos, p.role, p.label) for p in run.peaks]
            self._ctce_marks = marks
            self._ctce_overlay = {self._akey(run.path): run._species_overlay
                                  for r in results for run in r.runs
                                  if getattr(run, "_species_overlay", None)}
            for r in results:
                s = r.summary()
                tree.insert("", tk.END, iid=r.well, values=[s[c] for c in cols])
            how = {"template": "your picked peaks", "shape": "peak shape (one temperature)",
                   "change": "the change between temperatures",
                   "species": "the plate's Hom2/Hom3/het pattern (blocks with stutter)"}
            methods = {r.method for r in results}
            status.config(text=f"{len(results)} wells — peaks of interest from "
                               + ", ".join(how.get(m, m) for m in sorted(methods))
                               + ".  Double-click a well to show it.")
            self.redraw()

        def show_well(_e=None):
            sel = tree.selection()
            if not sel:
                return
            w = sel[0]
            files = [Path(f) / f"{w}.rsd" for f in state["folders"]]
            # show the runs the way the engine lined them up (session only;
            # the alignments saved in the library are not touched)
            res = next((r for r in state["results"] if r.well == w), None)
            if res is not None:
                ref = res.runs[0]
                for run in res.runs:
                    a = run.scale / ref.scale
                    b = (run.offset - ref.offset) / ref.scale
                    key = self._akey(run.path)
                    if abs(a - 1.0) < 1e-9 and abs(b) < 1e-9:
                        self._align.pop(key, None)
                    else:
                        self._align[key] = (a, b)
            self.overlay_runs.set(True)
            self.selected = files
            idx = {str(p.resolve()): i for i, p in enumerate(self.files)}
            self.file_list.selection_clear(0, tk.END)
            for f in files:
                i = idx.get(str(f.resolve()))
                if i is not None:
                    self.file_list.selection_set(i)
                    self.file_list.see(i)
            self.redraw()

        def save():
            if not state["results"]:
                return
            path = filedialog.asksaveasfilename(
                parent=win, title="Save CTCE table", defaultextension=".csv",
                filetypes=[("CSV", "*.csv"), ("All files", "*.*")])
            if not path:
                return
            import csv
            rows = [r.summary() for r in state["results"]]
            try:
                with open(path, "w", newline="", encoding="utf-8-sig") as fh:
                    wtr = csv.DictWriter(fh, fieldnames=list(rows[0]))
                    wtr.writeheader()
                    wtr.writerows(rows)
                status.config(text=f"Saved {len(rows)} wells to {Path(path).name}")
            except OSError as e:
                messagebox.showerror("Save CTCE table", str(e), parent=win)

        tree.bind("<Double-1>", show_well)
        tree.bind("<Return>", show_well)
        ttk.Button(bot, text="Close", command=win.destroy).pack(side=tk.RIGHT)
        ttk.Button(bot, text="Save table…", command=save).pack(side=tk.RIGHT, padx=4)
        ttk.Button(bot, text="Show well", command=show_well).pack(side=tk.RIGHT, padx=4)
        ttk.Button(bot, text="Analyse", command=analyse).pack(side=tk.RIGHT, padx=4)
        self._ctce_ui = SimpleNamespace(win=win, analyse=analyse, show_well=show_well, tree=tree,
                                        plate_var=plate_var, names=names, state=state,
                                        status=status, tpl_var=tpl_var)

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
        self.basecaller.set(data.get("basecaller", "mb1000_accuracy"))
        self.base_order_var.set(data.get("base_order", "ACTG"))
        self.bonus.set(data.get("channel_peak_bonus", 0.7))
        self.pullback.set(data.get("pullback_weight", 0.008))
        self._refresh_channel_buttons()
        self.status_var.set(f"Loaded settings from {path}")

    def show_help(self):
        """The user manual: contents on the left, the text on the right, and a
        search box (Ctrl+F, Enter = next match).  Text lives in user_manual.py."""
        import re as _re
        import user_manual as um
        old = getattr(self, "_help_win", None)
        try:
            if old is not None and old.winfo_exists():
                old.deiconify()
                old.lift()
                return
        except tk.TclError:
            pass
        win = tk.Toplevel(self)
        self._help_win = win
        win.title(um.TITLE)
        win.geometry("980x700")
        win.minsize(640, 420)
        win.transient(self)

        top = ttk.Frame(win)
        top.pack(fill=tk.X, padx=8, pady=(8, 4))
        ttk.Label(top, text="Search").pack(side=tk.LEFT)
        query = tk.StringVar()
        entry = ttk.Entry(top, textvariable=query, width=30)
        entry.pack(side=tk.LEFT, padx=(4, 4))
        found = ttk.Label(top, text="", foreground="#666")
        found.pack(side=tk.LEFT, padx=6)
        ttk.Button(top, text="Close", command=win.destroy).pack(side=tk.RIGHT)

        body = ttk.Panedwindow(win, orient=tk.HORIZONTAL)
        body.pack(fill=tk.BOTH, expand=True, padx=8, pady=(0, 8))
        toc = tk.Listbox(body, width=30, activestyle="none", exportselection=False,
                         font=("DejaVu Sans", 10), bg="#F4F6F9",
                         highlightthickness=0, borderwidth=0)
        body.add(toc, weight=0)
        box = scrolledtext.ScrolledText(body, wrap=tk.WORD, padx=16, pady=10,
                                        font=("DejaVu Sans", 10), bg="#FBFBF4",
                                        spacing1=1, spacing3=2)
        body.add(box, weight=1)
        box.tag_configure("h1", font=("DejaVu Sans", 15, "bold"), foreground="#1E4E8C",
                          spacing1=14, spacing3=6)
        box.tag_configure("title", font=("DejaVu Sans", 17, "bold"), spacing3=8)
        box.tag_configure("b", font=("DejaVu Sans", 10, "bold"))
        box.tag_configure("pre", font=("DejaVu Sans Mono", 9), wrap=tk.NONE,
                          background="#F1F1E6", lmargin1=12, lmargin2=12)
        box.tag_configure("bullet", lmargin1=8, lmargin2=24)
        box.tag_configure("hit", background="#FFE36E")
        box.tag_configure("cur", background="#FFB347")

        def add_inline(line, tags=()):
            parts = _re.split(r"(\*\*[^*]+\*\*)", line)
            for part in parts:
                if part.startswith("**") and part.endswith("**") and len(part) > 4:
                    box.insert(tk.END, part[2:-2], tags + ("b",))
                elif part:
                    box.insert(tk.END, part, tags)

        box.insert(tk.END, um.TITLE + "\n", ("title",))
        starts = []
        for i, (title, text) in enumerate(um.SECTIONS):
            toc.insert(tk.END, f"  {title}")
            starts.append(box.index("end-1c"))
            box.insert(tk.END, title + "\n", ("h1",))
            pre = False
            for line in text.split("\n"):
                if line.strip() == "```":
                    pre = not pre
                    continue
                if pre:
                    box.insert(tk.END, line + "\n", ("pre",))
                elif line.startswith("• "):
                    add_inline(line, ("bullet",))
                    box.insert(tk.END, "\n", ("bullet",))
                else:
                    add_inline(line)
                    box.insert(tk.END, "\n")
        box.insert(tk.END, "\n" * 40)        # so the last sections can scroll to the top
        box.configure(state=tk.DISABLED)

        def go(_e=None):
            sel = toc.curselection()
            if sel:
                box.yview(starts[sel[0]])
        toc.bind("<<ListboxSelect>>", go)

        state = {"hits": [], "i": -1, "q": ""}

        def search(_e=None):
            q = query.get().strip()
            if q != state["q"]:
                box.tag_remove("hit", "1.0", tk.END)
                state.update(hits=[], i=-1, q=q)
                if q:
                    idx = "1.0"
                    while True:
                        idx = box.search(q, idx, nocase=True, stopindex=tk.END)
                        if not idx:
                            break
                        end = f"{idx}+{len(q)}c"
                        box.tag_add("hit", idx, end)
                        state["hits"].append((idx, end))
                        idx = end
            box.tag_remove("cur", "1.0", tk.END)
            hits = state["hits"]
            if not q:
                found.config(text="")
                return "break"
            if not hits:
                found.config(text="not found")
                return "break"
            state["i"] = (state["i"] + 1) % len(hits)
            a_, b_ = hits[state["i"]]
            box.tag_add("cur", a_, b_)
            box.see(a_)
            found.config(text=f"{state['i'] + 1} of {len(hits)}")
            return "break"
        entry.bind("<Return>", search)
        entry.bind("<KP_Enter>", search)
        entry.bind("<KeyRelease>", lambda e: None if e.keysym in ("Return", "KP_Enter")
                   else (state.update(q=None), search()))
        for w in (win, box, toc):
            w.bind("<Control-f>", lambda e: (entry.focus_set(), "break")[1])
            w.bind("<Control-F>", lambda e: (entry.focus_set(), "break")[1])
        win.bind("<Escape>", lambda e: win.destroy())
        toc.selection_set(0)
        box.focus_set()
        self._help_parts = SimpleNamespace(box=box, toc=toc, entry=entry, query=query,
                                           found=found, search=search, go=go,
                                           starts=starts)

    def show_about(self):
        messagebox.showinfo(
            "About",
            "Limoncello CE Analyzer  —  CE trace viewer & base caller\n\n"
            "Serves all instruments that produce CE traces:\n"
            "Formats: .rsd, .scf, ABI .ab1, text/CSV traces\n"
            "Current trace: µA (RSD raw ÷10); fluorescence in Volts.\n"
            "Basecallers: mb1000/mb4000 × accuracy/length (pos_bonus07 and\n"
             "            pos_profile are legacy aliases of those; hz_soften,\n"
             "            raw_peaks remain), see the Basecall settings dialog.\n"
            "Multi-folder load, multi-graph, auto-tour, well sort/filter.\n"
            "ESD-style processed view under Base calling; manual peak picking\n"
            "genotyping (click peaks, +A tagging, sizes via internal standard,\n"
            "Excel/CSV/JSON export) under Genotyping.\n"
            "Uses our own tuned spacing tracker.",
        )


def main():
    ap = argparse.ArgumentParser(description="Limoncello CE Analyzer")
    ap.add_argument("--folder", action="append", type=Path, help="Data folder (repeatable)")
    ap.add_argument("--check", action="store_true",
                    help="report which dependencies are available and exit")
    args = ap.parse_args()
    if args.check:
        from analyzer_core import environment_report
        raise SystemExit(environment_report())
    app = LimoncelloAnalyzerApp(initial_folders=args.folder)
    app.mainloop()


if __name__ == "__main__":
    main()