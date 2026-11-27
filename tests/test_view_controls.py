"""View bar and plot mouse controls: graph count, Y scaling, zoom gestures,
overlay of runs, and the lower pane.  Built from the bundled M13 example
traces copied into two run folders of one plate, so it needs no external data;
skips where there is no display."""
import os
import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("tkinter")
pytest.importorskip("matplotlib")
if not os.environ.get("DISPLAY"):
    pytest.skip("no DISPLAY for the GUI test", allow_module_level=True)

import tkinter as tk  # noqa: E402
from matplotlib.backend_bases import MouseEvent  # noqa: E402

import sequence_analyzer as sa  # noqa: E402

M13 = Path(__file__).resolve().parents[1] / "example_data" / "M13"
WELLS = ("A01", "A02", "A03")


@pytest.fixture()
def app(tmp_path, monkeypatch):
    if not all((M13 / f"{w}.rsd").exists() for w in WELLS):
        pytest.skip("example M13 traces not present")
    import mark_library
    monkeypatch.setattr(mark_library, "LIBRARY_DIR", tmp_path / "library")
    root = tmp_path / "data"
    for run in ("PLATE_X_Run01", "PLATE_X_Run02"):
        (root / run).mkdir(parents=True)
        for w in WELLS:
            shutil.copy(M13 / f"{w}.rsd", root / run / f"{w}.rsd")
    try:
        a = sa.LimoncelloAnalyzerApp(initial_folders=[root])
    except tk.TclError as e:
        pytest.skip(f"cannot open a window: {e}")
    a.geometry("1400x900+0+0")
    a.overlay_runs.set(False)
    a.y_scale.set("shared")
    a.fit_y.set(True)
    a.always_table.set(False)
    a.update()
    yield a
    a.destroy()


def _select(app, idxs, n):
    app.n_graphs.set(n)
    app.file_list.selection_clear(0, tk.END)
    for i in idxs:
        app.file_list.selection_set(i)
    app.on_file_select()
    app.update()
    # let a pending margin refit (deferred after resizes) land now, so it
    # cannot move the graphs in the middle of a simulated drag
    job = getattr(app, "_relayout_job", None)
    if job is not None:
        app.after_cancel(job)
        app._relayout()
    app.update()


def _event(app, kind, ax, xd, yd=None, **kw):
    """A real matplotlib mouse event at scan *xd*; *yd* None = mid-height."""
    if yd is None:
        yd = sum(ax.get_ylim()) / 2.0
    x, y = ax.transData.transform((xd, yd))
    e = MouseEvent(kind, app.canvas, x, y, **kw)
    app.canvas.callbacks.process(kind, e)
    app.update()


def test_graph_count_is_clamped_and_survives_a_blank_box(app):
    app.n_graphs.set(99)
    assert app._graph_count() == sa.MAX_GRAPHS
    app.n_graphs.set(0)
    assert app._graph_count() == 1
    app.tk.globalsetvar(str(app.n_graphs), "")      # box being edited
    assert app._graph_count() == 1


def test_run_folders_of_one_plate_share_a_plate_name():
    f = sa.LimoncelloAnalyzerApp._plate_of
    assert f(Path("x/OY_rs1695_T9_270910Run01/A01.rsd")) == f(
        Path("x/OY_rs1695_T9_270910Run02/A01.rsd")) == "OY_rs1695_T9_270910"
    assert f(Path("x/M13/A01.rsd")) == "M13"


def test_sorting_by_well_puts_the_runs_of_a_well_together(app):
    names = [f"{p.parent.name}/{p.stem}" for p in app.files[:4]]
    assert names == ["PLATE_X_Run01/A01", "PLATE_X_Run02/A01",
                     "PLATE_X_Run01/A02", "PLATE_X_Run02/A02"]


def test_overlay_draws_one_graph_per_well_with_every_run(app):
    app.overlay_runs.set(True)
    _select(app, range(4), 2)
    assert len(app._plot_axes) == 2
    for ax in app._plot_axes:
        styles = {line.get_linestyle() for line in ax.get_lines()
                  if str(line.get_label()).startswith("Ch")}
        assert {"-", "--"} <= styles             # Run01 solid, Run02 dashed


def test_overlay_pages_by_well(app):
    app.overlay_runs.set(True)
    _select(app, [0, 1], 1)
    app._page_by(1)
    assert {p.stem for p in app.selected} == {"A02"}
    assert len(app.selected) == 2                # both runs of A02
    app._page_by(-1)
    assert {p.stem for p in app.selected} == {"A01"}


def test_more_selected_than_graphs_is_said_not_hidden(app):
    _select(app, range(4), 2)
    assert "showing 2 of 4" in app._page_lbl.cget("text")


def test_the_file_list_highlights_every_file_on_screen(app):
    """More graphs than files selected: the next files fill the graphs, and
    the list highlights them all, counted from the first selected down."""
    _select(app, [1], 3)
    assert list(app.file_list.curselection()) == [1, 2, 3]
    assert [str(p) for p in app.selected] == [str(app.files[i]) for i in (1, 2, 3)]
    assert len(app._plot_axes) == 3


def test_graphs_raised_after_selecting_extends_the_highlight(app):
    _select(app, [0], 1)
    assert list(app.file_list.curselection()) == [0]
    app.n_graphs.set(4)
    app._graphs_changed()
    app.update()
    assert list(app.file_list.curselection()) == [0, 1, 2, 3]


def test_last_files_fill_only_as_far_as_the_list_goes(app):
    n = len(app.files)
    _select(app, [n - 1], 3)
    assert list(app.file_list.curselection()) == [n - 1]


def test_typing_a_graph_count_extends_the_selection(app):
    """A number typed into the box (no Enter) takes effect after a pause."""
    import time
    _select(app, [0], 1)
    app.n_graphs.set(3)                          # as typing does
    t = time.time()
    while time.time() - t < 1.5:
        app.update()
        time.sleep(0.05)
    assert list(app.file_list.curselection()) == [0, 1, 2]
    assert len(app._plot_axes) == 3


def test_overlay_runs_fills_with_the_next_wells(app):
    app.overlay_runs.set(True)
    _select(app, [0], 1)
    first = app.files[0]
    app.n_graphs.set(2)
    app._graphs_changed()
    app.update()
    assert len(app._plot_axes) == 2
    wells = {p.stem for p in app.selected}
    assert first.stem in wells and len(wells) == 2
    nxt = (wells - {first.stem}).pop()
    runs = {p.parent.name for p in app.selected if p.stem == nxt}
    assert runs == {"PLATE_X_Run01", "PLATE_X_Run02"}    # every run of the added well
    sel = {str(app.files[i]) for i in app.file_list.curselection()}
    assert sel == {str(p) for p in app.selected}


def _click_pick(app, idxs):
    """Ctrl-click each of *idxs* in the file list, as the user does."""
    app.file_list.selection_clear(0, tk.END)
    for i in idxs:
        app._list_click = "pick"
        app.file_list.selection_set(i)
        app.on_file_select()
    app.update()


def test_ctrl_click_shows_exactly_the_wells_picked(app):
    """A01, A03 of Run01 and A02 of Run02, not in running order: exactly
    those, one graph each, nothing filled in between."""
    app.n_graphs.set(4)
    app._graphs_changed()
    _click_pick(app, [0, 2, 4])
    assert sorted(app.file_list.curselection()) == [0, 2, 4]
    assert [str(p) for p in app.selected] == [str(app.files[i]) for i in (0, 2, 4)]
    assert app._graph_count() == 3 and len(app._plot_axes) == 3


def test_plain_click_after_a_pick_fills_again(app):
    _click_pick(app, [0, 2])
    app.file_list.selection_clear(0, tk.END)
    app.file_list.selection_set(1)
    app._list_click = "plain"
    app.on_file_select()
    app.update()
    assert list(app.file_list.curselection()) == [1, 2]      # Graphs is 2 after the pick


def test_wheel_on_the_graph_zooms_x_about_the_cursor(app):
    _select(app, [0], 1)
    ax = app._plot_axes[0]
    lo, hi = ax.get_xlim()
    at = lo + 0.3 * (hi - lo)
    frac_before = (at - lo) / (hi - lo)
    _event(app, "scroll_event", ax, at, step=1)
    lo2, hi2 = ax.get_xlim()
    assert hi2 - lo2 < hi - lo                    # zoomed in
    # the scan under the cursor stays put (mouse positions are whole pixels)
    assert abs((at - lo2) / (hi2 - lo2) - frac_before) < 1e-3


def test_drag_zooms_to_the_dragged_scans(app):
    _select(app, [0], 1)
    ax = app._plot_axes[0]
    lo0, hi0 = ax.get_xlim()
    px = (hi0 - lo0) / ax.get_window_extent().width      # scans per pixel
    _event(app, "button_press_event", ax, 3000, button=1)
    _event(app, "motion_notify_event", ax, 4000)
    _event(app, "button_release_event", ax, 4000, button=1)
    lo, hi = ax.get_xlim()
    assert abs(lo - 3000) <= 2 * px and abs(hi - 4000) <= 2 * px
    assert abs(int(app.scan_from.get()) - lo) <= 1
    assert abs(int(app.scan_to.get()) - hi) <= 1


def test_scan_range_box_sets_the_view(app):
    _select(app, [0], 1)
    app.scan_from.set("2500")
    app.scan_to.set("2600")
    app._scan_entry_go()
    app.update()
    assert tuple(round(v) for v in app._plot_axes[0].get_xlim()) == (2500, 2600)


def test_each_graph_scale_differs_from_shared(app):
    _select(app, [0, 2], 2)                       # Run01/A01 and Run01/A02
    app.fit_y.set(False)
    app._apply_zoom()
    shared = [ax.get_ylim() for ax in app._plot_axes]
    assert shared[0] == shared[1]
    app.y_scale.set("each")
    app._apply_zoom()
    each = [ax.get_ylim() for ax in app._plot_axes]
    assert each[0] != each[1]


def test_fit_y_follows_the_scans_on_screen(app):
    _select(app, [0], 1)
    ax = app._plot_axes[0]
    full_top = ax.get_ylim()[1]
    app._zoom_to_scans(500, 900)                  # pre-signal baseline only
    assert ax.get_ylim()[1] < 0.5 * full_top


def test_lower_pane_hides_when_it_has_nothing_to_show(app):
    _select(app, [0], 1)
    assert app._below_shown is False
    app.always_table.set(True)
    app._below_refresh()
    assert app._below_shown is True


def test_ctrl_is_read_from_the_key_or_the_tk_state():
    assert sa._event_mods(SimpleNamespace(key="control", guiEvent=None)) == (True, False)
    assert sa._event_mods(SimpleNamespace(key=None,
                                          guiEvent=SimpleNamespace(state=0x0005))) == (True, True)
    assert sa._event_mods(SimpleNamespace(key=None, guiEvent=None)) == (False, False)


# --------------------------------------------------------------- dye order
def test_a_sequencing_file_declares_its_dye_order():
    import analyzer_core
    doc = analyzer_core.load_trace(M13 / "A01.rsd")
    assert analyzer_core.declared_dye_order(doc.footer_fields) == "TGCA"
    assert analyzer_core.declared_dye_order({"CHANNEL BASES": ["", "", "", ""]}) == ""
    assert analyzer_core.declared_dye_order({}) == ""


def test_declared_order_beats_the_general_setting_and_a_run_choice_beats_both(app):
    app.base_order_var.set("ACTG")
    _select(app, [0], 1)
    doc = app.docs[str(app.selected[0].resolve())]
    assert (doc.dye_order, doc.dye_order_source) == ("TGCA", "file")
    assert [cb.cget("text").strip() for cb in app._chan_cbs] == [
        "Ch1 T", "Ch2 G", "Ch3 C", "Ch4 A"]
    app._run_orders()["PLATE_X_Run01"] = "ACTG"          # chosen for this run
    app.docs.clear()
    _select(app, [0], 1)
    doc = app.docs[str(app.selected[0].resolve())]
    assert (doc.dye_order, doc.dye_order_source) == ("ACTG", "run")


# ------------------------------------------------------------- plate map
def test_well_and_call_parsing():
    rc = sa.LimoncelloAnalyzerApp._well_rc
    assert rc("A01") == ("A", 1) and rc("p24") == ("P", 24) and rc("Run01") is None
    cls = sa.LimoncelloAnalyzerApp._call_class
    assert [cls(c) for c in ("hom-1", "hom-2", "het", "no-call", "")] == [
        "hom1", "hom2", "het", "nocall", "unscored"]


def test_plate_map_click_shows_a_well_and_a_row(app):
    app.show_plate_map()
    app.update()
    app._plate_run.set("PLATE_X_Run02")
    app._plate_map_refresh(force=True)
    rows, ncol, cell, x0, y0 = app._plate_layout
    assert (rows, ncol) == ("ABCDEFGH", 12)
    ev = SimpleNamespace(x=int(x0 + 1.5 * cell), y=int(y0 + 0.5 * cell))   # A02
    app._plate_click(ev)
    assert [f"{p.parent.name}/{p.stem}" for p in app.selected] == ["PLATE_X_Run02/A02"]
    ev = SimpleNamespace(x=int(x0 / 2), y=int(y0 + 0.5 * cell))            # row A
    app._plate_click(ev)
    assert sorted(p.stem for p in app.selected) == ["A01", "A02", "A03"]
    assert app._graph_count() >= 3


# ------------------------------------------------------------ mode banner
def test_picking_shows_a_banner_and_esc_leaves(app):
    _select(app, [0], 1)
    app.toggle_genotyping_picking()
    app.update()
    assert app._mode_banner.winfo_ismapped()
    assert "picking" in app._mode_banner_lbl.cget("text").lower()
    app._on_escape()
    app.update()
    assert not app.genotyping_active
    assert not app._mode_banner.winfo_ismapped()


def test_scroll_bars_pan_a_zoomed_view(app):
    """The scroll bars are wired to _bar_cmd; dragging the thumb sends
    ("moveto", f) and the arrows/trough send ("scroll", n, units|pages)."""
    _select(app, [0], 1)
    ax = app._plot_axes[0]
    app._zoom_to_scans(1000, 1500)
    w0 = ax.get_xlim()[1] - ax.get_xlim()[0]
    app.xbar.cget("command")                      # the bar has a command
    app.tk.call(app.xbar.cget("command"), "moveto", "0.5")   # as Tk calls it
    app.update()
    lo, hi = ax.get_xlim()
    assert abs((hi - lo) - w0) < 1e-6             # same zoom width …
    assert lo > 1500                              # … slid to the middle
    app.tk.call(app.xbar.cget("command"), "scroll", "-1", "pages")
    app.update()
    assert ax.get_xlim()[0] < lo                  # paged back
    lo_y = ax.get_ylim()
    app._wheel_zoom("y", 1)
    app.tk.call(app.ybar.cget("command"), "moveto", "0.2")
    app.update()
    assert ax.get_ylim() != lo_y                  # the Y bar pans too


# ------------------------------------------------------------ alignment
def _two_runs_of_a01(app):
    _select(app, [0, 1], 2)            # Run01/A01 and Run02/A01: identical traces
    assert [p.parent.name for p in app.selected] == ["PLATE_X_Run01", "PLATE_X_Run02"]
    return app._plot_axes[0], app._plot_axes[1]


def test_aligning_a_graph_moves_only_its_own_view(app):
    ax0, ax1 = _two_runs_of_a01(app)
    app._zoom_to_scans(3000, 3500)
    key1 = app._axes_keys[1]
    app._set_alignment(key1, (1.0, 40.0))
    lo0, hi0 = ax0.get_xlim()
    lo1, hi1 = ax1.get_xlim()
    assert abs((lo0 - lo1) - 40) < 1e-6 and abs((hi0 - hi1) - 40) < 1e-6
    assert abs((hi0 - lo0) - 500) < 1e-6           # the zoom width is unchanged
    assert "+40 scans" in app._axes_titles[1][0].get_text()
    # the scan boxes show the bottom graph's own scans (its tick labels)
    assert (int(app.scan_from.get()), int(app.scan_to.get())) == (round(lo1), round(hi1))


def test_shift_drag_slides_one_graph(app):
    ax0, ax1 = _two_runs_of_a01(app)
    app._zoom_to_scans(3000, 3500)
    px = (3500 - 3000) / ax1.get_window_extent().width       # scans per pixel
    _event(app, "button_press_event", ax1, 3200, button=1, key="shift")
    _event(app, "motion_notify_event", ax1, 3260, key="shift")
    _event(app, "button_release_event", ax1, 3260, button=1, key="shift")
    a, b = app._align[app._axes_keys[1]]
    assert a == 1.0 and abs(b - 60) <= 3 * px                 # dragged ~60 scans right
    assert app._axes_keys[0] not in app._align               # the other graph stayed


def test_alt_nudge_moves_the_graph_under_the_mouse(app):
    _two_runs_of_a01(app)
    app._hover_idx = 1
    app._nudge(10)
    app._nudge(-1)
    assert app._align[app._axes_keys[1]] == (1.0, 9.0)
    app._reset_alignment("hover")
    assert app._axes_keys[1] not in app._align


def test_align_to_top_undoes_a_known_shift(app):
    _two_runs_of_a01(app)
    app._zoom_to_scans(3000, 4000)
    app._set_alignment(app._axes_keys[1], (1.0, 37.0))       # misalign the copy
    app.align_to_top()
    a, b = app._align.get(app._axes_keys[1], (1.0, 0.0))
    assert abs(b) <= 1.5                                      # back on top of the original


def test_align_on_standard_fits_offset_and_stretch(app):
    _two_runs_of_a01(app)
    k0, k1 = app._axes_keys
    r = [2100.0, 2200.0, 2400.0, 2500.0]
    s = [x * 0.99 - 12 for x in r]                            # the other run: shifted and squeezed
    app._gen_pickers[k0] = SimpleNamespace(std=[(x, "") for x in r])
    app._gen_pickers[k1] = SimpleNamespace(std=[(x, "") for x in s])
    app.align_on_standard()
    a, b = app._align[k1]
    for x_s, x_r in zip(s, r):
        assert abs(a * x_s + b - x_r) < 1e-6                  # each band lands on its partner


def test_shift_click_does_not_pick(app):
    _select(app, [0], 1)
    app.toggle_genotyping_picking()
    app.update()
    ax = app._plot_axes[0]
    before = sum(len(pk.records) for pk in app._gen_pickers.values())
    _event(app, "button_press_event", ax, 3000, button=1, key="shift")
    _event(app, "button_release_event", ax, 3000, button=1, key="shift")
    assert sum(len(pk.records) for pk in app._gen_pickers.values()) == before


def test_alignment_is_saved_in_the_library_and_restored(app, tmp_path):
    import mark_library
    _two_runs_of_a01(app)
    key1 = app._axes_keys[1]
    app._hover_idx = 1
    app._nudge(10)
    app._nudge(5)
    app._flush_alignment()                    # the batched save, without waiting
    wells = mark_library.load_run("PLATE_X_Run02")
    assert wells["A01"]["align"] == [1.0, 15.0]
    # a fresh session on the same library finds it again
    app._align.clear()
    app.docs.clear()
    app._library_cache.clear()
    _two_runs_of_a01(app)
    assert app._align.get(key1) == (1.0, 15.0)
    assert "+15 scans" in app._axes_titles[1][0].get_text()
    # resetting removes it from the library too
    app._hover_idx = 1
    app._reset_alignment("hover")
    app._flush_alignment()
    assert "align" not in (mark_library.load_run("PLATE_X_Run02").get("A01") or {})


def test_saving_an_alignment_keeps_the_wells_marks(app):
    import mark_library
    mark_library.save_run("PLATE_X_Run02", {"A01": {"std": [[2100, "HOM1"]], "std_source": None}})
    app._library_cache.clear()
    _two_runs_of_a01(app)
    app._hover_idx = 1
    app._nudge(3)
    app._flush_alignment()
    entry = mark_library.load_run("PLATE_X_Run02")["A01"]
    assert entry["std"] == [[2100, "HOM1"]] and entry["align"] == [1.0, 3.0]


def test_in_overlay_a_graph_moves_with_all_its_runs(app):
    app.overlay_runs.set(True)
    _select(app, [0, 1], 1)                       # Run01/A01 + Run02/A01 -> one graph
    assert len(app._plot_axes) == 1
    k1, k2 = (app._akey(p) for p in app.selected)
    app._align[k2] = (1.0, 30.0)                  # Run02 already aligned onto Run01
    app._hover_idx = 0
    app._nudge(5)
    assert app._align[k1] == (1.0, 5.0)           # the graph moved ...
    assert app._align[k2] == (1.0, 35.0)          # ... and its overlaid run with it
    app._reset_alignment("hover")
    assert k1 not in app._align and k2 not in app._align


def test_ctce_window_analyses_a_plate_and_shows_a_well(app):
    app.show_ctce()
    app.update()
    u = app._ctce_ui
    key = next(k for k in u.names if "PLATE_X" in k)
    u.plate_var.set(key)
    u.analyse()
    rows = u.tree.get_children()
    assert len(rows) == 3                                   # A01..A03 in both runs
    assert app._ctce_marks                                  # classes to mark in the graphs
    u.tree.selection_set(rows[0])
    u.show_well()
    app.update()
    assert app.overlay_runs.get()
    assert {p.parent.name for p in app.selected} == {"PLATE_X_Run01", "PLATE_X_Run02"}
    u.win.destroy()
