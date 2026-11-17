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
