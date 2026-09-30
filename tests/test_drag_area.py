"""Drag-to-measure: baseline-corrected areas from hand-placed spans.

``region_area`` is pure arithmetic so it is pinned with exact numbers rather
than tolerances -- a regression in the baseline would otherwise hide behind a
loose approx().
"""
import csv
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import genotyping as g

RUN = Path("/data/T9_run01/A01.rsd")


class _Doc:
    """Minimal stand-in for analyzer_core.TraceDocument."""

    def __init__(self, acgt, well="A01"):
        self.acgt = acgt
        self.n_scans = acgt.shape[0]
        self.base_order = "ACTG"
        self.well = well
        self.path = RUN
        self.sequence = ""
        self.qualities = []
        self.peak_positions = []
        self.raw = acgt


def _flat_trace(y, n_scans=None, col=0):
    n = n_scans or len(y)
    acgt = np.zeros((n, 4))
    acgt[:, col] = y
    return acgt


# --------------------------------------------------------------------------- #
# region_area
# --------------------------------------------------------------------------- #
def test_area_is_the_trace_above_a_flat_baseline():
    res = g.region_area([0, 0, 3, 0, 0], 0, 4)
    assert res["area"] == pytest.approx(3.0)
    assert res["midpoint"] == 2
    assert res["peak_scan"] == 2
    assert res["height"] == pytest.approx(3.0)


def test_baseline_follows_the_endpoint_voltages():
    """Endpoints lifted off zero: the baseline is the line through them, not
    the trace's own minimum."""
    res = g.region_area([2, 2, 5, 2, 2], 0, 4)
    assert res["area"] == pytest.approx(3.0)
    assert res["height"] == pytest.approx(3.0)
    assert res["baseline_left"] == pytest.approx(2.0)
    assert res["baseline_right"] == pytest.approx(2.0)


def test_sloped_baseline_is_interpolated_between_the_two_points():
    res = g.region_area([0, 1, 4, 3, 2], 0, 4)
    # baseline = linspace(0, 2, 5); differences = 0, .5, 3, 1.5, 0
    assert res["area"] == pytest.approx(5.0)
    assert res["peak_scan"] == 2
    assert res["height"] == pytest.approx(3.0)


def test_area_never_counts_area_below_the_baseline():
    """A dip under the line contributes nothing, it is not subtracted."""
    res = g.region_area([0, 3, 0, 0, 0], 0, 4)
    assert res["area"] == pytest.approx(3.0)
    assert res["peak_scan"] == 1


def test_reversed_span_is_normalised():
    assert g.region_area([0, 1, 4, 3, 2], 4, 0) == g.region_area(
        [0, 1, 4, 3, 2], 0, 4)


def test_span_shorter_than_the_minimum_is_refused():
    assert g.region_area([0, 0, 9, 0], 0, 2) is None


def test_flat_trace_has_no_area():
    assert g.region_area([2, 2, 2, 2, 2], 0, 4) is None


def test_span_is_clamped_to_the_trace():
    res = g.region_area([0, 1, 4, 3, 2], -50, 500)
    assert res is not None
    assert res["start"] == 0 and res["stop"] == 4


def test_empty_trace_is_refused():
    assert g.region_area([], 0, 4) is None


# --------------------------------------------------------------------------- #
# DragAreaPicker
# --------------------------------------------------------------------------- #
def _picker(acgt, **kw):
    return g.DragAreaPicker(_Doc(acgt), RUN, **kw)


def test_run_and_well_come_from_the_path_and_document():
    pk = _picker(_flat_trace([0, 0, 5, 0, 0]))
    assert pk.run_name == "T9_run01"
    assert pk.well == "A01"


def test_run_name_can_be_overridden():
    pk = _picker(_flat_trace([0, 0, 5, 0, 0]), run_name="custom run")
    assert pk.run_name == "custom run"


def test_add_records_midpoint_and_area():
    pk = _picker(_flat_trace([0, 0, 5, 0, 0]))
    rec = pk.add(0, 4)
    assert rec["midpoint"] == 2
    assert rec["area_Vscan"] == pytest.approx(5.0)
    assert rec["start_scan"] == 0 and rec["end_scan"] == 4
    assert rec["run"] == "T9_run01"
    assert rec["well"] == "A01"


def test_auto_channel_picks_the_tallest_hump():
    acgt = _flat_trace([0, 0, 1, 0, 0], col=0)
    acgt[:, 2] = [0, 0, 9, 0, 0]
    pk = _picker(acgt)
    assert pk.best_col(0, 4) == 2
    assert pk.add(0, 4)["base"] == "G"


def test_explicit_channel_overrides_the_auto_choice():
    acgt = _flat_trace([0, 0, 9, 0, 0], col=0)
    acgt[:, 2] = [0, 0, 1, 0, 0]
    pk = _picker(acgt)
    rec = pk.add(0, 4, col=2)
    assert rec["base"] == "G"
    assert rec["area_Vscan"] == pytest.approx(1.0)


def test_auto_channel_ignores_hidden_channels():
    acgt = _flat_trace([0, 0, 1, 0, 0], col=0)
    acgt[:, 2] = [0, 0, 9, 0, 0]
    pk = _picker(acgt, show=lambda c: c != 2)
    assert pk.best_col(0, 4) == 0


def test_add_refuses_a_degenerate_span():
    pk = _picker(_flat_trace([0, 0, 5, 0, 0]))
    assert pk.add(1, 2) is None
    assert pk.rows() == []


def test_add_refuses_an_out_of_range_channel():
    pk = _picker(_flat_trace([0, 0, 5, 0, 0]))
    assert pk.add(0, 4, col=9) is None


def test_undo_and_clear():
    pk = _picker(_flat_trace([0, 0, 5, 0, 0]))
    pk.add(0, 4)
    pk.add(0, 4)
    assert len(pk.rows()) == 2
    pk.undo_last()
    assert len(pk.rows()) == 1
    pk.undo_last()
    assert pk.rows() == []
    assert pk.undo_last() is None
    pk.add(0, 4)
    pk.clear_all()
    assert pk.rows() == []


def test_table_rows_keep_click_order_but_export_is_left_to_right():
    # two separate humps, so both spans have something to measure
    pk = _picker(_flat_trace([0, 0, 5, 0, 0, 0, 0, 0, 7, 0, 0]))
    pk.add(6, 10)
    pk.add(0, 4)
    assert [r["start_scan"] for r in pk.rows()] == [6, 0]
    assert [r["start_scan"] for r in pk.export_rows()] == [0, 6]


def test_rows_expose_the_columns_the_table_shows():
    pk = _picker(_flat_trace([0, 0, 5, 0, 0]))
    row = pk.add(0, 4)
    for key in ("run", "well", "midpoint", "area_Vscan", "start_scan",
                "end_scan", "base", "height_V"):
        assert key in row


def test_measurements_save_as_csv(tmp_path):
    pk = _picker(_flat_trace([0, 0, 5, 0, 0]))
    pk.add(0, 4)
    out = tmp_path / "areas.csv"
    g.save_table(out, pk.export_rows())
    with open(out, newline="", encoding="utf-8-sig") as fh:
        got = list(csv.DictReader(fh))
    assert len(got) == 1
    assert got[0]["well"] == "A01"
    assert got[0]["run"] == "T9_run01"
    assert got[0]["midpoint"] == "2"
    assert float(got[0]["area_Vscan"]) == pytest.approx(5.0)


def test_saving_nothing_still_writes_a_header_free_file(tmp_path):
    pk = _picker(_flat_trace([0, 0, 5, 0, 0]))
    out = tmp_path / "empty.csv"
    g.save_table(out, pk.export_rows())
    assert out.read_text(encoding="utf-8-sig").strip() == ""


# --------------------------------------------------------------------------- #
# drawing
# --------------------------------------------------------------------------- #
def test_overlay_draws_a_baseline_per_measurement():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from matplotlib.figure import Figure

    pk = _picker(_flat_trace([0, 0, 5, 0, 0]))
    pk.add(0, 4)
    fig = Figure()
    ax = fig.add_subplot(111)
    x = np.arange(5)
    for col, color in sorted(pk.col_color.items()):
        ax.plot(x, pk.doc.acgt[:, col], color=color, lw=0.7)
    pk.plot_overlay(ax)
    # one trace per channel, plus the baseline line and the two guides
    assert len(ax.lines) == len(pk.col_color) + 1 + 2
    assert len(ax.collections) == 1          # the shaded area


def test_overlay_shows_the_span_being_dragged():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from matplotlib.figure import Figure

    pk = _picker(_flat_trace([0, 0, 5, 0, 0]))
    fig = Figure()
    ax = fig.add_subplot(111)
    x = np.arange(5)
    for col, color in sorted(pk.col_color.items()):
        ax.plot(x, pk.doc.acgt[:, col], color=color, lw=0.7)
    pk.plot_overlay(ax, pending=(0, 4, 0))
    assert len(ax.lines) == len(pk.col_color) + 1 + 2


# --------------------------------------------------------------------------- #
# the main-window drag state machine
# --------------------------------------------------------------------------- #
# The tool now lives inside the Limoncello window, driven by the same handlers
# the peak-picking mode uses, so the state machine is exercised directly
# against a stub.  press/motion/release is new in this repo, and the GUI cannot
# be launched headlessly, so this is the only place it gets covered.
import sequence_analyzer as sa


class _Ev:
    def __init__(self, x, button=1, inaxes=True):
        self.xdata = x
        self.ydata = 0.0
        self.inaxes = object() if inaxes else None
        self.button = button
        self.x = None
        self.y = None


class _Status:
    def __init__(self):
        self.msgs = []
        self.text = ""

    def set(self, msg):
        self.msgs.append(msg)

    def config(self, **kw):          # stands in for a ttk.Label header
        if "text" in kw:
            self.text = kw["text"]


class _Canvas:
    def __init__(self):
        self.idles = 0

    def draw_idle(self):
        self.idles += 1

    def mpl_disconnect(self, cid):
        pass

    def pack_forget(self):
        pass

    def pack(self, **kw):
        pass

    def destroy(self):
        pass

    def mpl_connect(self, event, cb):
        self.ids = getattr(self, "ids", 0) + 1
        self.connected = getattr(self, "connected", [])
        self.connected.append((event, cb))
        return self.ids


class _Stub:
    """Enough of LimoncelloAnalyzerApp for the three area-mode handlers."""

    def __init__(self, acgt, path=RUN, area_mode=True, active=True):
        self.area_mode = area_mode
        self.genotyping_active = active
        self._area_drag = None
        self._gen_dursors = []
        self._area_cursors = []
        self._gen_pickers = {}
        self._area_pickers = {}
        self._gen_active_path = None
        self._paths = [path]
        self.status_var = _Status()
        self.canvas = _Canvas()
        self.redraws = 0
        self.syncs = 0
        self._doc = _Doc(acgt)

    def _gen_paths(self):
        return list(self._paths)

    def _ensure_area_picker(self, path):
        key = str(Path(path).resolve())
        pk = self._area_pickers.get(key)
        if pk is None:
            pk = g.DragAreaPicker(self._doc, path)
            self._area_pickers[key] = pk
        return pk

    def _gen_axes_hit(self, event):
        if event.xdata is None:
            return None
        return 0, None, event.xdata, 0.0

    def redraw(self):
        self.redraws += 1

    def _sync_pick_table(self):
        self.syncs += 1

    def records(self):
        pk = self._ensure_area_picker(self._paths[0])
        return pk.rows()


def _driver(stub):
    for name in ("_on_gen_pick", "_on_gen_motion", "_on_area_release",
                 "_area_preview"):
        setattr(stub, name,
                getattr(sa.LimoncelloAnalyzerApp, name).__get__(stub))


HUMP = _flat_trace([0, 0, 5, 0, 0])


def test_drag_records_one_measurement_on_release():
    st = _Stub(HUMP)
    _driver(st)
    st._on_gen_pick(_Ev(0))
    assert st.records() == []              # nothing stored while still dragging
    st._on_gen_motion(_Ev(2))
    st._on_area_release(_Ev(4))
    assert len(st.records()) == 1
    assert st.records()[0]["midpoint"] == 2
    assert st.syncs == 1                   # the table is refreshed on release


def test_drag_right_to_left_works():
    st = _Stub(HUMP)
    _driver(st)
    st._on_gen_pick(_Ev(4))
    st._on_area_release(_Ev(0))
    assert st.records()[0]["start_scan"] == 0


def test_release_uses_its_own_coordinate_when_motion_never_arrived():
    """A brisk drag can finish before any motion event; the release position is
    the reliable end of the span."""
    st = _Stub(HUMP)
    _driver(st)
    st._on_gen_pick(_Ev(0))
    st._on_area_release(_Ev(4))
    assert len(st.records()) == 1
    assert st.records()[0]["end_scan"] == 4


def test_click_without_drag_is_refused():
    st = _Stub(HUMP)
    _driver(st)
    st._on_gen_pick(_Ev(2))
    st._on_area_release(_Ev(2))
    assert st.records() == []
    assert any("at least" in m for m in st.status_var.msgs)
    assert st._area_drag is None           # the drag is finished either way


def test_motion_without_press_is_ignored():
    st = _Stub(HUMP)
    _driver(st)
    st._on_gen_motion(_Ev(3))
    st._on_area_release(_Ev(3))
    assert st.records() == []


def test_release_without_press_is_ignored():
    st = _Stub(HUMP)
    _driver(st)
    st._on_area_release(_Ev(4))
    assert st.records() == []


def test_right_button_does_not_start_a_drag():
    st = _Stub(HUMP)
    _driver(st)
    st._on_gen_pick(_Ev(0, button=3))
    assert st._area_drag is None
    st._on_area_release(_Ev(4, button=3))
    assert st.records() == []


def test_click_outside_the_axes_is_ignored():
    st = _Stub(HUMP)
    _driver(st)
    st._on_gen_pick(_Ev(0, inaxes=False))
    assert st._area_drag is None


def test_drag_over_flat_trace_explains_itself():
    st = _Stub(_flat_trace([2, 2, 2, 2, 2, 2, 2]))
    _driver(st)
    st._on_gen_pick(_Ev(0))
    st._on_area_release(_Ev(6))
    assert st.records() == []
    assert any("baseline" in m for m in st.status_var.msgs)


def test_tallest_channel_is_measured():
    acgt = _flat_trace([0, 0, 1, 0, 0], col=0)
    acgt[:, 1] = [0, 0, 8, 0, 0]
    st = _Stub(acgt)
    _driver(st)
    st._on_gen_pick(_Ev(0))
    st._on_area_release(_Ev(4))
    assert st.records()[0]["base"] == "C"


def test_two_drags_make_two_rows():
    st = _Stub(_flat_trace([0, 0, 5, 0, 0, 0, 0, 0, 7, 0, 0]))
    _driver(st)
    for a, b in ((0, 4), (6, 10)):
        st._on_gen_pick(_Ev(a))
        st._on_area_release(_Ev(b))
    assert [r["midpoint"] for r in st.records()] == [2, 8]


def test_peak_picking_mode_is_untouched_by_the_drag_handlers():
    """In the other mode a press picks a peak; it must not arm a span."""
    st = _Stub(HUMP, area_mode=False)
    _driver(st)
    st._on_gen_pick(_Ev(0))
    assert st._area_drag is None


def test_handlers_do_nothing_outside_the_mode():
    st = _Stub(HUMP, active=False)
    _driver(st)
    st._on_gen_pick(_Ev(0))
    assert st._area_drag is None
    st._on_area_release(_Ev(4))
    assert st.records() == []


# --------------------------------------------------------------------------- #
# the live preview
# --------------------------------------------------------------------------- #
def _preview_stub(acgt):
    """A stub with a real Agg figure behind it, so the preview artists exist."""
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from matplotlib.figure import Figure

    st = _Stub(acgt)
    _driver(st)
    fig = Figure()
    ax = fig.add_subplot(111)
    n = acgt.shape[0]
    for col, color in sorted(st._ensure_area_picker(RUN).col_color.items()):
        ax.plot(np.arange(n), acgt[:, col], color=color, lw=0.7)
    st._area_cursors = [sa.LimoncelloAnalyzerApp._make_area_markers(st, ax)]
    st._ax = ax
    return st


def test_preview_shows_the_span_being_dragged():
    st = _preview_stub(HUMP)
    st._on_gen_pick(_Ev(0))
    st._area_preview(_Ev(4))
    _ax, line, left, right, ann = st._area_cursors[0]
    assert line.get_visible()
    assert left.get_visible() and right.get_visible()
    assert list(left.get_xdata()) == [0, 0]
    assert list(right.get_xdata()) == [4, 4]
    assert "mid 2" in ann.get_text()
    assert "5.0" in ann.get_text()          # the live area


def test_preview_hides_again_when_the_span_is_too_short():
    st = _preview_stub(HUMP)
    st._on_gen_pick(_Ev(2))
    st._area_preview(_Ev(3))
    for art in st._area_cursors[0][1:]:
        assert not art.get_visible()


def test_release_clears_the_preview():
    st = _preview_stub(HUMP)
    st._on_gen_pick(_Ev(0))
    st._area_preview(_Ev(4))
    st._on_area_release(_Ev(4))
    for art in st._area_cursors[0][1:]:
        assert not art.get_visible()
    assert len(st.records()) == 1


def test_motion_drives_the_preview_without_a_full_redraw():
    st = _preview_stub(HUMP)
    st._on_gen_pick(_Ev(0))
    st._on_gen_motion(_Ev(4))
    assert st.redraws == 0                  # preview artists, not a rebuild
    _ax, line, _l, right, _a = st._area_cursors[0]
    assert list(right.get_xdata()) == [4, 4]


class _Var:
    """Stands in for a tk.BooleanVar: .get() reports the current flag."""

    def __init__(self):
        self.v = False

    def set(self, value):
        self.v = bool(value)

    def get(self):
        return self.v


def _mode_stub():
    """Bind the real mode handlers to a stub carrying the few attributes they
    touch, so the menu indicators can be watched across real transitions."""
    s = _Stub(HUMP, area_mode=False, active=False)
    s._mode_pick = _Var()
    s._mode_area = _Var()
    s.seq_hdr = _Status()
    s.seq_text = _Canvas()
    s.selected = [RUN]          # the "pick some wells first" guard
    s.pick_table = None
    s._build_pick_table = lambda area=False: None
    for name in ("enter_genotyping_picking", "enter_area_picking",
                 "exit_genotyping_picking", "toggle_genotyping_picking",
                 "toggle_area_picking", "_on_gen_pick", "_on_gen_motion",
                 "_on_area_release"):
        setattr(s, name, getattr(sa.LimoncelloAnalyzerApp, name).__get__(s))
    return s


def _indicators(s):
    return s._mode_pick.get(), s._mode_area.get()


def test_neither_mode_is_active_before_anything_is_chosen():
    assert _indicators(_mode_stub()) == (False, False)


def test_entering_peak_picking_ticks_only_its_own_box():
    s = _mode_stub()
    s.enter_genotyping_picking()
    assert _indicators(s) == (True, False)
    assert s.area_mode is False


def test_entering_drag_mode_ticks_only_its_own_box():
    s = _mode_stub()
    s.enter_area_picking()
    assert _indicators(s) == (False, True)
    assert s.area_mode is True


def test_the_two_modes_can_never_both_be_ticked():
    s = _mode_stub()
    for step in (s.toggle_genotyping_picking, s.toggle_area_picking,
                 s.toggle_area_picking, s.toggle_genotyping_picking,
                 s.toggle_genotyping_picking, s.toggle_area_picking,
                 s.exit_genotyping_picking):
        step()
        assert sum(_indicators(s)) <= 1, (step, _indicators(s))


def test_switching_modes_leaves_the_old_one_unticked():
    s = _mode_stub()
    s.toggle_genotyping_picking()
    assert _indicators(s) == (True, False)
    s.toggle_area_picking()          # ask for drag while picking
    assert _indicators(s) == (False, True)
    assert s.seq_hdr.text == "Measured areas"   # the mode really did switch
    s.toggle_genotyping_picking()    # and back again
    assert _indicators(s) == (True, False)
    assert s.seq_hdr.text == "Picked peaks"


def test_unticking_the_live_mode_returns_to_the_plain_viewer():
    s = _mode_stub()
    s.toggle_area_picking()
    assert _indicators(s) == (False, True)
    s.toggle_area_picking()          # same one clicked a second time
    assert _indicators(s) == (False, False)
    assert s.genotyping_active is False and s.area_mode is False


def test_exiting_clears_both_boxes():
    s = _mode_stub()
    s.toggle_genotyping_picking()
    s.exit_genotyping_picking()
    assert _indicators(s) == (False, False)
