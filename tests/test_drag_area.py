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
# the drag state machine
# --------------------------------------------------------------------------- #
# The GUI cannot be launched headlessly, but the handlers only touch a handful
# of attributes, so they are driven directly against a stub.  This is the part
# with no precedent in the repo -- button_release_event is used nowhere else --
# so the press/motion/release sequence deserves to be pinned explicitly.
class _Ev:
    def __init__(self, x, button=1, inaxes=True):
        self.xdata = x
        self.ydata = 0.0
        self.inaxes = object() if inaxes else None
        self.button = button


class _Toolbar:
    mode = ""


class _Var:
    def __init__(self, v):
        self.v = v

    def get(self):
        return self.v


class _Stub:
    """Just enough of ManualAreaEditor for the three mouse handlers."""

    def __init__(self, acgt, chan="auto", mode=""):
        self.pk = g.DragAreaPicker(_Doc(acgt), RUN)
        self.toolbar = _Toolbar()
        self.toolbar.mode = mode
        self.chan = _Var(chan)
        self._drag_from = None
        self._drag_to = None
        self.redraws = 0
        self.statuses = []
        self._syncs = 0

    def redraw(self):
        self.redraws += 1

    def _draw_pending(self):
        if self._drag_from is None or self._drag_to is None:
            return
        if abs(self._drag_to - self._drag_from) >= g.DRAG_MIN_SPAN:
            self.redraws += 1

    def _sync_table(self):
        self._syncs += 1

    def _status(self, msg):
        self.statuses.append(msg)

    def _chan_col(self):
        v = self.chan.get()
        return None if v == "auto" else g.CHANNEL_ORDER.index(v)


def _driver(stub):
    for name in ("_on_press", "_on_motion", "_on_release"):
        setattr(stub, name, getattr(g.ManualAreaEditor, name).__get__(stub))


HUMP = _flat_trace([0, 0, 5, 0, 0])


def test_drag_records_one_measurement_on_release():
    st = _Stub(HUMP)
    _driver(st)
    st._on_press(_Ev(0))
    assert st.pk.rows() == []          # nothing recorded while still dragging
    st._on_motion(_Ev(2))
    st._on_release(_Ev(4))
    assert len(st.pk.rows()) == 1
    assert st.pk.rows()[0]["midpoint"] == 2
    assert st._syncs == 1              # the table is refreshed on release


def test_dragging_right_to_left_works():
    st = _Stub(HUMP)
    _driver(st)
    st._on_press(_Ev(4))
    st._on_release(_Ev(0))
    assert st.pk.rows()[0]["start_scan"] == 0


def test_click_without_drag_is_refused():
    st = _Stub(HUMP)
    _driver(st)
    st._on_press(_Ev(2))
    st._on_release(_Ev(2))
    assert st.pk.rows() == []
    assert any("at least" in m for m in st.statuses)
    assert st._drag_from is None        # the drag is finished either way


def test_motion_without_press_is_ignored():
    st = _Stub(HUMP)
    _driver(st)
    st._on_motion(_Ev(3))
    st._on_release(_Ev(3))
    assert st.pk.rows() == []


def test_release_without_press_is_ignored():
    st = _Stub(HUMP)
    _driver(st)
    st._on_release(_Ev(4))
    assert st.pk.rows() == []


def test_right_button_does_not_start_a_drag():
    st = _Stub(HUMP)
    _driver(st)
    st._on_press(_Ev(0, button=3))
    assert st._drag_from is None
    st._on_release(_Ev(4, button=3))
    assert st.pk.rows() == []


def test_toolbar_pan_and_zoom_win_over_dragging():
    for mode in ("pan/zoom", "zoom rect"):
        st = _Stub(HUMP, mode=mode)
        _driver(st)
        st._on_press(_Ev(0))
        assert st._drag_from is None, mode
        st._on_release(_Ev(4))
        assert st.pk.rows() == [], mode


def test_click_outside_the_axes_is_ignored():
    st = _Stub(HUMP)
    _driver(st)
    st._on_press(_Ev(0, inaxes=False))
    assert st._drag_from is None


def test_drag_over_flat_trace_reports_why():
    st = _Stub(_flat_trace([2, 2, 2, 2, 2, 2, 2]))
    _driver(st)
    st._on_press(_Ev(0))
    st._on_release(_Ev(6))
    assert st.pk.rows() == []
    assert any("baseline" in m for m in st.statuses)


def test_chosen_channel_is_used():
    acgt = _flat_trace([0, 0, 1, 0, 0], col=0)
    acgt[:, 1] = [0, 0, 8, 0, 0]
    st = _Stub(acgt, chan="C")
    _driver(st)
    st._on_press(_Ev(0))
    st._on_release(_Ev(4))
    assert st.pk.rows()[0]["base"] == "C"


def test_two_drags_make_two_rows():
    acgt = _flat_trace([0, 0, 5, 0, 0, 0, 0, 0, 7, 0, 0])
    st = _Stub(acgt)
    _driver(st)
    for a, b in ((0, 4), (6, 10)):
        st._on_press(_Ev(a))
        st._on_release(_Ev(b))
    assert [r["midpoint"] for r in st.pk.rows()] == [2, 8]
