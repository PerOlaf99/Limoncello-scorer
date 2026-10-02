"""The channel bar must always read Ch1, Ch2, Ch3, Ch4 left to right.

The widgets are created in matrix-column order and widget k drives
chan_show[k], which is indexed by physical channel, so the widgets are already
in the right places.  What went wrong was the label: it was keyed off the
matrix column instead of off the channel, so with the degenerate ACTG dye
order the bar read "Ch1 A, Ch2 C, Ch4 G, Ch3 T" and Ch4 sat to the left of
Ch3.  These pin the labels, and the invariant that the channel numbers always
ascend whatever the dye order is.
"""
import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import sequence_analyzer as sa  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))


class _Var:
    """Stands in for a tk.StringVar holding a dye order."""

    def __init__(self):
        self.v = ""

    def set(self, value):
        self.v = value

    def get(self):
        return self.v


class _Bar:
    """Just enough app to call channel_bar_order without a display."""

    def __init__(self, order="ACTG"):
        self.base_order_var = _Var()
        self.base_order_var.set(order)

    channel_bar_order = sa.LimoncelloAnalyzerApp.channel_bar_order
    _col_to_chan = sa.LimoncelloAnalyzerApp._col_to_chan


def labels(order):
    return [t.strip() for _ci, t in _Bar(order).channel_bar_order()]


def channels(order):
    return [ci for ci, _t in _Bar(order).channel_bar_order()]


def test_megabace_actg_reads_in_channel_order():
    assert labels("ACTG") == ["Ch1 A", "Ch2 C", "Ch3 T", "Ch4 G"]


def test_t9_actg_puts_the_t_internal_standard_on_channel_three():
    # The T9 DyeSet2 traces validate against ACTG: the internal standard is the
    # T column (lowest CV across the 96 wells), which ACTG places on Ch3.
    assert labels("ACTG")[2] == "Ch3 T"


def test_dyenamic_tgca_reads_in_channel_order_too():
    assert labels("TGCA") == ["Ch1 T", "Ch2 G", "Ch3 C", "Ch4 A"]


def test_gcat_order_reads_in_channel_order():
    assert labels("GCAT") == ["Ch1 G", "Ch2 C", "Ch3 A", "Ch4 T"]


@pytest.mark.parametrize("order", ["ACTG", "ACGT", "TGCA", "GATC", "CTAG"])
def test_channel_numbers_always_ascend(order):
    assert channels(order) == [0, 1, 2, 3]


@pytest.mark.parametrize("order", ["ACTG", "ACGT", "TGCA", "GATC", "CTAG"])
def test_the_labelled_ch_is_the_one_it_toggles(order):
    bar = _Bar(order)
    for ci, text in bar.channel_bar_order():
        m = re.search(r"Ch(\d+) (\w)", text)
        assert m and int(m.group(1)) == ci + 1


def test_no_order_given_falls_back_to_actg():
    bar = _Bar("ACTG")
    bar.base_order_var.set("")
    assert [t.strip() for _c, t in bar.channel_bar_order()] == \
        ["Ch1 A", "Ch2 C", "Ch3 T", "Ch4 G"]


def test_col_to_chan_maps_matrix_columns_onto_channels():
    # ACTG: matrix is always A,C,G,T but the dyes sit A,C,T,G on Ch1..Ch4.
    assert _Bar("ACTG")._col_to_chan() == {0: 0, 1: 1, 3: 2, 2: 3}
    assert _Bar("ACGT")._col_to_chan() == {0: 0, 1: 1, 2: 2, 3: 3}


def test_the_bar_is_not_repacked_away_from_its_header():
    """Relabelling must not touch packing.  The checkboxes share their parent
    with side=RIGHT siblings, so pack_forget + pack would move the whole group
    to the end of the toolbar, away from the "Channels:" label."""
    src = sa.LimoncelloAnalyzerApp._refresh_channel_labels.__doc__ + "\n" + \
        __import__("inspect").getsource(sa.LimoncelloAnalyzerApp._refresh_channel_labels)
    assert "pack_forget" not in src
    assert "pack(" not in src
