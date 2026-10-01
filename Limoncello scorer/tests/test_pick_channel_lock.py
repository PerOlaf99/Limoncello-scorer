"""A click must be able to say which channel it meant.

The four channels are drawn overlaid in a single axes, so a click carries no
channel information.  `pick()` therefore inferred one, ranking candidates by
nearest apex scan first and falling back to voltage.  On the ABCC2 plate the
internal standard appears on two channels a few scans apart, and the taller
one won: a click aimed at Ch3 was recorded as Ch4 (it wrote ``channel=4`` into
the exported peak table for wells the operator had clicked on Ch3).

`only_col` makes the caller state the channel outright.  These tests pin that
lock, and that without it the old behaviour is still the documented fallback.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import genotyping as g

N_SCANS = 4000


class _Doc:
    def __init__(self, acgt, base_order="ACTG"):
        self.acgt = acgt
        self.n_scans = acgt.shape[0]
        self.base_order = base_order
        self.path = Path("/tmp/synthetic.rsd")
        self.well = "A01"
        self.sequence = ""
        self.qualities = []
        self.peak_positions = []
        self.raw = acgt


def _two_channel_trace():
    """A peak on Ch1 (acgt col 0) sitting 12 scans off the one on Ch2 (col 1).

    The gap and the height ratio mirror the real ABCC2 case, where the second
    channel's apex lands exactly on the click scan while the first channel's
    apex is a dozen scans away.
    """
    acgt = np.full((N_SCANS, 4), 50.0)
    x = np.arange(N_SCANS)
    sigma = 4.0
    # col 0: apex at 2388.  col 1: apex at 2400, taller.
    acgt[:, 0] += 600.0 * np.exp(-0.5 * ((x - 2388) / sigma) ** 2)
    acgt[:, 1] += 950.0 * np.exp(-0.5 * ((x - 2400) / sigma) ** 2)
    return acgt


def _picker():
    return g.PeakPicker(_Doc(_two_channel_trace()), Path("/tmp/synthetic.rsd"))


def test_locking_to_a_column_picks_that_column():
    pk = _picker()
    rec = pk.pick(2400, vol=600.0, only_col=0)
    assert rec is not None, "clicking at 2400 on col 0 should find its peak"
    # The record's column is the acgt column; channel is col+1.
    assert rec["col"] == 0
    assert rec["scan"] == 2388


def test_locking_refuses_to_cross_to_the_other_channel():
    """A lock must not be overridden by a nearer or taller peak elsewhere."""
    pk = _picker()
    rec = pk.pick(2400, vol=600.0, only_col=0)
    assert rec["col"] == 0, "must not drift onto col 1 despite its apex at 2400"


def test_unlocked_pick_prefers_the_nearest_apex():
    """Without a lock the old nearest-apex behaviour still applies.

    This is the behaviour that produced the wrong channel in the export, so it
    is pinned deliberately: changing it later would change every unmarked well.
    """
    pk = _picker()
    rec = pk.pick(2400, vol=600.0)
    assert rec["col"] == 1 and rec["scan"] == 2400


def test_locking_reports_no_peak_when_that_channel_has_none():
    acgt = np.full((N_SCANS, 4), 50.0)
    x = np.arange(N_SCANS)
    acgt[:, 0] += 900.0 * np.exp(-0.5 * ((x - 2400) / 4.0) ** 2)
    pk = g.PeakPicker(_Doc(acgt), Path("/tmp/synthetic.rsd"))
    assert pk.pick(2400, vol=900.0, only_col=3) is None
    assert pk._reject == "none"


def test_exported_channel_reflects_the_locked_column():
    """The exported `channel` is what went wrong in the operator's CSV."""
    pk = _picker()
    pk.pick(2400, vol=600.0, only_col=0)
    row = pk.export_rows()[0]
    assert row["channel"] == 1 and row["base"] == "A"