"""Internal-standard quartet detection (``genotyping.find_is_quartet``).

The geometry constraints are unit-tested on synthetic traces so the suite needs
no trace files.  The T9 plate cases are exercised too, but skipped when the
cached traces are not present -- they live outside the repository because the
plate is far too large to ship.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

import genotyping as g

CACHE = Path("/home/tv/t9_rs1695_analysis/data/t9raw.npz")
N_SCANS = 4000
BASELINE = 100.0


def synth(peaks, height=1000.0, n_scans=N_SCANS, noise=0.0, seed=0):
    """A trace with Gaussian peaks at *peaks* over a flat baseline."""
    rng = np.random.default_rng(seed)
    y = np.full(n_scans, BASELINE) + rng.normal(0.0, noise, n_scans)
    x = np.arange(n_scans)
    for p in peaks:
        y += height * np.exp(-0.5 * ((x - p) / 3.0) ** 2)
    return y


# d1=30, d2=85, d3=35 -> dm=32.5, |d1-d3|=5 <= 24.4, 85 in [55.25, 107.25]
GOOD = [2100, 2130, 2215, 2250]


def test_finds_a_conforming_quartet():
    found = g.find_is_quartet(synth(GOOD))
    assert found is not None
    assert found[0] == GOOD


def test_quartet_needs_both_outer_spacings():
    """Dropping the first peak apart must not still read as a standard."""
    assert g.find_is_quartet(synth([2100, 2110, 2215, 2250])) is None


def test_quartet_needs_the_right_middle_gap():
    """d2 is ~2.5x dm; a tight middle gap is a different fragment."""
    assert g.find_is_quartet(synth([2100, 2130, 2160, 2190])) is None


def test_quartet_needs_four_peaks():
    assert g.find_is_quartet(synth([2100, 2130, 2215])) is None


def test_equimolar_prior_rejects_a_lopsided_group():
    """One short peak among three tall ones is not a standard."""
    y = synth(GOOD)
    y = y - 1000.0 * np.exp(-0.5 * ((np.arange(N_SCANS) - 2250) / 3.0) ** 2)
    y += 60.0 * np.exp(-0.5 * ((np.arange(N_SCANS) - 2250) / 3.0) ** 2)
    assert g.find_is_quartet(y) is None


def test_plus_a_tail_is_dropped_not_mistaken_for_a_member():
    """A short +A product 10-32 scans behind its parent must not fill a slot."""
    peaks = GOOD + [2240]
    y = synth(peaks)
    # shrink the intruder at 2240 well below the equimolar floor
    x = np.arange(N_SCANS)
    y = y - 1000.0 * np.exp(-0.5 * ((x - 2240) / 3.0) ** 2)
    y += 100.0 * np.exp(-0.5 * ((x - 2240) / 3.0) ** 2)
    found = g.find_is_quartet(y)
    assert found is not None
    assert 2240 not in found[0]


def test_unequal_outer_spacings_rejected():
    """d1 and d3 must be comparable; 25 next to 120 is not a standard."""
    assert g.find_is_quartet(synth([2100, 2125, 2300, 2420])) is None


def test_flat_trace_returns_none():
    assert g.find_is_quartet(np.full(N_SCANS, BASELINE)) is None


def test_quartet_below_the_cut_is_ignored():
    """The standard migrates late; an early copy must not be picked up."""
    early = [p - 1900 for p in GOOD]
    found = g.find_is_quartet(synth(early + GOOD))
    assert found is not None and found[0] == GOOD


# --------------------------------------------------------------------------- #
# channel mapping
# --------------------------------------------------------------------------- #
def test_channel_to_acgt_index_on_a_megabace_plate():
    # "ACTG" dye order: Ch1=A, Ch2=C, Ch3=T, Ch4=G, all of which are already
    # in ACGT order here.
    assert g.acgt_index_for_channel("ACTG", 1) == 0
    assert g.acgt_index_for_channel("ACTG", 3) == 3
    # A different dye order must not silently return the same column.
    # "TGCA": Ch1=T, Ch2=G, Ch3=C, Ch4=A.
    assert g.acgt_index_for_channel("TGCA", 1) == 3
    assert g.acgt_index_for_channel("TGCA", 2) == 2
    assert g.acgt_index_for_channel("TGCA", 3) == 1
    assert g.acgt_index_for_channel("TGCA", 4) == 0


def test_channel_mapping_rejects_nonsense():
    with pytest.raises(ValueError):
        g.acgt_index_for_channel("ACTX", 1)
    with pytest.raises(ValueError):
        g.acgt_index_for_channel("ACTG", 5)


# --------------------------------------------------------------------------- #
# real plate traces
# --------------------------------------------------------------------------- #
needs_cache = pytest.mark.skipif(not CACHE.exists(),
                                 reason="T9 trace cache is not in the repository")


@pytest.mark.parametrize("well,expected", [
    ("A01", [2135, 2215, 2426, 2521]),
    ("A12", [2113, 2170, 2358, 2481]),
    ("H04", [2455, 2572, 2781, 2850]),
])
@needs_cache
def test_t9_plate_quadrets(well, expected):
    with np.load(CACHE) as z:
        found = g.find_is_quartet(z[well].astype(float)[:, 1])
    assert found is not None and found[0] == expected


@needs_cache
def test_every_t9_well_has_a_quartet():
    with np.load(CACHE) as z:
        missing = [w for w in z.files
                   if g.find_is_quartet(z[w].astype(float)[:, 1]) is None]
    assert missing == []
