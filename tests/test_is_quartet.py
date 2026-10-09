"""Internal-standard detection (``genotyping.find_is_quartet``).

Detection comes in two tiers, and both are unit-tested here on synthetic traces
so the suite needs no trace files.

**With a learned geometry** the quartet must match a shape fitted from the
operator's marks.  This is what runs in production, via
:class:`genotyping.PlateISModel`.

**Without one**, only equimolarity and "all the bands are one fragment, so they
cluster" are asserted -- see ``_search_equimolar_group``.  No spacing is
assumed, because spacing is exactly what varies between fragments.  A previous
version hard-coded T9's d1=25, d3=25, d2=2.5*d1 here as a fallback and it was
wrong on every other assay measured: ABCC2's d3 is 12 and its d2/d1 is 3.36.
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


# A shape fitted from marks, not a constant: d1=30, d2=85, d3=35.
GOOD = [2100, 2130, 2215, 2250]
LEARNED = g.ISGeometry(30.0, 85.0, 35.0, tol1=8.0, tol2=20.0, tol3=12.0)


# --------------------------------------------------------------------------- #
# with a learned geometry: the shape is what decides
# --------------------------------------------------------------------------- #
def test_finds_a_quartet_matching_the_learned_shape():
    found = g.find_is_quartet(synth(GOOD), geometry=LEARNED)
    assert found is not None
    assert found[0] == GOOD


def test_quartet_needs_both_outer_spacings():
    """Dropping the first peak apart must not still read as a standard."""
    assert g.find_is_quartet(synth([2100, 2110, 2215, 2250]),
                             geometry=LEARNED) is None


def test_quartet_needs_the_right_middle_gap():
    """d2 is a gap of this fragment's own; a tight one is a different fragment."""
    assert g.find_is_quartet(synth([2100, 2130, 2160, 2190]),
                             geometry=LEARNED) is None


def test_unequal_outer_spacings_rejected():
    """d1 and d3 must match the learned shape; 25 next to 120 is not it."""
    assert g.find_is_quartet(synth([2100, 2125, 2300, 2420]),
                             geometry=LEARNED) is None


def test_a_different_fragment_is_found_by_its_own_shape():
    """The point of learning the shape: a fragment T9 rejected is findable.

    ABCC2's real spacings (d1=76, d2=257, d3=12) fail every T9 constant, which
    is why no constant is used any more.
    """
    abcc2 = g.ISGeometry(76.0, 257.0, 12.0, tol1=15.0, tol2=45.0, tol3=6.0)
    peaks = [2400, 2476, 2733, 2745]
    found = g.find_is_quartet(synth(peaks), geometry=abcc2)
    assert found is not None and found[0] == peaks


def test_equimolar_prior_drops_the_short_peak_not_the_standard():
    """One short peak among three tall ones: it is dropped, not promoted.

    The three survivors are still a valid read -- a standard whose
    heteroduplexes co-migrate is three bands, which is the normal case. The
    failure to guard against is the *short* peak being accepted as a fourth.
    """
    x = np.arange(N_SCANS)
    y = synth(GOOD)
    y = y - 1000.0 * np.exp(-0.5 * ((x - 2250) / 3.0) ** 2)
    y += 60.0 * np.exp(-0.5 * ((x - 2250) / 3.0) ** 2)
    found = g.find_is_quartet(y, geometry=LEARNED)
    assert found is not None
    assert 2250 not in found[0]


def test_shape_free_search_drops_the_short_peak_too():
    x = np.arange(N_SCANS)
    y = synth(GOOD)
    y = y - 1000.0 * np.exp(-0.5 * ((x - 2250) / 3.0) ** 2)
    y += 60.0 * np.exp(-0.5 * ((x - 2250) / 3.0) ** 2)
    found = g.find_is_quartet(y)
    assert found is not None
    assert 2250 not in found[0]
    assert len(found[0]) == 3


def test_plus_a_tail_is_dropped_not_mistaken_for_a_member():
    """A short +A product 10-32 scans behind its parent must not fill a slot."""
    y = synth(GOOD + [2240])
    x = np.arange(N_SCANS)
    y = y - 1000.0 * np.exp(-0.5 * ((x - 2240) / 3.0) ** 2)
    y += 100.0 * np.exp(-0.5 * ((x - 2240) / 3.0) ** 2)
    found = g.find_is_quartet(y, geometry=LEARNED)
    assert found is not None
    assert 2240 not in found[0]


def test_three_bands_are_a_complete_match_when_merged():
    """Co-migrating heteroduplexes are the normal case, not a degraded one."""
    merged = g.ISGeometry(30.0, 85.0, None, tol1=8.0, tol2=20.0)
    found = g.find_is_quartet(synth([2100, 2130, 2215]), geometry=merged)
    assert found is not None and found[0] == [2100, 2130, 2215]


def test_flat_trace_returns_none():
    assert g.find_is_quartet(np.full(N_SCANS, BASELINE), geometry=LEARNED) is None


def test_quartet_below_the_cut_is_ignored():
    """The standard migrates late; an early copy must not be picked up."""
    early = [p - 1900 for p in GOOD]
    found = g.find_is_quartet(synth(early + GOOD), geometry=LEARNED)
    assert found is not None and found[0] == GOOD


def test_center_window_filters_to_the_expected_capillary():
    """A window around the learned centre keeps the search on one capillary."""
    decoy = [p + 600 for p in GOOD]
    found = g.find_is_quartet(synth(GOOD + decoy), geometry=LEARNED,
                              center=sum(GOOD) / 4.0, center_tol=80.0)
    assert found is not None and found[0] == GOOD


# --------------------------------------------------------------------------- #
# without a learned geometry: equimolar and clustered, nothing else
# --------------------------------------------------------------------------- #
def test_shape_free_search_finds_an_equimolar_cluster_of_any_spacing():
    """No spacing assumed -- which is what lets a new fragment be read at all."""
    for peaks in ([2100, 2130, 2215, 2250],       # d1=30  d2=85  d3=35
                  [2400, 2476, 2733, 2745],       # d1=76  d2=257 d3=12
                  [3000, 3012, 3400]):             # merged heteroduplexes
        found = g.find_is_quartet(synth(peaks))
        assert found is not None, peaks
        assert found[0] == peaks


def test_shape_free_search_never_returns_an_over_wide_span():
    """The ceiling is an invariant, not a filter on the input.

    Given four equimolar peaks the search may legitimately choose a
    three-band subset rather than all four, so the input alone does not decide.
    What must hold is that whatever comes back belongs to one fragment.
    """
    peaks = [2100, 2300, 2500, 2900]
    found = g.find_is_quartet(synth(peaks))
    assert found is not None
    assert found[0][-1] - found[0][0] <= g.IS_GROUP_SPAN


def test_shape_free_search_gives_up_when_no_three_bands_cluster():
    """Two tight pairs a long way apart are two fragments, not one standard."""
    assert g.find_is_quartet(synth([2100, 2200, 3400, 3500])) is None


def test_shape_free_search_needs_three_bands():
    assert g.find_is_quartet(synth([2100, 2130])) is None


def test_shape_free_search_returns_none_on_a_flat_trace():
    assert g.find_is_quartet(np.full(N_SCANS, BASELINE)) is None


# --------------------------------------------------------------------------- #
# PlateISModel.find: every stage runs, and the best match is what comes back
# --------------------------------------------------------------------------- #
# ``LEARNED`` says d1=30 with a tolerance of 8.  TRUE is a standard whose d1
# has drifted to 55 -- twice past what the fitted tolerance allows, so only the
# widest relaxed stage reaches it, and it sits inside the position window the
# marks establish.  GROUP matches the fitted shape exactly but is two hundred
# scans early and has no fourth band: exactly the kind of match the old
# first-hit rule reported, because the anchored stage searched without a window
# before any relaxed stage ever looked inside one.
TRUE = [2400, 2455, 2540, 2575]
GROUP = [1900, 1930, 2015]
# The same drift, from a model fitted from marks rather than handed LEARNED:
# d1=43 against a two-mark tolerance of 4.
DRIFT = [2400, 2443, 2528, 2563]
SEEDS = {"A01": [2400, 2430, 2515, 2550],
         "A02": [2404, 2434, 2519, 2554]}


def plate_model(centers=(2440.0, 2520.0)):
    return g.PlateISModel(LEARNED, list(centers), len(centers))


def test_find_reports_the_standard_inside_the_window_not_the_early_group():
    found = plate_model().find(synth(TRUE + GROUP, noise=4.0, seed=1))
    assert found is not None, "the standard must still be placed"
    assert list(found[0]) == TRUE
    assert found[2] == "relaxed-4"


def test_find_reports_the_group_when_it_is_all_the_well_has():
    found = plate_model().find(synth(GROUP, noise=4.0, seed=1))
    assert found is not None
    assert list(found[0]) == GROUP


def test_find_keeps_a_three_band_standard_when_the_fourth_band_is_absent():
    """A three-band read is half a match, not a rejected one."""
    found = plate_model().find(synth(TRUE[:3], noise=4.0, seed=1))
    assert found is not None
    assert list(found[0]) == TRUE[:3]


def test_find_gives_up_on_a_flat_trace():
    assert plate_model().find(np.full(N_SCANS, BASELINE)) is None


def _acgt(y):
    """One well's four columns, standard on column 3 (T), the rest flat."""
    acgt = np.full((N_SCANS, 4), BASELINE)
    acgt[:, g.acgt_index_for_channel("ACTG", 3)] = y
    return acgt


def test_semi_auto_plate_places_the_drifted_standard_not_the_early_group():
    """The whole workflow, and the failure that was reported against it: a
    well's real standard was ignored because an unwindowed three-band match
    came first in the stage order."""
    traces = {w: _acgt(synth(marks, noise=4.0, seed=i))
              for i, (w, marks) in enumerate(SEEDS.items())}
    traces["B08"] = _acgt(synth(DRIFT + GROUP, noise=4.0, seed=9))
    model, found = g.semi_auto_plate_model(
        traces, SEEDS, is_col=g.acgt_index_for_channel("ACTG", 3))
    assert model is not None
    assert list(found["B08"][0]) == DRIFT
    assert found["B08"][1] == "relaxed-4"
    assert found["B08"][2] in ("strong", "fair")


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


@needs_cache
def test_every_t9_well_has_an_equimolar_standard():
    """With the constants gone, the T9 plate is read on equimolarity alone."""
    with np.load(CACHE) as z:
        missing = [w for w in z.files
                   if g.find_is_quartet(z[w].astype(float)[:, 1]) is None]
    assert missing == []
