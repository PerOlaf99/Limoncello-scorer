"""Genotyping rests on H1 + H2; heteroduplex resolution must not gate it.

The standard exists to establish the nature of the two homoduplexes.  Whether
the mismatched heteroduplexes separate is a different question, and their
answer does not change the call.  These tests pin that separation of concerns:

* a standard that merges its heteroduplexes is a complete standard, not a
  degraded one, so it must still produce a call;
* a merged-het well is still a heterozygote when both alleles are present;
* HET2, when present, is attached as an extra measurement rather than being
  required for the match;
* the flag distinguishes the two so an allelic-imbalance read knows when the
  HET1/HET2 areas it would compare are not actually separable.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import genotyping as g
import scorer as s

N_SCANS = 4000
# Areas must clear s.T9_MIN_TOTAL_AREA (3000) to be scored at all.
HOM1, HOM2, HET = 4000.0, 4000.0, 3600.0


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


def _channel(peaks, sigma=4.0, n_scans=N_SCANS):
    """One channel carrying equal-height Gaussian duplexes at *peaks*."""
    y = np.full(n_scans, 50.0)
    x = np.arange(n_scans)
    for p in peaks:
        y += 900.0 * np.exp(-0.5 * ((x - p) / sigma) ** 2)
    return y


# Geometry fitted from a plate whose marks were all resolved four-band wells.
SPLIT = g.ISGeometry(80.0, 260.0, 12.0, 25.0, 78.0, 9.0)


# --------------------------------------------------------------------------- #
# the shape, not the split
# --------------------------------------------------------------------------- #
def test_merged_standard_is_a_complete_match():
    """Three bands satisfy the template outright -- not 'close enough'."""
    assert SPLIT.matches(2400, 2480, 2740, None)
    assert not SPLIT.merged
    assert SPLIT.matches_hom(2400, 2480)
    assert not SPLIT.matches_hom(2400, 2600)


def test_h1_h2_match_independently_of_any_heteroduplex():
    """The minimal test a call depends on needs nothing past H2."""
    assert SPLIT.matches_hom(2400, 2480)
    assert SPLIT.matches_hom(2400, 2472)          # 72, inside tol1 = 25
    assert not SPLIT.matches_hom(2400, 2540)      # 140, outside tol1


def test_het2_is_attached_only_when_it_really_resolves():
    """HET1 -> HET2 = 12 scans, the fitted d3: report all four."""
    y = _channel([2400, 2480, 2740, 2752])
    pk = g._is_candidates(y, 1900, geometry=SPLIT)
    scans, _heights = g._search_core(pk, SPLIT)
    assert scans == [2400, 2480, 2740, 2752]


def test_unresolved_heteroduplex_still_yields_the_core():
    """No fourth band: the three-band core is returned, unharmed."""
    y = _channel([2400, 2480, 2740])
    pk = g._is_candidates(y, 1900, geometry=SPLIT)
    scans, _heights = g._search_core(pk, SPLIT)
    assert scans == [2400, 2480, 2740]


def test_core_prefers_the_band_that_admits_a_partner():
    """HET2 alone can satisfy the core, so the score must include the partner.

    Regression: scoring the bare core took whichever of HET1/HET2 was taller,
    which then searched for a fifth band and dropped a resolved pair back to
    three.  ABCC2 fell from 56/56 resolved standards to 0/56 that way.
    """
    y = np.full(N_SCANS, 50.0)
    x = np.arange(N_SCANS)
    for p, h in ((2400, 900.0), (2480, 900.0), (2740, 820.0), (2752, 980.0)):
        y += h * np.exp(-0.5 * ((x - p) / 4.0) ** 2)
    pk = g._is_candidates(y, 1900, geometry=SPLIT)
    scans, _heights = g._search_core(pk, SPLIT)
    assert scans == [2400, 2480, 2740, 2752], scans


# --------------------------------------------------------------------------- #
# the call
# --------------------------------------------------------------------------- #
def test_merged_het_allele_pair_is_called_het():
    """Both alleles in the merged band -> het, not a fallback to hom-1."""
    call, frac, _flags = s.t9_call(HOM1, HOM2, HET, None)
    assert call == "het"
    assert 0.4 < frac < 0.6, frac


def test_merged_het_absent_is_hom():
    call, _frac, _flags = s.t9_call(HOM1, HOM2, 0.0, None)
    assert call == "hom-1"
    call, _frac, _flags = s.t9_call(0.0, HOM2, 0.0, None)
    assert call == "hom-2"


def test_merged_het_below_noise_is_not_het():
    """The relaxed test still needs the band above the significance floor."""
    call, _frac, _flags = s.t9_call(HOM1, HOM2, HET, None,
                                    sigmas=[50.0, 50.0, 0.1, 0.0])
    assert call != "het"


def test_split_sample_test_is_unchanged():
    """One band where two heteroduplexes belong is a shoulder, not a het."""
    assert s.t9_call(HOM1, HOM2, HET, HET)[0] == "het"
    assert s.t9_call(HOM1, HOM2, HET, 0.0)[0] == "hom-1"
    # Strongly asymmetric pair: still not a clean heterozygote.
    assert s.t9_call(HOM1, HOM2, HET, 20.0)[0] == "hom-1"


def test_empty_merged_sample_is_no_call():
    assert s.t9_call(0.0, 0.0, 0.0, None)[0] == "no-call"


def test_no_product_still_reads_hom_1():
    """The user's rule: no product is hom-1 whether the het merged or not."""
    assert s.t9_call(HOM1, HOM2, 0.0, None)[0] == "hom-1"
    assert s.t9_call(HOM1, HOM2, 0.0, 0.0)[0] == "hom-1"