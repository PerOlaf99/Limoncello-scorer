"""A substantial satellite must not silently report as a homozygote.

Regression tests for the case where a peak one base downstream of a main is
swallowed as "+A".  _shoulders() has to pick one label and it picks "+A", so a
second real allele ends up outside clusters() and clust_frac() reports 0.0 -- a
heterozygote reported as a confident homozygote.

These tests pin the *mechanism* (ambiguity is recorded and propagates to an
uncertain call) rather than a particular threshold, because the threshold is
uncalibrated and is expected to move once known minor-allele samples exist.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import genotyping as g
import scorer as s

N_SCANS = 4000


class _Doc:
    """Minimal stand-in for analyzer_core.TraceDocument."""

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


def _trace(second_frac, apex=2000, spacing=14, sigma=3.0, n_scans=N_SCANS):
    """Channel 0 with a main peak and a same-channel peak one spacing later.

    *second_frac* is the second peak's height relative to the main.
    """
    acgt = np.full((n_scans, 4), 50.0)
    x = np.arange(n_scans)
    for a, h in ((apex, 1000.0), (apex + spacing, 1000.0 * second_frac)):
        acgt[:, 0] += h * np.exp(-0.5 * ((x - a) / sigma) ** 2)
    return acgt


def _picker(second_frac, **kw):
    acgt = _trace(second_frac, **kw)
    p = g.PeakPicker(_Doc(acgt), Path("/tmp/synthetic.rsd"))
    p.pick(2000)
    return p


# --------------------------------------------------------------------------- #
# the picker records the ambiguity
# --------------------------------------------------------------------------- #
def test_major_secondary_peak_is_recorded_as_ambiguous():
    p = _picker(0.5)
    sats = p.ambiguous_satellites()
    assert len(sats) == 1
    assert sats[0]["kind"] == "+A"
    assert sats[0]["ambiguous"] is True


def test_tiny_plus_a_tail_is_not_flagged():
    """A normal +A product stays quiet; otherwise the flag is noise."""
    p = _picker(0.05)
    assert p.ambiguous_satellites() == []


def test_ambiguity_reports_the_measured_height_fraction():
    """Overlap compresses a real 50% peak well below 0.5 -- the flag must not
    assume the picked height is the true allele ratio."""
    p = _picker(0.5)
    frac = p.ambiguous_satellites()[0]["height_frac"]
    assert 0.0 < frac < 0.5


def test_export_rows_carries_the_ambiguity():
    p = _picker(0.5)
    rows = p.export_rows()
    sat = [r for r in rows if r["kind"] == "+A"][0]
    assert sat["ambiguous"] is True
    assert sat["height_frac"] == p.ambiguous_satellites()[0]["height_frac"]
    main = [r for r in rows if r["kind"] == "main"][0]
    assert main["ambiguous"] is False


# --------------------------------------------------------------------------- #
# the bug itself: the peak is dropped and the fraction reads 0.0
# --------------------------------------------------------------------------- #
def test_unflagged_would_have_been_a_silent_homozygote():
    """Document why this matters: the raw picker state is a confident 0.0."""
    p = _picker(0.5)
    mains = [r for r in p.records if r["kind"] == "main"]
    assert [len(c) for c in p.clusters()] == [1]
    assert [p.clust_frac(r) for r in mains] == [0.0]


# --------------------------------------------------------------------------- #
# it propagates to the scorer
# --------------------------------------------------------------------------- #
# Note on the scorer's model: group_alleles() never mixes channels, and
# zygosity comes from the CTC duplex labels (any HET peak proves a het) rather
# than from pairing alleles across channels.  So the damaging case is a
# homozygote carrying a large unresolved same-channel peak one base later: the
# cluster is a lone HOM peak and scores as a confident "hom-major" with
# fraction 0.0.


def _hom_rows(second_frac=0.5, **kw):
    """Rows for a HOM duplex with an unresolved same-channel peak behind it."""
    rows = [dict(r) for r in _picker(second_frac, **kw).export_rows()]
    for r in rows:
        if r["kind"] == "main":
            r["duplex"] = "HOM1"
    return rows


def test_scorer_downgrades_to_uncertain_instead_of_homozygote():
    """The regression: a lone HOM peak scoring confidently at 0.0."""
    rows = _hom_rows()
    clusters = s.group_alleles(rows, s.median_peak_spacing(rows))
    # The raw cluster really does look like a clean homozygote...
    assert s.genotype_call(clusters[0]) == "hom-major"
    assert s.minor_fraction(clusters[0]) == 0.0
    # ...and that is exactly what used to be reported as fact.
    scored = s.score_rows(rows)
    main = [r for r in scored if r["kind"] == "main"][0]
    assert main["call"] == "uncertain"
    assert main["confidence"] == 0.0
    assert "minor allele" in main["ambiguous_reason"]


def test_quiet_plus_a_still_scores_confidently():
    scored = s.score_rows(_hom_rows(second_frac=0.05))
    mains = [r for r in scored if r["kind"] == "main"]
    assert all(r["call"] == "hom-major" for r in mains)
    assert not any("ambiguous_reason" in r for r in mains)


def test_well_summary_reports_the_uncertainty():
    rows = [dict(r, file="synthetic.rsd", well="A01") for r in _hom_rows()]
    w = s.well_summary(rows)[0]
    assert w["n_uncertain"] == 1
    assert w["genotype"] == "uncertain"


def test_a_clean_homozygote_still_reports_hom():
    """An ordinary hom well must not start reporting uncertain."""
    rows = [dict(r, file="synthetic.rsd", well="A01")
            for r in _picker(0.0).export_rows()]
    for r in rows:
        if r["kind"] == "main":
            r["duplex"] = "HOM1"
    w = s.well_summary(rows)[0]
    assert w["n_uncertain"] == 0
    assert w["genotype"] == "hom-major"


def test_helper_detects_ambiguity():
    rows = _hom_rows()
    sat = [r for r in rows if r["kind"] == "+A"][0]
    assert s.has_ambiguous_satellite([sat]) is True
    assert s.has_ambiguous_satellite([r for r in rows if r["kind"] == "main"]) is False
