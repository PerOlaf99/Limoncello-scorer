"""rs1695 / T9 calls against the 96-well manual score.

The plate is a real CTC-CE run (``OY_rs1695_T9_270910Run01``), so this is the
project's only ground truth.  The traces themselves are far too large to keep in
the repository, so the fixtures hold the *measured* per-position duplex areas and
peak significances instead; the calling rule is what gets regression-tested.

    rs1695_expected.csv   the manual score (well, expected call, notes)
    rs1695_measured.csv   measured hom1/hom2/het1/het2 areas + sigmas per well

One disagreement is expected and documented: G09.
"""
import csv
from pathlib import Path

import pytest

import scorer as s

DATA = Path(__file__).parent / "data"


def _load(name):
    with open(DATA / name, newline="") as fh:
        return {r["well"]: r for r in csv.DictReader(fh)}


EXPECTED = _load("rs1695_expected.csv")
MEASURED = _load("rs1695_measured.csv")

# G09 is scored "Het (AI) Big H2" and marked Borderline.  It is called hom-2 here
# on purpose: its heteroduplex (18.5/28.8 sigma) sits at the detection floor and
# C07 outranks it on both measures while being a real homo2, so no threshold
# separates them.  Relaxing T9_MIN_HET_SIGMA to catch it costs 7 true homozygotes.
KNOWN_DIVERGENCE = {"G09"}


def _score(well):
    m = MEASURED[well]
    areas = [float(m[f"hom{i}"]) for i in (1, 2)]
    hets = [float(m[f"het{i}"]) for i in (1, 2)]
    sig = [float(m[f"snr{i}"]) for i in (1, 2, 3, 4)]
    return s.t9_call(areas[0], areas[1], hets[0], hets[1], sig)


@pytest.fixture(scope="module")
def calls():
    return {w: _score(w) for w in MEASURED}


def test_fixtures_cover_the_whole_plate():
    assert len(EXPECTED) == 96
    assert set(EXPECTED) == set(MEASURED)


def test_calls_match_manual_score(calls):
    wrong = {}
    for well, exp in EXPECTED.items():
        want = exp["expected"]
        if want == "none":
            want = "no-call"
        elif want == "hom1":
            want = "hom-1"
        elif want == "hom2":
            want = "hom-2"
        got = calls[well][0]
        if got != want and well not in KNOWN_DIVERGENCE:
            wrong[well] = (want, got)
    assert not wrong, f"unexplained call differences: {wrong}"


def test_only_known_divergences(calls):
    differing = {w for w, e in EXPECTED.items()
                 if calls[w][0] != {"none": "no-call", "hom1": "hom-1",
                                     "hom2": "hom-2"}.get(e["expected"], e["expected"])}
    assert differing == KNOWN_DIVERGENCE


def test_heteroduplex_beats_homoduplex_ratio():
    """A het with severe allelic imbalance must not read as a homozygote.

    A12 is scored het; its homoduplex ratio alone (89860 vs 85700) says
    homozygous, and only the heteroduplex pair gives it away.
    """
    call, frac, flags = _score("A12")
    assert call == "het"
    assert "ai" in flags
    assert frac < s.T9_AI_DEVIATION          # the imbalance that makes it AI


def test_homozygote_is_not_read_as_het_by_a_one_sided_tail():
    """B12 is homo2 with a 119/24 sigma one-sided heteroduplex.

    An asymmetric tail is an artifact, not evidence of heterozygosity.
    """
    call, _frac, _flags = _score("B12")
    assert call == "hom-2"


def test_ai_flag_only_on_imbalanced_hets(calls):
    flagged = {w for w, c in calls.items() if "ai" in c[2]}
    # G09 is scored AI but is called homo-2 here, so it is the only omission.
    assert flagged == {"A12", "G10"}


def test_no_product_well_is_not_called(calls):
    call, frac, _flags = calls["H01"]
    assert call == "no-call"
    assert frac == 0.0


def test_allele_fraction_matches_reported_value():
    m = MEASURED["A01"]
    got = s.t9_allele_fraction(*(float(m[k]) for k in
                                 ("hom1", "hom2", "het1", "het2")))
    assert got == pytest.approx(0.508, abs=0.002)


def test_allele_fraction_is_zero_without_sample():
    assert s.t9_allele_fraction(0, 0, 0, 0) == 0.0


def test_area_only_path_needs_both_heteroduplexes():
    """With no significances, area comparability stands in for the sigma test."""
    # Areas are past T9_MIN_TOTAL_AREA so the no-call gate does not fire first.
    assert s.t9_call(2000, 2000, 1500, 1500)[0] == "het"   # both, comparable
    assert s.t9_call(2000, 2000, 1500, 0)[0] == "hom-1"    # one-sided -> not het


def test_symmetric_heterozygote_is_exactly_one_half():
    """With hom1 == hom2 and both heteroduplexes present the fraction is 0.5.

    Each heteroduplex is half A, so it cannot bias the fraction -- only an
    unbalanced homoduplex can, which is what makes AI detectable.
    """
    _call, frac, flags = s.t9_call(1000, 1000, 200, 200, (900, 900, 300, 300))
    assert frac == pytest.approx(0.5)
    assert "ai" not in flags


def test_ai_flag_needs_an_unbalanced_homoduplex():
    # Almost all H2 homoduplex, no H1, but both heteroduplexes clearly present:
    # the "Het, big H2" pattern.
    call, frac, flags = s.t9_call(0, 1000, 200, 200, (50, 5000, 300, 300))
    assert call == "het"
    assert 0.0 < frac < s.T9_AI_DEVIATION
    assert "ai" in flags
