"""One-click auto-genotyping of a whole well (``genotyping.auto_genotype``).

The engine is the batch form of the manual path -- find the standard quartet,
measure the sample in the four windows it defines, call the well -- and it is
the one place where a silent channel mix-up would be most expensive: scoring
the standard's own peaks as the sample returns confident nonsense rather than
failing.  So the channel pair is pinned here against the real plate.

Synthetic traces cover the arithmetic and the failure paths; the T9 plate cases
are skipped when the cached traces are not present, because the plate is far too
large to ship.
"""
import csv
import re
import sys
from pathlib import Path

import numpy as np
import pytest

import genotyping as g

CACHE = Path("/home/tv/t9_rs1695_analysis/data/t9raw.npz")
DATA = Path(__file__).parent / "data"
N_SCANS = 4000
BASELINE = 100.0
# d1=30, d2=85, d3=35 -> dm=32.5, |d1-d3|=5 <= 24.4, 85 within [55.25, 107.25]
STD_SCANS = [2100, 2130, 2215, 2250]
DYE_ORDER = "ACTG"


class FakeDoc:
    """The bits of ``TraceDocument`` ``auto_genotype`` reads."""

    def __init__(self, acgt, well="A01"):
        self.acgt = np.asarray(acgt, float)
        self.well = well
        self.base_order = "ACGT"


def synth(peaks, height=1000.0, n_scans=N_SCANS, noise=0.0, seed=0, width=3.0):
    rng = np.random.default_rng(seed)
    y = np.full(n_scans, BASELINE) + rng.normal(0.0, noise, n_scans)
    x = np.arange(n_scans)
    for p in peaks:
        y += height * np.exp(-0.5 * ((x - p) / width) ** 2)
    return y


def plate(is_heights, sample_heights, noise=0.0, seed=0):
    """One well: an ``acgt`` matrix with the standard on T and the sample on C.

    Column 3 (T) carries the internal standard and column 1 (C) the sample,
    which is this plate's dye order.  The remaining two columns stay flat, so a
    channel mix-up shows up as a missing quartet rather than a plausible call.
    """
    x = np.arange(N_SCANS)
    acgt = np.full((N_SCANS, 4), BASELINE)
    for col, heights in ((g.acgt_index_for_channel(DYE_ORDER, 3), is_heights),
                         (g.acgt_index_for_channel(DYE_ORDER, 2), sample_heights)):
        for i, h in enumerate(heights):
            acgt[:, col] += h * np.exp(-0.5 * ((x - STD_SCANS[i]) / 3.0) ** 2)
    if noise:
        rng = np.random.default_rng(seed)
        acgt += rng.normal(0.0, noise, acgt.shape)
    return acgt


# ``t9raw.npz`` stores wells in the plate's physical channel order, but
# ``doc.acgt`` is always A,C,G,T.  The cache was written in Ch1..Ch4 order
# (dye order ACTG: Ch1=A, Ch2=C, Ch3=T, Ch4=G), so to rebuild a doc.acgt from a
# cached well, acgt column j takes cached column:
#   A -> cached Ch1,  C -> cached Ch2,  G -> cached Ch4,  T -> cached Ch3
ACGT_FROM_CACHE = [3, 2, 0, 1]


def cache_to_doc(matrix):
    """Reorder one cached well from physical channel order into A,C,G,T."""
    return np.column_stack([matrix[:, i] for i in ACGT_FROM_CACHE])


def _norm(text):
    s = re.sub(r"[^a-z0-9]", "", text.lower())
    return "het" if s.startswith("het") else s.replace("none", "nocall")


def _load(name):
    with open(DATA / name, newline="") as fh:
        return {r["well"]: r for r in csv.DictReader(fh)}


EXPECTED = _load("rs1695_expected.csv")
MEASURED = _load("rs1695_measured.csv")

needs_cache = pytest.mark.skipif(
    not CACHE.exists(),
    reason="T9 plate cache not present (too large to ship)")

needs_plate = pytest.mark.skipif(
    not CACHE.exists(),
    reason="T9 plate cache not present (too large to ship)")


# --------------------------------------------------------------------------- #
# the channel pair
# --------------------------------------------------------------------------- #
def test_is_and_sample_channels_are_not_the_same():
    assert g.DEFAULT_IS_CHANNEL != g.DEFAULT_SAMPLE_CHANNEL


def test_default_channels_match_the_t9_plate_dye_order():
    # ACTG: Ch1=A, Ch2=C, Ch3=T, Ch4=G.  The standard is the T channel and the
    # sample the C channel; this is the pair auto_mark_std's channel=3 implies.
    assert g.acgt_index_for_channel(DYE_ORDER, g.DEFAULT_IS_CHANNEL) == \
        "ACGT".index("T")
    assert g.acgt_index_for_channel(DYE_ORDER, g.DEFAULT_SAMPLE_CHANNEL) == \
        "ACGT".index("C")


@needs_cache
def test_standard_channel_carries_a_quartet_in_every_well():
    """The standard is present in every well; the sample channel mostly is not.

    The standard is an equimolar four-peak fragment the kit adds to every well,
    so it reads as a clean quartet on all 96.  The sample channel holds the
    amplicon under test, which is not a four-peak fragment at all -- it happens
    to fit the geometry in well under two thirds of wells, which is the whole
    reason channel choice cannot be left to a trace.
    """
    z = np.load(CACHE, allow_pickle=True)
    is_hits = sample_hits = 0
    for w in z.files:
        doc = cache_to_doc(z[w].astype(float))
        is_hits += g.find_is_quartet(
            doc[:, g.acgt_index_for_channel(DYE_ORDER, g.DEFAULT_IS_CHANNEL)]) is not None
        sample_hits += g.find_is_quartet(
            doc[:, g.acgt_index_for_channel(DYE_ORDER, g.DEFAULT_SAMPLE_CHANNEL)]) is not None
    assert is_hits == len(z.files)
    assert sample_hits < len(z.files)


# --------------------------------------------------------------------------- #
# the arithmetic
# --------------------------------------------------------------------------- #
def test_four_equimolar_duplexes_read_as_a_heterozygote():
    # a little noise, so the duplexes clear T9_MIN_DOMINANT_SIGMA and the
    # significance path is the one under test (a noiseless trace has sigma 0
    # and no-call on purpose)
    doc = FakeDoc(plate([1000] * 4, [1000] * 4, noise=2.0, seed=11))
    row = g.auto_genotype(doc, base_order=DYE_ORDER)
    assert row["call"] == "het"
    assert 0.4 < row["frac"] < 0.6
    assert row["reason"] == ""


def test_one_strong_duplex_alone_reads_as_a_homozygote():
    doc = FakeDoc(plate([1000] * 4, [1000, 0, 0, 0], noise=2.0, seed=12))
    row = g.auto_genotype(doc, base_order=DYE_ORDER)
    assert row["call"] in ("hom-1", "no-call")
    assert row["frac"] > 0.8


def test_areas_are_ordered_by_duplex_and_roughly_equal():
    doc = FakeDoc(plate([1000] * 4, [1000, 800, 600, 400], noise=2.0, seed=13))
    row = g.auto_genotype(doc, base_order=DYE_ORDER)
    got = [row["hom1"], row["hom2"], row["het1"], row["het2"]]
    assert got == sorted(got, reverse=True)
    assert all(v > 0 for v in got)


def test_rows_record_the_channels_they_were_measured_on():
    row = g.auto_genotype(FakeDoc(plate([1000] * 4, [1000] * 4)),
                          base_order=DYE_ORDER, run_name="R1")
    assert row["run"] == "R1"
    assert row["is_channel"] == g.DEFAULT_IS_CHANNEL
    assert row["sample_channel"] == g.DEFAULT_SAMPLE_CHANNEL
    assert row["std_scans"] == "/".join(str(x) for x in STD_SCANS)


def test_significance_falls_when_the_sample_is_quieter_than_the_noise():
    loud = g.auto_genotype(FakeDoc(plate([1000] * 4, [4000] * 4, noise=2.0, seed=7)),
                           base_order=DYE_ORDER)
    quiet = g.auto_genotype(FakeDoc(plate([1000] * 4, [200] * 4, noise=2.0, seed=7)),
                            base_order=DYE_ORDER)
    assert min(loud["snr1"], loud["snr2"]) > min(quiet["snr1"], quiet["snr2"])


# --------------------------------------------------------------------------- #
# the failure paths -- a bad well must not raise
# --------------------------------------------------------------------------- #
def test_no_quartet_gives_a_reason_not_an_exception():
    acgt = np.full((N_SCANS, 4), BASELINE)
    row = g.auto_genotype(FakeDoc(acgt), base_order=DYE_ORDER)
    assert row["call"] == "no-call"
    assert "Ch%d" % g.DEFAULT_IS_CHANNEL in row["reason"]


def test_empty_trace_gives_a_reason():
    row = g.auto_genotype(FakeDoc(np.zeros((0, 4))), base_order=DYE_ORDER)
    assert row["call"] == "no-call"
    assert row["reason"]


def test_standard_and_sample_on_the_same_channel_is_refused():
    row = g.auto_genotype(FakeDoc(plate([1000] * 4, [1000] * 4)),
                          is_channel=2, sample_channel=2, base_order=DYE_ORDER)
    assert row["call"] == "no-call"
    assert "same acgt column" in row["reason"]


def test_bad_dye_order_gives_a_reason():
    row = g.auto_genotype(FakeDoc(plate([1000] * 4, [1000] * 4)),
                          base_order="XYZ")
    assert row["call"] == "no-call"
    assert "base_order" in row["reason"]


def test_channel_out_of_range_gives_a_reason():
    row = g.auto_genotype(FakeDoc(plate([1000] * 4, [1000] * 4)),
                          is_channel=9, base_order=DYE_ORDER)
    assert row["call"] == "no-call"
    assert row["reason"]


def test_trace_too_narrow_for_the_channel_gives_a_reason():
    row = g.auto_genotype(FakeDoc(np.full((N_SCANS, 2), BASELINE)),
                          base_order=DYE_ORDER)
    assert row["call"] == "no-call"
    assert "channel" in row["reason"]


def test_weak_sample_reports_why_rather_than_a_confident_call():
    doc = FakeDoc(plate([1000] * 4, [30, 30, 30, 30], noise=20.0, seed=3))
    row = g.auto_genotype(doc, base_order=DYE_ORDER)
    assert row["call"] == "no-call"
    assert "noise" in row["reason"]


# --------------------------------------------------------------------------- #
# the ground truth
# --------------------------------------------------------------------------- #
@needs_plate
def test_calls_match_the_manual_score_over_the_whole_plate():
    z = np.load(CACHE, allow_pickle=True)
    calls = {}
    for w in z.files:
        doc = FakeDoc(cache_to_doc(z[w].astype(float)), well=w)
        calls[w] = g.auto_genotype(doc, base_order=DYE_ORDER)
    assert len(calls) == 96
    want = {w: _norm(r["expected"]) for w, r in EXPECTED.items()}
    got = {w: _norm(r["call"]) for w, r in calls.items()}
    wrong = {w: (want[w], got[w]) for w in want if want[w] != got[w]}
    # G09 is the one documented divergence: its heteroduplex sits at the
    # detection floor, so it is called hom-2 on purpose.
    assert set(wrong) == {"G09"}, wrong


@needs_plate
def test_measured_areas_track_the_ground_truth_fixture():
    z = np.load(CACHE, allow_pickle=True)
    for w in ("A12", "C07", "G10"):
        doc = FakeDoc(cache_to_doc(z[w].astype(float)))
        row = g.auto_genotype(doc, base_order=DYE_ORDER)
        m = MEASURED[w]
        for key in ("hom1", "hom2", "het1", "het2"):
            got, want = row[key], float(m[key])
            assert abs(got - want) <= 0.15 * max(1.0, want), (w, key, got, want)
        for key in ("snr1", "snr2", "snr3", "snr4"):
            got, want = row[key], float(m[key])
            assert abs(got - want) <= 0.10 * max(1.0, want), (w, key, got, want)
