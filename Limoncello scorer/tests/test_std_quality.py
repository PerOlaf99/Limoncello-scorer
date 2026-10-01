"""A weak standard is flagged, not rejected; a half-formed one is rejected.

Two different things can go wrong with the internal standard, and they are not
the same problem:

* **weak** -- the pattern is there but faint.  The operator reads these in
  context, judging a peak against the wells above and below it, so the call
  stands and only a flag is raised.  ABCC2 D07 reads 18 sigma against a plate
  median of 112 and is still an unambiguous hom-1.
* **half-formed** -- the heteroduplexes fail to resolve on a plate whose marked
  wells always resolve them.  That is not a weak standard, it is a different
  shape, and the sample windows read off it are misplaced, so the well is not
  genotypeable at all.  Every one of the operator's 56 ABCC2_N10 marks is
  four-band, and its only two three-band reads (D01, E01) are both no-calls.

An absolute signal cutoff cannot tell these apart: D01 measures 27 sigma and
D07 18, yet D01 is unusable and D07 is readable.  Judging against the plate's
own marked wells can.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import genotyping as g

N_SCANS = 4000
# ACTG: Ch1=A(col 0), Ch2=C(col 1), Ch3=T(col 3), Ch4=G(col 2).
# The standard rides Ch3 and the sample Ch2, as on the real plates.
STD_COL = g.acgt_index_for_channel("ACTG", 3)
SAMP_COL = g.acgt_index_for_channel("ACTG", 2)


class _Doc:
    def __init__(self, acgt, base_order="ACTG", well="A01"):
        self.acgt = acgt
        self.n_scans = acgt.shape[0]
        self.base_order = base_order
        self.path = Path("/tmp/synthetic.rsd")
        self.well = well
        self.sequence = ""
        self.qualities = []
        self.peak_positions = []
        self.raw = acgt


def _acgt(std, samp=(), n_scans=N_SCANS, noise=4.0, std_amp=900.0,
          samp_amp=700.0, seed=0):
    """Standard duplexes on Ch3, sample duplexes on Ch2, over a noise floor."""
    a = np.full((n_scans, 4), 50.0)
    x = np.arange(n_scans)
    if noise:
        a = a + np.random.default_rng(seed).normal(0.0, noise, a.shape)
    for p in std:
        a[:, STD_COL] += std_amp * np.exp(-0.5 * ((x - p) / 4.0) ** 2)
    for p in samp:
        a[:, SAMP_COL] += samp_amp * np.exp(-0.5 * ((x - p) / 4.0) ** 2)
    return a


# Plate geometry fitted from four-band marks: d1=80, d2=260, d3=12.
SOLVE = g.ISGeometry(80.0, 260.0, 12.0, 25.0, 78.0, 9.0)
FOUR = [2400, 2480, 2740, 2752]
THREE = [2400, 2480, 2740]


def _model(merged_seen, std_snr=None):
    m = g.PlateISModel(SOLVE, [sum(FOUR) / 4.0], n_marked=20,
                       merged_seen=merged_seen)
    m.std_snr = list(std_snr or [])
    return m


# --------------------------------------------------------------------------- #
# what a merged detection means depends on the plate
# --------------------------------------------------------------------------- #
def test_merged_geometry_records_whether_marks_showed_one():
    only_four = g.PlateISModel.from_marks([FOUR] * 5)
    assert only_four.merged_seen is False
    with_merged = g.PlateISModel.from_marks([FOUR] * 5 + [THREE] * 2)
    assert with_merged.merged_seen is True


def test_three_band_standard_is_off_model_on_a_never_merged_plate():
    doc = _Doc(_acgt(THREE, samp=[2400, 2480]))
    row = g.auto_genotype(doc, is_channel=3, sample_channel=2,
                          base_order="ACTG", is_model=_model(False))
    assert row["call"] == "no-call"
    assert "unresolved heteroduplex" in row["reason"]
    assert row["het_resolved"] is False


def test_three_band_standard_is_normal_on_a_merged_plate():
    doc = _Doc(_acgt(THREE, samp=[2400, 2480, 2740]))
    row = g.auto_genotype(doc, is_channel=3, sample_channel=2,
                          base_order="ACTG", is_model=_model(True))
    assert row["call"] != "no-call", row["reason"]
    assert "het-merged" in row["flags"]


def test_manual_marks_override_the_off_model_rejection():
    """The operator's word stands where the shape test would refuse."""
    doc = _Doc(_acgt(THREE, samp=[2400, 2480]))
    row = g.auto_genotype(doc, is_channel=3, sample_channel=2, base_order="ACTG",
                          std_scans_manual=THREE, is_model=_model(False))
    assert "unresolved heteroduplex" not in row["reason"]


def test_four_band_standard_is_unaffected_by_the_shape_test():
    doc = _Doc(_acgt(FOUR, samp=[2400, 2480]))
    row = g.auto_genotype(doc, is_channel=3, sample_channel=2,
                          base_order="ACTG", is_model=_model(False))
    assert "unresolved heteroduplex" not in row["reason"]
    assert row["het_resolved"] is True


# --------------------------------------------------------------------------- #
# weakness is judged against the plate
# --------------------------------------------------------------------------- #
def test_weak_standard_is_flagged_and_the_call_stands():
    strong = _acgt(FOUR, samp=[2400, 2480])
    faint = _acgt(FOUR, samp=[2400, 2480], std_amp=60.0)
    model = _model(False, std_snr=[400.0, 420.0, 390.0, 410.0])
    ok = g.auto_genotype(_Doc(strong), is_channel=3, sample_channel=2,
                         base_order="ACTG", is_model=model)
    faint_row = g.auto_genotype(_Doc(faint), is_channel=3, sample_channel=2,
                                base_order="ACTG", is_model=model)
    assert "std-weak" not in ok["flags"]
    assert "std-weak" in faint_row["flags"]
    # The call itself is unchanged: a faint pattern still reads in context.
    assert faint_row["call"] == ok["call"]


def test_threshold_needs_a_reference_before_it_flags():
    """With too few marked wells there is no plate to judge against."""
    assert _model(False, std_snr=[100.0, 200.0]).weak_std_threshold() is None
    assert _model(False, std_snr=[100.0] * 5).weak_std_threshold() == 25.0


def test_manual_wells_are_still_flagged():
    """Positions are trusted; the measured strength is still reported."""
    faint = _acgt(FOUR, samp=[2400, 2480], std_amp=60.0)
    row = g.auto_genotype(_Doc(faint, well="D07"), is_channel=3, sample_channel=2,
                          base_order="ACTG", std_scans_manual=FOUR,
                          is_model=_model(False, std_snr=[400.0] * 5))
    assert "std-weak" in row["flags"]
    assert row["call"] == "hom-1"


def test_weakest_band_snr_reads_the_standard_channel():
    a = _acgt(FOUR)
    strong = g.weakest_std_band_snr(a, STD_COL, FOUR)
    assert strong > 100.0
    # Nothing at the standard positions -> no signal above the noise.
    assert g.weakest_std_band_snr(_acgt([]), STD_COL, FOUR) < strong