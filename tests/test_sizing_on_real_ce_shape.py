"""Fragment sizing on a ladder with MegaBACE-measured trace shape.

``tests/test_fragment_sizing.py`` sizes a clean synthetic trace: no noise, one
peak width, no baseline.  That is enough to pin the algebra, and it would pass
just as happily on an instrument with a different peak shape.  What it cannot
tell us is whether the sizer survives the *shape* of a real MegaBACE run.

So the ladder and the sample here are built from statistics measured off the
machine's own M13 traces rather than invented.  Every constant below is a
measurement, with the run it came from named alongside:

* peak FWHM 10-12 scans (median 11), so consecutive bands are narrow enough to
  overlap wherever the migration packs them close together;
* baseline noise SD 6-21 per point in the injection region (median ~8), against
  peak heights of 800-2100 -- a signal-to-noise ratio of 100-500x, which is what
  sets how many spurious peaks the detector sees;
* an injection baseline of 170-350 rather than zero;
* a migration slope of ~824 scans per ln(bp), fitted to the base-call positions
  of three real reads.

The slope matters as much as the width.  At 824 scans/ln(bp) a GeneScan 500
ladder spreads over ~2180 scans, so neighbouring fragments sit 145 scans apart
on average -- but only 16 scans apart at its tightest (139/150 bp), and
Geneflo 1000 closes to 21.  Those are the pairs that decide whether the
alignment really skips a band or just shifts every length.

This is still a synthetic ladder: it tests the sizing code against real
instrument shape, not real ladder chemistry.  The accuracy figure it reports is
a property of the model used here.  Only a run with a real size standard in its
own channel gives the accuracy of the sizer itself -- read the leave-one-out RMS
off that run, never the fit residual.
"""
import numpy as np
import pytest

from analyzer_core import acgt_index_for_channel
import fragment_sizing as fs

# --------------------------------------------------------------------------- #
# measured off ~/Dokumenter/MB4000_DEMO_DATA (four wells, four channels)
# --------------------------------------------------------------------------- #
FWHM_SCANS = 11.0          # median full width at half maximum of a real peak
NOISE_SD = 8.0             # per-point baseline noise SD (injection region)
PEAK_HEIGHT = 1300.0       # median height of a called base peak
BASELINE = 180.0           # median injection-region baseline level
BASE_ORDER = "ACTG"
LADDER_CHANNEL = 4         # ROX
SAMPLE_CHANNEL = 2         # FAM

N_SCANS = 11000

# Migration slope in scans per ln(bp), fitted to the base-call positions of the
# machine's own M13 reads (801, 851 and 824 scans/ln(bp) on A01, A02 and C03).
# Taken from real reads rather than chosen: it sets how far apart consecutive
# ladder bands land, and that gap is what decides whether the alignment has to
# work for its answer.
MIGRATION_SLOPE = 824.0


def scan_of(bp, *, offset=1200.0):
    """The CE migration: scan grows with log(bp)."""
    return offset + MIGRATION_SLOPE * np.log(bp)


def _band(trace, x, scan, height, fwhm=FWHM_SCANS):
    """Add one Gaussian band, sigma derived from the measured FWHM."""
    sigma = fwhm / 2.3548
    trace += height * np.exp(-0.5 * ((x - scan) / sigma) ** 2)


def _channel(scans, *, height=PEAK_HEIGHT, noise=NOISE_SD, seed=0,
             fwhm=FWHM_SCANS):
    """A channel holding Gaussian bands on a noisy baseline."""
    rng = np.random.default_rng(seed)
    x = np.arange(N_SCANS, dtype=float)
    y = np.full(N_SCANS, BASELINE) + rng.normal(0.0, noise, N_SCANS)
    for s in scans:
        _band(y, x, float(s), height, fwhm)
    return y


class _Doc:
    """The bits of ``TraceDocument`` that ``size_trace`` reads."""

    def __init__(self, acgt, well="A01", path="A01.rsd"):
        self.acgt = np.asarray(acgt, float)
        self.well = well
        self.path = path
        self.base_order = BASE_ORDER


def _doc(ladder_bp, sample_bp, **kw):
    acgt = np.zeros((N_SCANS, 4))
    acgt[:, acgt_index_for_channel(BASE_ORDER, LADDER_CHANNEL)] = _channel(
        [scan_of(b) for b in ladder_bp], height=PEAK_HEIGHT * 1.2, seed=11, **kw)
    acgt[:, acgt_index_for_channel(BASE_ORDER, SAMPLE_CHANNEL)] = _channel(
        [scan_of(b) for b in sample_bp], seed=12, **kw)
    return _Doc(acgt)


# --------------------------------------------------------------------------- #
# the bundled GeneScan 500 ROX table
# --------------------------------------------------------------------------- #
GENESCAN500 = (35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400,
               450, 490, 500)


# --------------------------------------------------------------------------- #
# tests
# --------------------------------------------------------------------------- #
def test_genescans500_is_recovered_from_a_realistic_trace():
    """Every bundled fragment is sized back to its own length.

    This is the end-to-end claim of the module: given a ladder on its own
    channel and a sample on another, the lengths come back.  Held to the shape
    of the machine's own traces rather than to a noise-free idealisation.
    """
    sample_bp = (62, 118, 175, 211, 275, 306, 420, 470)
    result = fs.size_trace(_doc(GENESCAN500, sample_bp),
                           ladder="genescan500_rox",
                           ladder_channel=LADDER_CHANNEL,
                           sample_channel=SAMPLE_CHANNEL,
                           base_order=BASE_ORDER)
    assert not result.warnings, result.warnings
    assert result.ladder.name.lower().startswith("genescan")
    assert result.anchors, "no ladder band was matched"

    got = [r["length_bp"] for r in result.rows]
    assert len(got) == len(sample_bp)
    for want, have in zip(sample_bp, got):
        assert abs(have - want) < 5.0, f"{want} bp sized as {have:.1f}"


def test_sizing_error_stays_under_a_few_bp():
    """The reported leave-one-out RMS is small, and it is not the fit residual.

    The PCHIP passes through its anchors, so its residual is identically zero
    and says nothing.  Leave-one-out is the honest number and is what a user
    should read off a real run.
    """
    sample_bp = (80, 140, 190, 230, 290, 330, 380, 430, 480)
    result = fs.size_trace(_doc(GENESCAN500, sample_bp),
                           ladder="genescan500_rox",
                           ladder_channel=LADDER_CHANNEL,
                           sample_channel=SAMPLE_CHANNEL,
                           base_order=BASE_ORDER)
    assert result.quality["rms_error_bp"] is not None
    assert result.quality["rms_error_bp"] < 5.0, result.quality


def test_a_dropped_ladder_band_does_not_shift_the_others():
    """Lose one fragment from the ladder and the rest keep their lengths.

    This is the failure a length table is prone to on a real run -- one band
    under the threshold, or one below the injection window.  Pairing peaks to
    lengths position by position instead of skipping the gap would push every
    later length off by a whole fragment.
    """
    sample_bp = (62, 118, 175, 275, 380, 470)
    dropped = tuple(bp for bp in GENESCAN500 if bp != 250)
    result = fs.size_trace(_doc(dropped, sample_bp),
                           ladder="genescan500_rox",
                           ladder_channel=LADDER_CHANNEL,
                           sample_channel=SAMPLE_CHANNEL,
                           base_order=BASE_ORDER)
    got = [r["length_bp"] for r in result.rows]
    assert len(got) == len(sample_bp)
    for want, have in zip(sample_bp, got):
        assert abs(have - want) < 6.0, f"{want} bp sized as {have:.1f}"


def test_peaks_outside_the_ladder_are_flagged_not_guessed():
    """A band beyond the ladder's span is marked out of range.

    Extrapolating a PCHIP past its anchors is cheap and looks like data; the
    flag is what stops a caller reading it as a length.
    """
    sample_bp = (60, 250, 700)          # 700 bp is past the 500 bp top band
    result = fs.size_trace(_doc(GENESCAN500, sample_bp),
                           ladder="genescan500_rox",
                           ladder_channel=LADDER_CHANNEL,
                           sample_channel=SAMPLE_CHANNEL,
                           base_order=BASE_ORDER)
    out = [r for r in result.rows if not r["in_range"]]
    assert out, "a 700 bp band against a 35-500 bp ladder was not flagged"
    assert all(r["length_bp"] > 500.0 for r in out)


@pytest.mark.parametrize("lengths", [(400, 425, 450, 475, 500,
                                      525, 550, 575, 600, 625, 650, 675,
                                      700, 725, 750, 775, 800, 825, 850,
                                      875, 900, 925, 950, 975, 1000)])
def test_geneflo1000_is_recovered_from_a_realistic_trace(lengths):
    """The 400-1000 bp table sizes correctly on the same trace shape.

    Runs the other end of the migration, where the log curve is flattest and
    bands sit closest together -- the hardest part of the range to resolve.
    """
    sample_bp = (455, 530, 610, 705, 790, 880, 960)
    result = fs.size_trace(_doc(lengths, sample_bp),
                           ladder="geneflo1000_rox",
                           ladder_channel=LADDER_CHANNEL,
                           sample_channel=SAMPLE_CHANNEL,
                           base_order=BASE_ORDER)
    assert not result.warnings, result.warnings
    got = [r["length_bp"] for r in result.rows]
    assert len(got) == len(sample_bp)
    for want, have in zip(sample_bp, got):
        assert abs(have - want) < 6.0, f"{want} bp sized as {have:.1f}"


def test_a_run_with_no_ladder_is_refused_loudly():
    """An empty ladder channel cannot produce a length.

    Sizing must fail visibly rather than return a curve fitted to noise: a
    confidently wrong bp value is worse than no value at all.
    """
    doc = _doc(GENESCAN500, (62, 118, 175))
    acgt = np.asarray(doc.acgt, float).copy()
    # blank the ladder channel to bare baseline plus noise
    rng = np.random.default_rng(3)
    acgt[:, acgt_index_for_channel(BASE_ORDER, LADDER_CHANNEL)] = (
        BASELINE + rng.normal(0.0, NOISE_SD, N_SCANS))
    result = fs.size_trace(_Doc(acgt), ladder="genescan500_rox",
                           ladder_channel=LADDER_CHANNEL,
                           sample_channel=SAMPLE_CHANNEL,
                           base_order=BASE_ORDER)
    assert result.warnings, "a ladderless run produced no warning"
    assert not result.rows or all(not r["in_range"] for r in result.rows)


def test_swapped_channels_are_called_out():
    """Pointing the sizer at the wrong channel warns.

    Reading a ladder off the sample channel (or sizing a ladder) mis-sizes
    every peak while looking entirely normal, so this has to be loud.
    """
    sample_bp = (62, 118, 175, 275, 380, 470)
    result = fs.size_trace(_doc(GENESCAN500, sample_bp),
                           ladder="genescan500_rox",
                           ladder_channel=SAMPLE_CHANNEL,
                           sample_channel=SAMPLE_CHANNEL,
                           base_order=BASE_ORDER)
    assert any("same" in w for w in result.warnings), result.warnings