"""Fragment sizing against the MegaBACE Genetic Profiler 2.2 demo plates.

The GP 2.2 CD ships two 96-well MegaBACE 1000 runs with an ET400-R ladder in
ET-ROX on physical channel 1 (the footer names the dyes ET-ROX, FAM, NED, HEX
on channels 1-4) and, for the Genotyping Test Plate, the fragment sizes GP
expects in each dye (``CS_Ladder_info.tdb``, same channel numbering).

The alignment tests below use ladder-channel peak scans taken from three of
those wells, so they always run.  The whole-plate tests need the demo data
itself (not shipped: GE software content); point ``LIMONCELLO_GP22_DEMO`` at the
folder holding ``GPScoreCard_DemoData``, ``GP20_Demo_Data`` and ``gp_tables``.
"""
import os
from pathlib import Path

import pytest

import fragment_sizing as fs

# ET400-R band table, from the GP 2.2 test-plate definition (test_plate.tdb)
ET400R = [60, 90, 100, 120, 150, 160, 170, 190, 200, 220, 250, 270, 290, 300,
          310, 330, 350, 360, 380, 400]


def test_et400r_is_bundled():
    assert list(fs.load_ladder("et400r").lengths) == ET400R


def _bands(anchors):
    return {int(scan): int(bp) for scan, bp in anchors}


# --------------------------------------------------------------------------- #
# alignment on real ladder-channel peaks (no data needed)
# --------------------------------------------------------------------------- #
def test_extra_peak_before_the_first_band_is_not_taken_as_60bp():
    # GPScoreCard_DemoData A01: one strong extra peak (2071) ahead of the
    # ladder, whose 60 bp band is at 2183
    scans = [2071, 2183, 2405, 2478, 2624, 2856, 2935, 3017, 3180, 3263, 3435,
             3700, 3880, 4064, 4158, 4252, 4439, 4632, 4728, 4918, 5110]
    got = _bands(fs.align_peaks(scans, ET400R))
    assert got == dict(zip(scans[1:], ET400R))


def test_primer_cluster_and_dye_bleed_between_bands():
    # GP20_Demo_Data A02: primer/dye peaks at 1806-2013 ahead of the ladder,
    # sample-dye bleed at 2671, 3473 and 4004 between bands
    scans = [1806, 1816, 1823, 1859, 1867, 1997, 2013, 2108, 2330, 2404, 2551,
             2671, 2787, 2868, 2951, 3117, 3202, 3379, 3473, 3649, 3834, 4004,
             4022, 4118, 4214, 4407, 4603, 4702, 4896, 5090]
    ladder = [2108, 2330, 2404, 2551, 2787, 2868, 2951, 3117, 3202, 3379, 3649,
              3834, 4022, 4118, 4214, 4407, 4603, 4702, 4896, 5090]
    assert _bands(fs.align_peaks(scans, ET400R)) == dict(zip(ladder, ET400R))


def test_bleed_peak_right_next_to_a_band():
    # GP20_Demo_Data A03: bleed at 2744 and 2961 sits close to the 150 and
    # 170 bp bands (2847, 3014)
    scans = [1838, 1849, 1856, 1867, 1884, 1891, 1901, 2031, 2050, 2151, 2379,
             2455, 2606, 2744, 2847, 2930, 2961, 3014, 3184, 3271, 3452, 3584,
             3728, 3917, 4109, 4207, 4305, 4501, 4701, 4802, 5000, 5199]
    ladder = [2151, 2379, 2455, 2606, 2847, 2930, 3014, 3184, 3271, 3452, 3728,
              3917, 4109, 4207, 4305, 4501, 4701, 4802, 5000, 5199]
    assert _bands(fs.align_peaks(scans, ET400R)) == dict(zip(ladder, ET400R))


def test_ladder_candidates_drop_weak_bleed():
    peaks = [fs.Peak(s, h, h) for s, h in
             [(1913, 150), (2071, 472), (2183, 712), (2252, 152), (2405, 720),
              (2541, 79), (2624, 810)]]
    kept = [p.scan for p in fs.ladder_candidates(peaks, 4)]
    assert kept == [2071, 2183, 2405, 2624]


# --------------------------------------------------------------------------- #
# whole plates (needs the GP 2.2 demo data)
# --------------------------------------------------------------------------- #
DEMO = Path(os.environ.get("LIMONCELLO_GP22_DEMO",
                           "/run/media/per/NVMe Disk/limoncello-data/"
                           "fragment_sizing/GP22_demo"))
needs_demo = pytest.mark.skipif(
    not (DEMO / "GPScoreCard_DemoData").is_dir(),
    reason="GP 2.2 demo data not present (set LIMONCELLO_GP22_DEMO)")

ORDER = "TGCA"       # any order will do, as long as loader and sizer share it


def _size(run, well, sample_channel=2):
    """Ladder on physical Ch1 (ET-ROX); sample default Ch2 (FAM)."""
    from analyzer_core import load_trace
    doc = load_trace(DEMO / run / f"{well}.rsd", base_order=ORDER)
    return fs.size_trace(doc, fs.ladder_from_lengths(ET400R), ladder_channel=1,
                         sample_channel=sample_channel, base_order=ORDER)


@needs_demo
@pytest.mark.parametrize("run,well", [
    ("GPScoreCard_DemoData", "B01"), ("GPScoreCard_DemoData", "F05"),
    ("GPScoreCard_DemoData", "H12"), ("GP20_Demo_Data", "A02"),
    ("GP20_Demo_Data", "A03"), ("GP20_Demo_Data", "B06")])
def test_demo_ladder_fits_to_a_fraction_of_a_base(run, well):
    res = _size(run, well)
    assert len(res.anchors) == len(ET400R)
    assert res.quality["rms_error_bp"] < 0.3
    assert not res.warnings


@needs_demo
@pytest.mark.parametrize("run,well", [
    ("GPScoreCard_DemoData", "A04"), ("GPScoreCard_DemoData", "A08"),
    ("GP20_Demo_Data", "C12"), ("GP20_Demo_Data", "G12")])
def test_demo_spectral_standard_wells_warn(run, well):
    # single-dye matrix-standard wells hold no ROX ladder
    assert _size(run, well).warnings


@needs_demo
@pytest.mark.parametrize("well", ["B01", "C07", "E10", "H06"])
def test_demo_test_plate_fragments_match_genetic_profiler(well):
    expected = {}
    for line in (DEMO / "gp_tables" / "CS_Ladder_info.tdb").read_text().split("\n"):
        parts = line.split()
        if len(parts) == 2 and int(parts[1]) != 1:      # 1 = the ladder itself
            expected.setdefault(int(parts[1]), []).append(float(parts[0]))
    for channel, sizes in expected.items():
        rows = _size("GPScoreCard_DemoData", well, channel).rows
        top = max(r["height"] for r in rows)
        for bp in sizes:
            near = [r for r in rows if abs(r["length_bp"] - bp) <= 2
                    and r["height"] >= 0.1 * top]
            assert near, f"{well} ch{channel}: no peak near {bp} bp"
            best = max(near, key=lambda r: r["height"])
            assert best["length_bp"] == pytest.approx(bp, abs=0.4)
