"""Fragment-length sizing (``fragment_sizing``).

A ladder run is the one place where getting the channel or the length table
wrong is silently expensive: the sample is still "sized", just to the wrong
curve.  These tests therefore pin the ladder registry, the peak<->length
alignment (including a missing ladder band) and the scanned-back accuracy on a
synthetic trace whose true lengths are known exactly.

Everything here is headless -- ``fragment_sizing`` must stay importable by the
scorer without Tk or matplotlib.
"""
import json

import numpy as np
import pytest

from analyzer_core import acgt_index_for_channel
import fragment_sizing as fs

N_SCANS = 6000
BASELINE = 5.0
# A monotone, CE-like migration: scan grows with log(bp).
SCAN_OF = lambda bp: 1200.0 + 260.0 * np.log(bp)  # noqa: E731


class FakeDoc:
    """The bits of ``TraceDocument`` that ``size_trace`` reads."""

    def __init__(self, acgt, well="A01", path="A01.rsd", base_order="ACTG"):
        self.acgt = np.asarray(acgt, float)
        self.well = well
        self.path = path
        self.base_order = base_order


def synth_channel(scans, height=1000.0, width=6.0, noise=0.0, seed=0):
    rng = np.random.default_rng(seed)
    x = np.arange(N_SCANS)
    y = np.full(N_SCANS, BASELINE) + rng.normal(0.0, noise, N_SCANS)
    for s in scans:
        y += height * np.exp(-0.5 * ((x - s) / width) ** 2)
    return y


def make_doc(sample_bp, ladder_bp, *, seed=0, base_order="ACTG"):
    acgt = np.zeros((N_SCANS, 4))
    acgt[:, acgt_index_for_channel(base_order, 2)] = synth_channel(
        [SCAN_OF(b) for b in sample_bp], height=800.0, noise=2.0, seed=seed)
    acgt[:, acgt_index_for_channel(base_order, 4)] = synth_channel(
        [SCAN_OF(b) for b in ladder_bp], height=1000.0, noise=2.0, seed=seed + 1)
    return FakeDoc(acgt, base_order=base_order)


# --------------------------------------------------------------------------- #
# the ladder registry
# --------------------------------------------------------------------------- #
def test_builtin_ladders_load_and_are_sorted():
    for key in fs.list_ladders():
        ladder = fs.load_ladder(key)
        assert list(ladder.lengths) == sorted(ladder.lengths)
        assert ladder.lengths[0] > 0
    geneflo = fs.load_ladder("geneflo1000_rox")
    assert geneflo.lengths == tuple(range(400, 1001, 25))
    assert geneflo.dye == "ROX"


def test_ladder_name_lookup_is_case_insensitive():
    assert fs.load_ladder("GeneScan500_ROX").name == fs.load_ladder(
        "genescan500_rox").name


def test_load_ladder_from_json(tmp_path):
    path = tmp_path / "my_ladder.json"
    path.write_text(json.dumps({
        "name": "My FAM ladder", "dye": "FAM",
        "lengths": [50, 100, 150, 200], "notes": "home-made",
    }), encoding="utf-8")
    ladder = fs.load_ladder(str(path))
    assert ladder.name == "My FAM ladder"
    assert ladder.dye == "FAM"
    assert ladder.lengths == (50.0, 100.0, 150.0, 200.0)


def test_load_ladder_treats_a_bare_list_as_lengths():
    ladder = fs.load_ladder("50, 100, 150,200")
    assert ladder.lengths == (50.0, 100.0, 150.0, 200.0)


def test_load_ladder_rejects_unknown_name():
    with pytest.raises(ValueError, match="Unknown ladder"):
        fs.load_ladder("not_a_real_ladder")


def test_ladder_sorts_unsorted_lengths_and_rejects_empty():
    ladder = fs.Ladder(name="x", lengths=(300, 100, 200))
    assert ladder.lengths == (100.0, 200.0, 300.0)
    with pytest.raises(ValueError):
        fs.Ladder(name="empty", lengths=())


# --------------------------------------------------------------------------- #
# peak detection and alignment
# --------------------------------------------------------------------------- #
def test_detect_peaks_finds_each_gaussian():
    scans = [1500, 1700, 2000, 2400]
    peaks = fs.detect_peaks(synth_channel(scans), min_prominence_frac=0.02,
                            min_height_frac=0.03)
    assert [p.scan for p in peaks] == scans


def test_align_peaks_pairs_directly_when_counts_match():
    lengths = [100, 200, 300, 400, 500]
    scans = [int(SCAN_OF(b)) for b in lengths]
    anchors = fs.align_peaks(scans, lengths)
    assert [round(b) for _, b in anchors] == lengths


def test_align_peaks_skips_a_missing_ladder_band():
    lengths = [100, 200, 300, 400, 500]
    dropped = 300
    scans = [int(SCAN_OF(b)) for b in lengths if b != dropped]
    anchors = fs.align_peaks(scans, lengths)
    assert [round(b) for _, b in anchors] == [b for b in lengths if b != dropped]


# --------------------------------------------------------------------------- #
# the sizing curve
# --------------------------------------------------------------------------- #
def test_size_curve_passes_through_its_anchors():
    scans = [1500, 1800, 2100, 2500]
    lengths = [100.0, 200.0, 350.0, 500.0]
    curve = fs.SizeCurve(scans, lengths)
    assert curve.bp(scans) == pytest.approx(lengths)


def test_size_curve_is_monotone_between_anchors():
    curve = fs.SizeCurve([1500, 1800, 2100, 2500], [100.0, 200.0, 350.0, 500.0])
    scan = np.linspace(1500, 2500, 200)
    bp = curve.bp(scan)
    assert np.all(np.diff(bp) > 0)


def test_size_curve_flags_extrapolated_points():
    curve = fs.SizeCurve([1500, 1800, 2100, 2500], [100.0, 200.0, 350.0, 500.0])
    in_range = curve.in_range([1400, 2000, 2600])
    assert list(in_range) == [False, True, False]
    # extrapolated values stay ordered
    assert curve.bp([1400])[0] < curve.bp([1500])[0]
    assert curve.bp([2600])[0] > curve.bp([2500])[0]


def test_ladder_fit_quality_is_undefined_below_four_anchors():
    assert fs.ladder_fit_quality([1, 2, 3], [100.0, 200.0, 300.0])["rms_error_bp"] is None
    q = fs.ladder_fit_quality([1, 2, 3, 4], [100.0, 200.0, 300.0, 400.0])
    assert q["rms_error_bp"] is not None


# --------------------------------------------------------------------------- #
# end to end on a synthetic well
# --------------------------------------------------------------------------- #
def test_size_trace_recovers_known_fragment_lengths():
    ladder = fs.load_ladder("genescan500_rox")
    truth = [90.0, 180.0, 275.0, 420.0]
    doc = make_doc(truth, ladder.lengths)
    res = fs.size_trace(doc, ladder, ladder_channel=4, sample_channel=2)
    assert len(res.rows) == len(truth)
    got = [r["length_bp"] for r in res.rows]
    assert got == pytest.approx(truth, abs=3.0)
    assert all(r["in_range"] for r in res.rows)
    assert res.quality["rms_error_bp"] is not None
    assert not res.warnings


def test_size_trace_flags_peaks_beyond_the_ladder():
    ladder = fs.load_ladder("genescan500_rox")
    doc = make_doc([120.0, 700.0], ladder.lengths)
    res = fs.size_trace(doc, ladder, ladder_channel=4, sample_channel=2)
    by_scan = {r["length_bp"]: r["in_range"] for r in res.rows}
    longest = max(by_scan)
    assert by_scan[longest] is False


def test_size_trace_warns_when_channels_are_swapped():
    ladder = fs.load_ladder("genescan500_rox")
    doc = make_doc([120.0], ladder.lengths)
    res = fs.size_trace(doc, ladder, ladder_channel=2, sample_channel=2)
    assert any("same" in w for w in res.warnings)


def test_size_trace_records_channels_and_base_letter():
    ladder = fs.load_ladder("genescan500_rox")
    doc = make_doc([120.0], ladder.lengths)
    res = fs.size_trace(doc, ladder, ladder_channel=4, sample_channel=2)
    assert res.rows[0]["channel"] == 2
    assert res.rows[0]["base"] == "C"      # Ch2 under ACTG
    summary = res.summary()
    assert summary["ladder_channel"] == 4 and summary["sample_channel"] == 2
