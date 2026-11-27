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
from pathlib import Path

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


# --------------------------------------------------------------------------- #
# the whole plate: automatic ladder/channel, plate-guided rescue
# --------------------------------------------------------------------------- #
ET550 = fs.load_ladder("et550r").lengths
# MegaBACE-like spread: ~9 scans/bp at the short end, ~6 at 550 bp, so the
# 300/310 pair is two peaks as on the real plates
MIG = lambda bp: 900.0 + 9.0 * bp - 0.003 * bp ** 2  # noqa: E731


def plate_doc(*, ladder_ch=1, bleed=0.12, seed=0, stretch=1.0, junk=(), shift=None):
    """ET550-R in *ladder_ch* (bleeding into the next channel), alleles in
    channel 3, optional *junk* scans (taller than the bands) in the ladder
    channel; *stretch* slows this capillary; *shift* ({bp: offset}) makes a
    band migrate as if it were that much longer (a ladder-lot quirk)."""
    acgt = np.zeros((N_SCANS, 4))
    put = lambda ch, y: acgt.__setitem__((slice(None), acgt_index_for_channel("ACTG", ch)), y)  # noqa: E731
    lad = [500 + stretch * (MIG(b + (shift or {}).get(b, 0.0)) - 500) for b in ET550]
    lc = synth_channel(lad, height=1200.0, noise=2.0, seed=seed)
    if junk:
        lc = lc + synth_channel(junk, height=2500.0, seed=seed + 7) - BASELINE
    put(ladder_ch, lc)
    nxt = ladder_ch % 4 + 1
    put(nxt, BASELINE + bleed * (lc - BASELINE) + np.random.default_rng(seed + 3).normal(0, 2, N_SCANS))
    others = [c for c in (1, 2, 3, 4) if c not in (ladder_ch, nxt)]
    put(others[0], synth_channel([MIG(b) for b in (112, 116, 180, 184, 231, 239)],
                                 height=900.0, noise=2.0, seed=seed + 1))
    put(others[1], synth_channel([], noise=2.0, seed=seed + 2))
    return FakeDoc(acgt, well=f"W{seed:02d}", path=f"W{seed:02d}.rsd")


def test_detect_ladder_finds_kit_and_channel_not_the_bleed():
    docs = [plate_doc(ladder_ch=1, seed=s) for s in range(4)]
    det = fs.detect_ladder(docs)
    assert det["ladder"] == "et550r"          # not ILS600/ILS500, which share many lengths
    assert det["channel"] == 1                # not channel 2, which holds 12 % bleed
    assert det["wells_ok"] == 4
    assert det["rivals"] == []                # a clean plate leaves no doubt about the kit


def test_detect_ladder_finds_a_ladder_in_channel_four():
    det = fs.detect_ladder([plate_doc(ladder_ch=4, seed=s) for s in range(3)])
    assert (det["ladder"], det["channel"]) == ("et550r", 4)


def test_detect_ladder_reports_none_without_a_size_standard():
    rng = np.random.default_rng(5)
    docs = []
    for s in range(3):
        acgt = np.zeros((N_SCANS, 4))
        for k in range(4):     # busy sample channels, no ladder anywhere
            acgt[:, k] = synth_channel(sorted(rng.uniform(1500, 5800, 60)),
                                       height=float(rng.uniform(200, 900)), noise=2.0, seed=s * 4 + k)
        docs.append(FakeDoc(acgt))
    assert fs.detect_ladder(docs)["ladder"] is None


def test_sample_peaks_in_the_ladder_channel_are_flagged():
    # a sample's peaks bleeding into the ladder channel: scattered, mixed heights
    rng = np.random.default_rng(11)
    doc = plate_doc(junk=sorted(rng.uniform(MIG(60), MIG(550), 40)))
    res = fs.size_trace(doc, "et550r", ladder_channel=1, sample_channels=[3])
    assert res.quality["unexplained_tall"] > fs.MAX_UNEXPLAINED
    assert any("does not explain" in w for w in res.warnings)


def test_template_rescues_a_ladder_under_taller_junk():
    good = [fs.size_trace(plate_doc(seed=s, stretch=1.0 + 0.02 * s), "et550r",
                          ladder_channel=1, sample_channels=[3]) for s in range(4)]
    tmpl = fs.plate_template(good)
    assert tmpl and len(tmpl) == len(ET550)
    lad = [500 + 1.06 * (MIG(b) - 500) for b in ET550]
    junk = [s + 0.5 * (t - s) for s, t in zip(lad, lad[1:])][::2]   # between bands, taller
    bad = plate_doc(seed=9, stretch=1.06, junk=junk)
    pk = fs.detect_peaks(bad.acgt[:, acgt_index_for_channel("ACTG", 1)],
                         min_prominence_frac=fs.DEFAULT_LADDER_PROMINENCE_FRAC,
                         min_height_frac=fs.DEFAULT_LADDER_HEIGHT_FRAC)
    got = fs.align_to_template([p.scan for p in pk], ET550, tmpl)
    assert len(got) >= 0.8 * len(ET550)
    for scan, bp in got:
        assert abs(scan - (500 + 1.06 * (MIG(bp) - 500))) < 4


def test_template_refuses_a_channel_of_equal_small_peaks():
    good = [fs.size_trace(plate_doc(seed=s), "et550r", ladder_channel=1,
                          sample_channels=[3]) for s in range(4)]
    tmpl = fs.plate_template(good)
    scans = list(range(1500, 5900, 9))         # a peak every 9 scans, all alike
    got = fs.align_to_template(scans, ET550, tmpl, {s: 6.0 for s in scans})
    assert got == []


# --- per-run band correction ------------------------------------------------

def test_local_band_offsets_are_zero_on_a_smooth_ladder():
    anchors = [(MIG(b), b) for b in ET550]
    offs = fs.local_band_offsets(anchors)
    assert set(offs) == set(ET550[2:-2])          # interpolation only, ends skipped
    assert max(abs(v) for v in offs.values()) < 0.05


def _sized_run(n, shift=None, band_correct=True):
    docs = {f"W{s:02d}.rsd": plate_doc(seed=s, stretch=1.0 + 0.01 * (s % 5), shift=shift)
            for s in range(n)}
    return fs.size_files(list(docs), ladder="et550r", ladder_channel=1, matrix="none",
                         load=lambda p: docs[str(p)], band_correct=band_correct)


def test_run_corrects_a_band_that_runs_long_in_every_well():
    res = _sized_run(10, shift={125: 2.5})
    corr = res[0].band_corrections
    assert list(corr) == [125.0] and abs(corr[125.0] - 127.5) < 0.6
    assert all(r.band_corrections == corr for r in res)
    assert "125 -> " in res[0].detection
    std = {r["band_bp"]: r for r in res[0].standards}
    assert std[125.0]["band_used_bp"] == corr[125.0]
    assert std[150.0]["band_used_bp"] == 150.0
    # the allele at 116 bp sits under the bent part of the curve: corrected,
    # it is sized closer to its true length than without the correction
    def allele(rr):
        return min((r["length_bp"] for r in rr[0].rows if r["channel"] == 3),
                   key=lambda v: abs(v - 116))
    off = _sized_run(10, shift={125: 2.5}, band_correct=False)
    assert abs(allele(res) - 116) < abs(allele(off) - 116)
    assert res[0].quality["rms_error_bp"] < off[0].quality["rms_error_bp"]


def test_run_leaves_a_ladder_that_matches_its_table_alone():
    res = _sized_run(10)
    assert res[0].band_corrections == {}
    assert "corrected" not in res[0].detection


def test_too_few_wells_give_no_correction():
    res = _sized_run(fs.BAND_CORRECT_MIN_WELLS - 1, shift={125: 2.5})
    assert res[0].band_corrections == {}


# --- the run's ladder: nested kits, profile, cache ----------------------------

LONG = lambda bp: 800.0 + 6.0 * bp - 0.002 * bp ** 2  # noqa: E731  (runs to ~950 bp)


def kit_doc(kit, *, seed=0, ladder_ch=1, drop=()):
    """*kit*'s bands (minus *drop*) in *ladder_ch*, a few alleles in channel 3."""
    lad = fs.load_ladder(kit)
    acgt = np.zeros((N_SCANS, 4))
    put = lambda ch, y: acgt.__setitem__((slice(None), acgt_index_for_channel("ACTG", ch)), y)  # noqa: E731
    put(ladder_ch, synth_channel([LONG(b) for b in lad.lengths if b not in drop],
                                 height=1200.0, noise=2.0, seed=seed))
    for c in (1, 2, 3, 4):
        if c != ladder_ch:
            put(c, synth_channel([LONG(b) for b in (130, 250)] if c == 3 else [],
                                 height=900.0, noise=2.0, seed=seed + c))
    return FakeDoc(acgt, well=f"W{seed:02d}", path=f"W{seed:02d}.rsd")


@pytest.mark.parametrize("true,nested", [("et900r", "et550r"),
                                         ("promega_ils600", "promega_ils500"),
                                         ("mapmarker1500", "mapmarker1000")])
def test_detect_prefers_the_full_kit_over_one_nested_in_it(true, nested):
    det = fs.detect_ladder([kit_doc(true, seed=s) for s in range(4)],
                           kits=(true, nested))
    assert det["ladder"] == true


def test_detect_takes_the_kit_without_missing_bands_when_two_explain_the_peaks():
    # MapMarker 400 is Beckman SS400 without its 60 and 420 bp bands
    det = fs.detect_ladder([kit_doc("mapmarker400", seed=s) for s in range(4)],
                           kits=("beckman_ss400", "mapmarker400"))
    assert det["ladder"] == "mapmarker400"
    assert det["rivals"] == []


def test_run_ladder_round_trips_through_json():
    prof = fs.RunLadder(ladder=fs.load_ladder("et550r"), channel=3, status="found",
                        note="n", template={60.0: 1700.0, 75.0: 1800.5},
                        corrections={125.0: 127.2}, wells=96, wells_clean=40)
    back = fs.RunLadder.from_json(json.loads(json.dumps(prof.to_json())))
    assert back.ladder == prof.ladder and back.channel == 3
    assert back.template == prof.template and back.corrections == prof.corrections
    assert back.from_cache and back.reliable


def _run_on_disk(tmp_path, n=10):
    docs = {}
    for s in range(n):
        p = tmp_path / f"W{s:02d}.rsd"
        p.write_bytes(b"x")
        docs[p] = plate_doc(ladder_ch=3, seed=s, stretch=1.0 + 0.01 * (s % 5))
    calls = []

    def load(p):
        calls.append(Path(p))
        return docs[Path(p)]
    return sorted(docs), load, calls


def test_a_selection_is_sized_with_the_whole_runs_ladder(tmp_path):
    paths, load, _ = _run_on_disk(tmp_path)
    profs = {}
    one = fs.size_files(paths[:1], matrix="none", load=load, profiles=profs)
    prof = profs[str(tmp_path)]
    assert (prof.ladder.name, prof.channel, prof.wells) == ("MegaBACE ET550-R", 3, 10)
    assert prof.wells_clean == 10 and prof.template
    full = fs.size_files(paths, matrix="none", load=load)
    assert one[0].rows == full[0].rows


def test_run_ladder_is_cached_and_redone_when_a_file_changes(tmp_path):
    paths, load, calls = _run_on_disk(tmp_path)
    first = fs.find_run_ladder(paths, matrix=None, load=load)
    assert not first.from_cache and calls
    calls.clear()
    again = fs.find_run_ladder(paths, matrix=None, load=load)
    assert again.from_cache and not calls          # nothing re-read
    assert (again.ladder, again.channel) == (first.ladder, first.channel)
    paths[0].write_bytes(b"changed")
    assert not fs.find_run_ladder(paths, matrix=None, load=load).from_cache
