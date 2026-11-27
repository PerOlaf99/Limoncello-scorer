"""SNP genotyping engine (snp_genotyping): synthetic plates.

The real-data check is tools/genotyping_qc/benchmark_snp.py against Per's hand
calls (12 plates, 95 %).  These synthetic plates pin the behaviour: the
standard is found with no marks, in whichever channel carries it, past a
channel the sample bleeds into; wells drift and their spacing varies; calls
cover hom1 / hom2 / het / AI het / merged heteroduplexes / no product.
"""
import numpy as np
import pytest

import snp_genotyping as sg

N = 4000
ORDER = "ACTG"                       # Ch1=A, Ch2=C, Ch3=T, Ch4=G


def _peak(t, c, h, w=2.5):
    return h * np.exp(-0.5 * ((t - c) / w) ** 2)


def _well(rng, genotype, gaps=(90, 210, 90), ai=1.0, product=8000.0,
          std_h=1500.0, sample_ch=2, std_ch=3, bleed_ch=1, offset=-12):
    """One well's acgt array.  *genotype* hom1/hom2/het/none."""
    t = np.arange(N, dtype=float)
    shift = rng.uniform(2000, 2350)
    jit = rng.normal(1.0, 0.06, 3)
    g = np.asarray(gaps, float) * jit
    pos = shift + np.concatenate([[0], np.cumsum(g)])
    std = np.zeros(N)
    for i, p in enumerate(pos):
        if gaps[2] == 0 and i == 3:
            continue
        std += _peak(t, p, std_h * rng.uniform(0.8, 1.2))
        std += _peak(t, p + 20, 0.35 * std_h)          # the satellite every band has
    samp = np.zeros(N)
    sp = pos + offset
    if genotype == "hom1":
        samp += _peak(t, sp[0], product)
    elif genotype == "hom2":
        samp += _peak(t, sp[1], product)
    elif genotype == "het":
        a1, a2 = product, product * ai
        samp += _peak(t, sp[0], a1) + _peak(t, sp[1], a2)
        het = 0.6 * np.sqrt(a1 * a2)
        samp += _peak(t, sp[2], het)
        if gaps[2] > 0:
            samp += _peak(t, sp[3], het)
    # primer/dye front: tall, uneven peaks early in every channel
    front = sum(_peak(t, 1650 + 17 * k, 1200 * (1 + k % 3)) for k in range(5))
    a = np.zeros((N, 4))
    col = {ch: sg.acgt_index_for_channel(ORDER, ch) for ch in (1, 2, 3, 4)}
    a[:, col[sample_ch]] = samp + 0.3 * front
    a[:, col[std_ch]] = std + front + 0.01 * samp
    a[:, col[bleed_ch]] = 0.25 * samp + 0.2 * front
    a += 100 + rng.normal(0, 5, a.shape)
    return a


def _plate(seed=1, genotypes=None, **kw):
    rng = np.random.default_rng(seed)
    genotypes = genotypes or (["hom1"] * 30 + ["hom2"] * 30 + ["het"] * 30 + ["none"] * 6)
    wells = [f"{r}{c:02d}" for c in range(1, 13) for r in "ABCDEFGH"]
    truth, arrs = {}, {}
    for w, gt in zip(wells, genotypes):
        truth[w] = gt
        arrs[w] = _well(rng, gt, **kw)
    return truth, arrs


def _agree(truth, rows):
    norm = {"no-call": "none"}
    return {w: (truth[w], norm.get(rows[w]["call"], rows[w]["call"]))
            for w in truth if norm.get(rows[w]["call"], rows[w]["call"]) != truth[w]}


def test_finds_the_standard_with_no_marks():
    truth, arrs = _plate()
    ps = sg.find_plate_standard(arrs, ORDER)
    assert ps is not None
    assert ps.is_col == sg.acgt_index_for_channel(ORDER, 3)
    assert ps.sample_col == sg.acgt_index_for_channel(ORDER, 2)
    assert abs(ps.gaps[0] - 90) < 12 and abs(ps.gaps[1] - 210) < 20
    assert abs(ps.gaps[2] - 90) < 15
    assert ps.wells_agree >= 60


def test_calls_every_genotype():
    truth, arrs = _plate()
    rows = sg.genotype_plate(arrs, ORDER)
    assert _agree(truth, rows) == {}


def test_standard_in_another_channel_and_a_bleed_channel():
    truth, arrs = _plate(seed=2, std_ch=4, bleed_ch=3)
    rows = sg.genotype_plate(arrs, ORDER)
    assert "standard on Ch4" in next(iter(rows.values()))["plate_note"]
    assert _agree(truth, rows) == {}


def test_merged_heteroduplexes():
    truth, arrs = _plate(seed=3, gaps=(60, 230, 0))
    rows = sg.genotype_plate(arrs, ORDER)
    assert "merged" in next(iter(rows.values()))["plate_note"]
    assert _agree(truth, rows) == {}
    assert all("merged-het" in rows[w]["flags"] for w in truth if truth[w] != "none")


def test_allelic_imbalance_is_het_and_flagged():
    gts = ["hom1"] * 30 + ["hom2"] * 30 + ["het"] * 36
    rng = np.random.default_rng(4)
    wells = [f"{r}{c:02d}" for c in range(1, 13) for r in "ABCDEFGH"]
    arrs, truth = {}, {}
    for i, (w, gt) in enumerate(zip(wells, gts)):
        ai = 0.15 if (gt == "het" and i % 2) else 1.0
        arrs[w] = _well(rng, gt, ai=ai)
        truth[w] = (gt, ai)
    rows = sg.genotype_plate(arrs, ORDER)
    for w, (gt, ai) in truth.items():
        assert rows[w]["call"] == gt, (w, gt, rows[w])
        if gt == "het":
            assert ("ai" in rows[w]["flags"].split(",")) == (ai < 1.0), (w, rows[w])
            assert 0.0 < rows[w]["frac2"] < 1.0


def test_directional_fraction():
    truth, arrs = _plate(seed=5)
    rows = sg.genotype_plate(arrs, ORDER)
    for w, gt in truth.items():
        if gt == "hom1":
            assert rows[w]["frac2"] < 0.1
        elif gt == "hom2":
            assert rows[w]["frac2"] > 0.9
        elif gt == "het":
            assert 0.3 < rows[w]["frac2"] < 0.7


def test_no_product_is_a_no_call_with_a_reason():
    truth, arrs = _plate(seed=6)
    rows = sg.genotype_plate(arrs, ORDER)
    for w, gt in truth.items():
        if gt == "none":
            assert rows[w]["call"] == "no-call"
            assert "no product" in rows[w]["reason"]


def test_plate_without_a_standard():
    rng = np.random.default_rng(7)
    arrs = {f"A{c:02d}": rng.normal(100, 5, (N, 4)) for c in range(1, 13)}
    rows = sg.genotype_plate(arrs, ORDER)
    assert all(r["call"] == "no-call" for r in rows.values())
    assert all("no internal standard" in r["reason"] for r in rows.values())


@pytest.mark.parametrize("h1,h2,he,want", [
    (5000, 30, 20, "hom1"),
    (30, 5000, 20, "hom2"),
    (3000, 2500, 1500, "het"),
    (3000, 2500, 0, "het"),          # heteroduplexes lost, two homoduplexes still say two alleles
    (3000, 300, 900, "het"),         # AI het
    (3000, 60, 2000, "hom1"),        # a stray band without a second homoduplex
    (20, 10, 5, "no-call"),
])
def test_call_rule(h1, h2, he, want):
    assert sg.call_well(h1, h2, [he])[0] == want


# --------------------------------------------------------------- alignment
def test_align_map_methods():
    std, ref = [2100, 2190, 2400, 2490], [2000, 2100, 2300, 2400]
    pw = sg.align_map(std, ref, "piecewise")
    assert np.allclose(pw(std), ref)                       # exact through the bands
    assert pw([2145])[0] == pytest.approx(2050)            # straight in between
    sc = sg.align_map(std, ref, "scale")                   # Per's sum normalisation
    assert sc([sum(std) / 4])[0] == pytest.approx(sum(ref) / 4)
    af = sg.align_map(std, ref, "affine")
    assert np.abs(af(std) - ref).max() < 15
    with pytest.raises(ValueError):
        sg.align_map(std, ref, "warp")


def test_plate_reference_uses_hetero2_only_where_present():
    ref = sg.plate_reference({"A": [10, 20, 40, 50], "B": [12, 22, 42], "C": [11, 21, 41, 51],
                              "D": [10, 20, 40, 50], "E": None})
    assert ref[:3] == [10.5, 20.5, 40.5] and ref[3] == 50


def test_alignment_lines_up_the_sample_across_wells():
    truth, arrs = _plate(seed=8)
    rows = sg.genotype_plate(arrs, ORDER)
    aligned, ref = sg.align_plate(arrs, rows)
    col = sg.acgt_index_for_channel(ORDER, 2)
    raw = [int(np.argmax(arrs[w][1500:, col])) + 1500 for w in aligned if truth[w] == "hom1"]
    lined = [int(np.argmax(aligned[w][1500:, col])) + 1500 for w in aligned if truth[w] == "hom1"]
    assert np.std(raw) > 50                                # capillaries drift ~100 scans
    assert np.std(lined) < 4                               # one plate axis
    assert abs(np.median(lined) - (ref[0] - 12)) < 4       # beside the reference Homo1


def test_failed_pcr_plate_is_all_no_calls():
    """Standard and dye front in every well, no product anywhere: the
    standard is found without the sample and nothing is called."""
    truth, arrs = _plate(seed=12, genotypes=["none"] * 96)
    # the unused primers carry the sample's dye: the tallest peaks of the plate
    col = sg.acgt_index_for_channel(ORDER, 2)
    t = np.arange(N, dtype=float)
    for a in arrs.values():
        a[:, col] += _peak(t, 1580, 30000.0) + _peak(t, 1600, 20000.0)
    rows = sg.genotype_plate(arrs, ORDER)
    assert "no sample product" in next(iter(rows.values()))["plate_note"]
    assert all(r["call"] == "no-call" for r in rows.values())
    assert sum(bool(r["std_scans"]) for r in rows.values()) >= 90


def test_a_placement_far_from_the_plate_is_placed_again():
    """A well whose bands fit the spacing on the dye front as well as on the
    standard is placed where the rest of the plate has its standard."""
    truth, arrs = _plate(seed=13)
    rows = sg.genotype_plate(arrs, ORDER)
    h1 = np.array([int(r["std_scans"].split("/")[0]) for r in rows.values() if r["std_scans"]])
    assert np.all(np.abs(h1 - np.median(h1)) <= sg.POSITION_WINDOW)
