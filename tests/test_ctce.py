"""CTCE engine on synthetic runs: what moves with temperature is a melting
duplex, what holds still is static, the run order follows the recorded
temperature, satellites ride with their parent, and areas are measured."""
import sys
from pathlib import Path

import numpy as np
from scipy.stats import exponnorm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import ctce  # noqa: E402

N = 4000


def emg(x, mu, sigma, tau, area):
    """Exponentially modified Gaussian: sharp front, exponential tail."""
    return area * exponnorm.pdf(x, tau / sigma, loc=mu, scale=sigma)


def gauss(x, mu, sigma, area):
    return area / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def trace(shift, static_drift=0.0, satellite=True, scale=1.0):
    x = np.arange(N, dtype=float)
    y = 200.0 + gauss(x, 1900, 3.0, 2.0e5)                            # primer / anchor
    y += gauss(x, 2600 + static_drift, 6.0, 1.0e5)                    # ssDNA: symmetric, static
    y += emg(x, 2250 + shift, 2.5, 8.0, 2.4e5)                        # melting duplex
    if satellite:
        y += emg(x, 2268 + shift, 2.5, 8.0, 0.8e5)                    # its stutter / +A
    return y * scale


def run(temp, shift, **kw):
    y = trace(shift, **kw)
    y = ctce.channel_trace(np.column_stack([y * 0.05, y, y * 0.05, y * 0.05]), "ACGT", 2)
    rt = ctce.RunTrace(path=Path(f"x/PL_Run0{int(temp) - 50}/A01.rsd"), run=f"PL_Run0{int(temp) - 50}",
                       plate="PL", well="A01", temperature=temp, start="", channel=2, y=y)
    rt.peaks = ctce.detect_peaks(y)
    rt.anchor = ctce.anchor_scan(rt)
    return rt


def analyse(runs, monkeypatch):
    it = iter(runs)
    monkeypatch.setattr(ctce, "load_run_trace", lambda p, ch=None: next(it))
    return ctce.analyse_well(["a", "b", "c"][:len(runs)])


def test_the_duplex_moves_later_and_the_ssdna_holds_still(monkeypatch):
    res = analyse([run(53, 0, static_drift=0), run(54, 40, static_drift=-5),
                   run(55, 130, static_drift=-10)], monkeypatch)
    assert [s > 15 for s in res.shifts] == [True, True]
    for r in res.runs:
        roles = {p.role for p in r.peaks if p.height > 0.08 * r.y.max()}
        assert {"anchor", "static", "melting", "satellite"} <= roles
        duplex = max((p for p in r.peaks if p.role == "melting"), key=lambda p: p.height)
        assert duplex.asym > 1.3                      # sharp front, long tail
        ss = next(p for p in r.peaks if p.role == "static")
        assert ss.asym < 1.3                          # symmetric


def test_runs_are_ordered_by_temperature_not_by_run_number(monkeypatch):
    # the hottest run was done first: the direction must still come out right
    res = analyse([run(55, 130), run(53, 0), run(54, 40)], monkeypatch)
    assert [r.temperature for r in res.runs] == [53, 54, 55]
    assert all(s > 15 for s in res.shifts)


def test_satellite_area_is_added_to_its_parent(monkeypatch):
    res = analyse([run(53, 0), run(54, 40), run(55, 130)], monkeypatch)
    assert len(res.areas) == 1                        # one species, its satellite merged
    best = res.runs[res.best]
    sat = next(p for p in best.peaks if p.role == "satellite")
    par = best.peaks[sat.parent]
    assert abs(res.areas[0] - (sat.area + par.area)) < 1e-6
    assert any("satellite" in f for f in res.flags)


def test_valley_areas_are_close_to_the_true_areas():
    y = trace(0, satellite=False) - 200.0
    peaks = ctce.detect_peaks(y)
    duplex = min(peaks, key=lambda p: abs(p.pos - 2255))
    assert abs(duplex.area - 2.4e5) / 2.4e5 < 0.03
    ss = min(peaks, key=lambda p: abs(p.pos - 2600))
    assert abs(ss.area - 1.0e5) / 1.0e5 < 0.03


def test_one_temperature_is_flagged(monkeypatch):
    res = analyse([run(55, 0), run(55, 0)], monkeypatch)
    assert any("one temperature" in f for f in res.flags)


def test_four_peaks_give_a_mass_action_minor_fraction():
    # two homoduplexes then two heteroduplexes, p = 0.2: hom 0.04, 0.64; het 0.16 each
    r = ctce.WellResult(plate="P", well="A01", runs=[], shifts=[], regimes=[], best=None)
    areas = [0.64, 0.04, 0.16, 0.16]
    total = sum(areas)
    p = (areas[0] + (areas[2] + areas[3]) / 2) / total
    r.minor_fraction = min(p, 1 - p)
    assert abs(r.minor_fraction - 0.2) < 1e-9


def test_plate_of_strips_the_run_number():
    assert ctce.plate_of("040723_mt_nucl_#19_20_dil_curveRun02") == "040723_mt_nucl_#19_20_dil_curve"


# ------------------------------------------------------------- user guidance
def test_template_finds_the_clicked_pattern_shifted_in_another_well():
    tpl = ctce.Template(offsets={54.0: [350.0, 368.0]}, n_wells=1)   # duplex + its partner
    rt = run(54, 25)                      # this well's pattern sits 25 scans later
    assert ctce._apply_template(rt, tpl.at(54)) == 2
    labelled = sorted((p for p in rt.peaks if p.label), key=lambda p: p.pos)
    assert [p.label for p in labelled] == ["P1", "P2"]
    assert abs(labelled[0].raw_pos - (rt.anchor + 375)) <= 3


def test_clicked_peaks_are_never_merged_as_a_satellite(monkeypatch):
    tpl = ctce.Template(offsets={53.0: [350.0, 368.0], 54.0: [390.0, 408.0]}, n_wells=1)
    runs = [run(53, 0), run(54, 40)]
    it = iter(runs)
    monkeypatch.setattr(ctce, "load_run_trace", lambda p, ch=None: next(it))
    res = ctce.analyse_well(["a", "b"], template=tpl)
    assert res.method == "template"
    assert len(res.areas) == 2            # both clicked peaks measured on their own


def test_one_temperature_uses_shape_and_replicates(monkeypatch):
    runs = [run(55, 0), run(55, 0, scale=1.05)]
    it = iter(runs)
    monkeypatch.setattr(ctce, "load_run_trace", lambda p, ch=None: next(it))
    res = ctce.analyse_well(["a", "b"])
    assert res.method == "shape"
    best = res.runs[res.best]
    duplex = max((p for p in best.peaks if p.role == "melting"), key=lambda p: p.height)
    assert duplex.asym >= ctce.DUPLEX_ASYM
    assert all(p.role != "melting" for p in best.peaks if p.asym < ctce.SSDNA_ASYM)
    assert res.fraction_sd and all(s < 0.01 for s in res.fraction_sd)   # identical replicates


# ------------------------------------------------------------ species model
def _species_well(k, p=0.15, rng=None):
    """A well with Homo1 / Homo2 / Hetero1 / Hetero2 (each with its stutter), a
    static ssDNA peak and the early peak, on a well-specific migration scale."""
    rng = rng or np.random.default_rng(k)
    x = np.arange(N, dtype=float)
    anchor = 1850 + 12 * k
    total = 4.0e5 * (1 + 0.3 * (k % 3))
    shares = [(1 - p) ** 2, p ** 2, p * (1 - p), p * (1 - p)]
    h2 = 1.31 * anchor * (1 + 0.004 * ((k % 5) - 2))
    pos = [h2, h2 * 1.046, h2 * 1.095, h2 * 1.107]
    y = gauss(x, anchor, 3.0, 2.0e5) + gauss(x, 1.226 * anchor, 4.0, 6.0e4)
    for c, s in zip(pos, shares):
        y += emg(x, c, 2.0, 4.0, 0.7 * s * total) + emg(x, c + 16, 2.0, 4.0, 0.3 * s * total)
    y += rng.normal(0, 5, N)
    rt = ctce.RunTrace(path=Path(f"x/PL_Run01/W{k:02d}.rsd"), run="PL_Run01", plate="PL",
                       well=f"W{k:02d}", temperature=55.0, start="", channel=2, y=y)
    rt.peaks = ctce.detect_peaks(y)
    rt.anchor = ctce.anchor_scan(rt)
    truth = [s * total for s in shares]
    return rt, [int(round(c)) for c in pos], truth


def test_species_are_found_as_a_plate_pattern_and_measured_as_blocks():
    wells = [_species_well(k) for k in range(12)]
    results = [ctce.WellResult(plate="PL", well=rt.well, runs=[rt], shifts=[], regimes=[],
                               best=0, method="shape") for rt, _, _ in wells]
    pattern = ctce.quantify_species(results)
    assert pattern.r3 is not None and abs(pattern.r3 - 1.046) < 0.006
    assert any(abs(s - 1.226) < 0.012 for s in pattern.static)       # ssDNA learnt as static
    for res, (rt, pos, truth) in zip(results, wells):
        assert res.method == "species"
        assert res.species == ["Homo1", "Homo2", "Hetero1", "Hetero2"]
        marked = {p.label: p.raw_pos for p in rt.peaks if p.label}
        for name, c in zip(ctce.SPECIES_NAMES, pos):
            assert abs(marked[name] - c) <= 4                         # the main peak, not the stutter
        for a, t in zip(res.areas[:2], truth[:2]):
            assert abs(a - t) / t < 0.05                              # block = main + stutter
        # each heteroduplex is measured above its own line, split at the
        # valley between them, as an operator drags them: what lies below the
        # valley level is not counted, so the sum reads low when they overlap
        het, het_true = sum(res.areas[2:]), sum(truth[2:])
        assert 0.6 < het / het_true < 1.05
        assert abs(res.minor_fraction - 0.15) < 0.03


def test_a_clipped_peak_is_measured_in_another_channel():
    from types import SimpleNamespace
    x = np.arange(N, dtype=float)
    true = 300.0 + emg(x, 2400, 2.0, 4.0, 1.5e6) + emg(x, 2520, 2.0, 4.0, 6.0e4)
    raw = np.column_stack([0.12 * true + 40, 0.02 * true + 30, np.minimum(true, 65535.0), 0.05 * true + 20])
    rt = ctce.RunTrace(path=Path("x/PL_Run01/A01.rsd"), run="PL_Run01", plate="PL", well="A01",
                       temperature=55.0, start="", channel=3,
                       y=ctce.channel_trace(raw, "ACGT", 3))
    rt._doc, rt._order = SimpleNamespace(acgt=raw), "ACGT"
    assert true.max() > 65535                                         # the sample channel clips
    y = ctce.quant_trace(rt)
    assert rt._desat[0] == "A"                                        # the strongest clean channel
    big = ctce.species_areas(y, [2400])[0]
    assert abs(big - 1.5e6) / 1.5e6 < 0.03
    clipped = ctce.species_areas(ctce.channel_trace(raw, "ACGT", 3, ctce.SPECIES_BASELINE), [2400])[0]
    assert clipped < 0.8 * 1.5e6                                      # what clipping would have cost


def test_homo1_marked_static_by_the_other_temperature_is_still_homo1():
    # between close temperatures Homo1 barely moves and the runs pair it
    # with itself as "static"; a far bigger peak is still Homo1, not the
    # small peak before it
    rt, pos, _truth = _species_well(1)
    y = rt.y + gauss(np.arange(N, dtype=float), pos[0] - 18, 2.5, 6.0e3)   # a small peak just before
    m = ctce.find_species(y, rt.anchor, static_scans=[pos[0], pos[0] + 16])
    assert m is not None and abs(m[0] - pos[0]) <= 3


def test_a_bump_on_the_rise_of_homo2_is_not_homo2():
    rt, pos, _truth = _species_well(2)
    x = np.arange(N, dtype=float)
    bump = pos[1] - 11
    y = rt.y + gauss(x, bump, 2.0, 1.5e3)                                    # a small bump before Homo2
    # ... where the plate expects Homo2: its footprint borrows Homo2's area
    m = ctce.find_species(y, rt.anchor, ctce.SpeciesPattern(r3=bump / pos[0]))
    assert m is not None and abs(m[1] - pos[1]) <= 3


def test_merging_heteroduplexes_are_measured_as_one_area():
    # two heteroduplexes so close the valley between them stays well above
    # the floor: they count as one area under one line (each half goes to the
    # homoduplexes, so the split does not matter), not two cut at the valley
    x = np.arange(N, dtype=float)
    y = 300.0 + emg(x, 2600, 4.0, 5.0, 1.0e5) + emg(x, 2620, 4.0, 5.0, 1.0e5)
    c1, c2 = 2602, 2622
    v = c1 + int(np.argmin(y[c1:c2]))
    assert (y[v] - 300) > ctce.MERGED_VALLEY * (y[c2] - 300)          # they have begun to merge
    blocks = ctce.species_blocks(y, [c1, c2], stutter=16)
    (lo1, hi1, a1, b1), (lo2, hi2, a2, b2) = blocks
    slope = (b2 - a1) / (hi2 - lo1)
    assert lo2 == hi1 + 1 and abs(a2 - b1 - slope) < 1e-6             # one line, split at the valley
    total = sum(ctce.block_area(y, *bl) for bl in blocks)
    assert abs(total - 2.0e5) / 2.0e5 < 0.08


def test_heteroduplex_doublets_are_paired_in_order():
    # Hetero1 and Hetero2, each with its +A stutter 13 scans on, the second
    # 50 scans after the first: the stutter (as tall as its main) is not
    # Hetero2, and the doublets are not split by a gap rule
    wells = [_species_well(k) for k in range(12)]
    results = [ctce.WellResult(plate="PL", well=rt.well, runs=[rt], shifts=[], regimes=[],
                               best=0, method="shape") for rt, _, _ in wells]
    pattern = ctce.quantify_species(results)
    x = np.arange(N, dtype=float)
    rng = np.random.default_rng(99)
    anchor, h1 = 1900, 2400
    h2 = int(round(h1 * pattern.r3))
    y = gauss(x, anchor, 3.0, 2.0e5) + rng.normal(0, 5, N)
    for c, a in [(h1, 3.0e5), (h2, 3.0e4)]:
        y += emg(x, c, 2.0, 4.0, 0.7 * a) + emg(x, c + 16, 2.0, 4.0, 0.3 * a)
    het1 = h2 + 95
    het2 = het1 + 50
    for c in (het1, het2):
        y += emg(x, c, 2.0, 3.0, 1.2e4) + emg(x, c + 13, 2.0, 3.0, 1.2e4)
    mains = ctce.find_species(y, anchor, pattern)
    assert mains is not None and len(mains) == 4
    assert abs(mains[2] - het1) <= 3 and abs(mains[3] - het2) <= 3


def test_width_stops_at_a_valley():
    # two peaks that never fall to half height between them: each is narrow
    x = np.arange(400, dtype=float)
    y = gauss(x, 200, 3.0, 1.0e4) + gauss(x, 209, 3.0, 1.0e4)
    assert y[204] > 0.5 * y[200]                                       # merged above half height
    assert ctce._fwhm(y, 200) <= 12
    broad = gauss(x, 200, 25.0, 1.0e4)
    assert ctce._fwhm(broad, 200) > 40


def test_a_block_starts_at_the_foot_of_a_sharp_peak():
    # a sharp peak whose 5-scan mean tops out a scan before the raw apex
    # (clean Pfu products): the block must still start at the foot of the
    # rise, not at the peak itself -- that cut heteroduplexes to a tenth
    ys = np.zeros(200)
    c = 100
    ys[c - 6:c + 1] = [50, 200, 500, 800, 900, 950, 1000]
    ys[c + 1:c + 4] = [100, 50, 20]
    assert ctce._front(ys, c) <= c - 5
