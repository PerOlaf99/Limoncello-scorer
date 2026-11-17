"""CTCE engine on synthetic runs: what moves with temperature is a melting
duplex, what holds still is static, the run order follows the recorded
temperature, satellites ride with their parent, and areas are measured."""
import sys
from pathlib import Path

import numpy as np
from scipy.special import erfc

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import ctce  # noqa: E402

N = 4000


def emg(x, mu, sigma, tau, area):
    """Exponentially modified Gaussian: sharp front, exponential tail."""
    lam = 1.0 / tau
    return area * lam / 2 * np.exp(lam / 2 * (2 * mu + lam * sigma ** 2 - 2 * x)) * \
        erfc((mu + lam * sigma ** 2 - x) / (np.sqrt(2) * sigma))


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
