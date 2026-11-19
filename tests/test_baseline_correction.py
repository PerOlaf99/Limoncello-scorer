"""Baseline correction of the display: the slow background goes, peaks stay,
and a drag measures exactly the corrected trace that is shown."""
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import analyzer_core as ac  # noqa: E402
from genotyping import DragAreaPicker  # noqa: E402


def _trace():
    x = np.arange(3000, dtype=float)
    hump = 800.0 * np.exp(-0.5 * ((x - 1500) / 400) ** 2) + 200.0     # broad background
    peak = 1000.0 * np.exp(-0.5 * ((x - 1500) / 3.0) ** 2)            # area 1000*3*sqrt(2pi)
    return x, hump, peak


def test_background_goes_and_the_peak_stays():
    _x, hump, peak = _trace()
    y = ac.baseline_correct(hump + peak, 101)
    true = 1000 * 3 * np.sqrt(2 * np.pi)
    # the opening sits a little below a curved hump's crest, so a small peak
    # right on top of a big hump keeps a few % of it (here ~7 %)
    assert abs(y[1440:1561].sum() - true) / true < 0.08
    assert np.abs(y[:1300]).max() < 20 and np.abs(y[1700:]).max() < 20
    assert (y >= 0).all()


def test_display_trace_follows_the_setting():
    _x, hump, peak = _trace()
    acgt = np.column_stack([hump, hump + peak, hump, hump])
    doc = SimpleNamespace(acgt=acgt, raw=acgt)
    off = ac.display_trace(doc, ac.AnalysisSettings())
    on = ac.display_trace(doc, ac.AnalysisSettings(baseline_correct=True))
    assert off is acgt
    assert on[:, 1].max() < off[:, 1].max() and np.abs(on[:1300]).max() < 20


def test_a_drag_measures_the_corrected_trace_and_is_measured_again(tmp_path):
    _x, hump, peak = _trace()
    acgt = np.column_stack([hump, hump + peak, hump, hump])
    doc = SimpleNamespace(acgt=acgt, well="A01")
    state = {"on": False}

    def trace():
        return (ac.baseline_correct(acgt, 101), "window 101") if state["on"] else (acgt, "none")

    pk = DragAreaPicker(doc, tmp_path / "A01.rsd", trace=trace)
    rec = pk.add(1440, 1560, col=1)      # the drag's ends sit high on the hump
    raw_area = rec["area_Vscan"]
    assert rec["baseline_correction"] == "none"
    state["on"] = True
    assert pk.remeasure() == 1
    true = 1000 * 3 * np.sqrt(2 * np.pi)
    assert rec["baseline_correction"] == "window 101"
    assert abs(rec["area_Vscan"] - true) / true < 0.08
    assert rec["area_Vscan"] != raw_area
