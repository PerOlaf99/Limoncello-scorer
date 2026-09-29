"""Loader and basecaller tests.

These exercise the real example traces that ship with the repo, so they also
guard the vendored basecaller against accidental breakage."""
from pathlib import Path

import pytest

import analyzer_core
from analyzer_core import AnalysisSettings, load_trace, run_basecall

ROOT = Path(__file__).resolve().parent.parent
M13 = ROOT / "example_data" / "M13"


pytestmark = pytest.mark.skipif(
    not analyzer_core.BASECALLER_AVAILABLE,
    reason="vendored basecaller not importable",
)


def _first_rsd() -> Path:
    files = sorted(M13.glob("*.rsd"))
    if not files:
        pytest.skip("example_data/M13 has no .rsd files")
    return files[0]


def test_load_trace_shape():
    doc = load_trace(_first_rsd())
    assert doc.acgt.ndim == 2 and doc.acgt.shape[1] == 4
    assert doc.acgt.shape[0] == doc.n_scans > 0
    assert doc.source == "rsd"


@pytest.mark.parametrize("preset", ["mb1000_accuracy", "mb1000_length", "raw_peaks"])
def test_run_basecall_calls_bases(preset):
    doc = load_trace(_first_rsd())
    run_basecall(doc, AnalysisSettings(basecaller=preset))
    assert len(doc.sequence) == len(doc.qualities)
    assert len(doc.sequence) > 100
    assert set(doc.sequence) <= set("ACGTN")
    assert len(doc.peak_positions) == len(doc.sequence)


def test_unknown_preset_falls_back_with_warning():
    """An unknown preset must not crash: it warns and uses the UI knobs."""
    doc = load_trace(_first_rsd())
    with pytest.warns(RuntimeWarning, match="Unknown basecaller preset"):
        run_basecall(doc, AnalysisSettings(basecaller="does-not-exist"))
    assert len(doc.sequence) > 0


def test_display_trace_modes():
    doc = load_trace(_first_rsd())
    raw = analyzer_core.display_trace(doc, AnalysisSettings(view_mode="raw"))
    base = analyzer_core.display_trace(doc, AnalysisSettings(view_mode="baseline"))
    assert raw.shape == doc.raw.shape
    assert base.shape == doc.acgt.shape
