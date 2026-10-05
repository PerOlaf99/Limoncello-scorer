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


# --------------------------------------------------------------------------- #
# RSD run fields (Comments -> Run info)
# --------------------------------------------------------------------------- #
def _footer(label: str, value: str) -> bytes:
    """Build one ``<label> NUL 0x05 <len> <value>`` footer field."""
    val = value.encode("latin-1") + b"\x00"
    return label.encode("latin-1") + b"\x00\x05" + bytes([len(val)]) + val


def test_parse_rsd_footer_reads_run_fields():
    raw = (_footer("R CODE", "run7_080801") + _footer("SAMPLE NAME", "A01")
           + _footer("WELL ID", "A01") + _footer("MACHINE ID", "GRYPHON2")
           + _footer("BASE CALLER", "Cimarron 1.53"))
    fields = analyzer_core.parse_rsd_footer(raw)
    assert fields["R CODE"] == "run7_080801"
    assert fields["SAMPLE NAME"] == "A01"
    assert fields["WELL ID"] == "A01"
    assert fields["MACHINE ID"] == "GRYPHON2"
    assert fields["BASE CALLER"] == "Cimarron 1.53"


def test_parse_rsd_footer_ignores_padding_between_fields():
    """Inter-field type/length bytes must not be glued onto the next label."""
    raw = (_footer("PLATE ID", "P1") + b"\x01\x12" + _footer("PLATE ID", "P1")
           + b"\x05\x0c" + _footer("BASE CALLER", "Cimarron"))
    fields = analyzer_core.parse_rsd_footer(raw)
    assert "PLATE ID" in fields
    assert "BASE CALLER" in fields
    for key in fields:
        assert key == key.lstrip(chr(0) or "").strip("\t\r\n ")
        assert key.encode("latin-1").decode("latin-1") == key


def test_parse_rsd_footer_tolerates_missing_or_garbage():
    assert analyzer_core.parse_rsd_footer(b"") == {}
    assert analyzer_core.parse_rsd_footer(None) == {}
    assert analyzer_core.parse_rsd_footer(b"no structure here at all") == {}
    # Truncated length byte must not raise.
    assert isinstance(analyzer_core.parse_rsd_footer(b"X\x00\x05\x40ab"), dict)


def test_footer_summary_is_readable_not_raw_bytes():
    """Run info used to show b'\\x00\\x00\\x00\\x00' and a truncated footer."""
    raw = (_footer("R CODE", "run7") + _footer("SAMPLE NAME", "A01")
           + _footer("WELL ID", "A01"))
    summary = analyzer_core.rsd_footer_summary(analyzer_core.parse_rsd_footer(raw))
    assert summary.startswith("RSD |")
    assert "SAMPLE NAME=A01" in summary
    assert "\\x" not in summary and "b'" not in summary


def test_loaded_rsd_exposes_decoded_run_fields():
    doc = load_trace(_first_rsd())
    if doc.source != "rsd":
        pytest.skip("not an RSD file")
    fields = doc.footer_fields
    # A real RSD records at least its run and well.
    assert fields.get("WELL ID") or fields.get("SAMPLE NAME")
    assert doc.well == (fields.get("WELL ID") or fields.get("SAMPLE NAME")
                        or doc.path.stem)
    assert "\\x" not in doc.meta
