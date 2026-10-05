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


# ---------------------------------------------------------------------------
# Run conditions from the instrument's own text export
# ---------------------------------------------------------------------------
# Run Voltage / time / injection / temperature / PMT are NOT in the .rsd.
# They only exist in the MegaBACE raw-scan export, so these tests pin the two
# parsing traps: the inconsistent ": " separator and the free-text fields
# that continue onto following lines.

_EXPORT = """Run Info for the file : D:\\Data\\plate01\\A01.rsd
Sample name : A01
Plate name : plate01
Comment : Hel plate colo 829 DNA
Grad 50_65C
Inject 10KV, 60 sec, run 9kv, CTCE (53-50)x20
Chemistry name : ET Terminators
Base Caller : Cimarron 3.12
Run Voltage : 9
Run time : 60
Injection time : 25
Injection voltage : 10
Temperature : 53
PMT Voltage1: 750
PMT Voltage2: 750
Base order : TGCA
Number of lines : 12
Scan\tChannel1\tChannel2\tChannel3\tChannel4\tCurrent
0\t1\t2\t3\t4\t5
1\t6\t7\t8\t9\t10
"""


def _write_export(tmp_path, text, encoding="utf-8"):
    run = tmp_path / "plate01"
    (run / "Text").mkdir(parents=True)
    path = run / "Text" / "A01.txt"
    path.write_bytes(text.encode(encoding))
    return run


def test_run_export_parses_the_seven_headline_fields(tmp_path):
    run = _write_export(tmp_path, _EXPORT)
    f = analyzer_core.parse_run_export(run / "Text" / "A01.txt")
    assert f["Run Voltage"] == "9"
    assert f["Run time"] == "60"
    assert f["Injection time"] == "25"
    assert f["Injection voltage"] == "10"
    assert f["Temperature"] == "53"
    assert f["PMT Voltage1"] == "750"   # note: no space before the colon
    assert f["PMT Voltage2"] == "750"


def test_run_export_keeps_colon_in_windows_path(tmp_path):
    """'Run Info for the file : D:\\Data\\...' must split on the FIRST colon."""
    run = _write_export(tmp_path, _EXPORT)
    f = analyzer_core.parse_run_export(run / "Text" / "A01.txt")
    assert f["Run Info for the file"] == "D:\\Data\\plate01\\A01.rsd"


def test_run_export_joins_multiline_comment(tmp_path):
    run = _write_export(tmp_path, _EXPORT)
    f = analyzer_core.parse_run_export(run / "Text" / "A01.txt")
    assert "CTCE (53-50)x20" in f["Comment"]
    assert "Grad 50_65C" in f["Comment"]


def test_run_export_stops_before_the_numeric_table(tmp_path):
    """Data rows must not be mistaken for metadata keys."""
    run = _write_export(tmp_path, _EXPORT)
    f = analyzer_core.parse_run_export(run / "Text" / "A01.txt")
    assert "0" not in f
    assert "Channel1" not in f


@pytest.mark.parametrize("encoding", ["utf-8", "utf-16-le", "utf-16"])
def test_run_export_detects_encoding(tmp_path, encoding):
    """utf-16 without a BOM must still decode, not silently become mojibake."""
    run = _write_export(tmp_path, _EXPORT, encoding=encoding)
    f = analyzer_core.parse_run_export(run / "Text" / "A01.txt")
    assert f.get("Temperature") == "53"
    assert f.get("Run Voltage") == "9"


def test_find_run_export_and_load_run_params(tmp_path):
    run = _write_export(tmp_path, _EXPORT)
    found = analyzer_core.find_run_export(run, "A01")
    assert found is not None and found.name == "A01.txt"
    rp = analyzer_core.load_run_params(run / "A01.rsd")  # file need not exist
    assert rp["Temperature"] == "53"


def test_no_text_export_yields_empty_params(tmp_path):
    """The common case (OY/ runs have none) must degrade quietly."""
    run = tmp_path / "plate01"
    run.mkdir()
    assert analyzer_core.find_run_export(run, "A01") is None
    assert analyzer_core.load_run_params(run / "A01.rsd") == {}
    assert analyzer_core.run_params_summary({}) == []


def test_run_params_summary_puts_headline_fields_first(tmp_path):
    run = _write_export(tmp_path, _EXPORT)
    rows = analyzer_core.run_params_summary(
        analyzer_core.parse_run_export(run / "Text" / "A01.txt"))
    labels = [k for k, _ in rows]
    assert labels[:7] == list(analyzer_core.RUN_CONDITION_FIELDS)
    assert "Run Voltage" in labels and "Temperature" in labels


# ---------------------------------------------------------------------------
# Run conditions decoded from the .rsd binary
# ---------------------------------------------------------------------------
# The settings ARE in the .rsd: after the TLV text footer sits a fixed-layout
# binary block whose first five floats are run voltage, run time, injection
# time, injection voltage and temperature, with the two PMT voltages 24 bytes
# later. These tests pin both the decoding and the range validation that keeps
# arbitrary bit patterns in trace data from being reported as settings.

def _rsd_with_settings(values=(9.0, 60.0, 25.0, 10.0, 53.0),
                       pmt=(750.0, 750.0), pad_before=b""):
    import struct as _struct
    block = b"".join(_struct.pack("<f", v) for v in values)
    block += b"\x02\x00\x00\x00"          # 4-byte tag, as in real files
    block += b"".join(_struct.pack("<f", v) for v in pmt)
    return pad_before + block


def test_rsd_run_settings_decodes_all_seven_fields():
    got = analyzer_core.parse_rsd_run_settings(_rsd_with_settings())
    assert got["Run Voltage"] == "9"
    assert got["Run time"] == "60"
    assert got["Injection time"] == "25"
    assert got["Injection voltage"] == "10"
    assert got["Temperature"] == "53"
    assert got["PMT Voltage1"] == "750"
    assert got["PMT Voltage2"] == "750"


def test_rsd_run_settings_found_at_any_offset():
    """The block is located by content, not by a hard-coded file offset."""
    for pad in (b"", b"\x00" * 3, b"\x00" * 17, b"abc\x00\x00"):
        got = analyzer_core.parse_rsd_run_settings(_rsd_with_settings(pad_before=pad))
        assert got.get("Temperature") == "53", pad


def test_rsd_run_settings_rejects_implausible_values():
    """Trace data can hold any bit pattern; bogus settings must be refused."""
    # A run cannot be 3400 kV, or 900 s at 0 V, or 900 C.
    assert analyzer_core.parse_rsd_run_settings(
        _rsd_with_settings(values=(3400.0, 60.0, 25.0, 10.0, 53.0))) == {}
    assert analyzer_core.parse_rsd_run_settings(
        _rsd_with_settings(values=(9.0, 60.0, 25.0, 0.0, 53.0))) == {}
    assert analyzer_core.parse_rsd_run_settings(
        _rsd_with_settings(values=(9.0, 60.0, 25.0, 10.0, 900.0))) == {}
    # Plausible voltages but nonsense PMT must not be accepted either.
    assert analyzer_core.parse_rsd_run_settings(
        _rsd_with_settings(pmt=(3.0, 0.5))) == {}


def test_rsd_run_settings_ignores_nan_and_empty():
    import struct as _struct
    nan = b"".join(_struct.pack("<f", float("nan")) for _ in range(5))
    assert analyzer_core.parse_rsd_run_settings(nan) == {}
    assert analyzer_core.parse_rsd_run_settings(b"") == {}
    assert analyzer_core.parse_rsd_run_settings(b"\x00" * 8) == {}


def test_rsd_run_settings_allows_long_injection_time():
    """One real plate used a 220 s injection -- the bound must not be tight."""
    got = analyzer_core.parse_rsd_run_settings(
        _rsd_with_settings(values=(9.0, 75.0, 220.0, 10.0, 53.0)))
    assert got["Injection time"] == "220"


def test_real_rsd_files_expose_run_conditions():
    """Guards the decoder against the actual files on disk."""
    cases = [
        ("/media/per/SIDS/mt_nucl/130323_mt_nucl_#8_grad_50_65CRun01/A01.rsd",
         {"Run Voltage": "9", "Run time": "60", "Temperature": "53",
          "PMT Voltage1": "750"}),
        ("/media/per/78B0C7DE1FA7081C/OY/OY_rs1695_N2_210910Run01/A04.rsd",
         {"Run Voltage": "9", "Run time": "75", "Temperature": "56"}),
    ]
    for spec, want in cases:
        p = Path(spec)
        if not p.exists():
            pytest.skip(f"sample data not present: {spec}")
        doc = load_trace(p)
        for key, val in want.items():
            assert doc.run_params.get(key) == val, f"{p.name}:{key}"
