"""Genotyping tests: peak-picker clustering and table export.

The genotyping module imports tkinter + matplotlib at module load (the picker
engine is headless but shares the file), so these skip cleanly where a GUI
toolkit is unavailable."""
from pathlib import Path

import pytest

pytest.importorskip("tkinter")
pytest.importorskip("matplotlib")

import genotyping as g  # noqa: E402
from analyzer_core import load_trace  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
M13 = ROOT / "example_data" / "M13"


def _rec(scan, col, area=100.0, height=1.0, kind="main", gid=1):
    return g._Record(file="x.rsd", well="A01", scan=scan, channel=col + 1,
                     base=g.CHANNEL_ORDER[col], kind=kind, height=height,
                     area=area, left=scan - 1, right=scan + 1, onset=scan - 1,
                     end=scan + 1, color="#000000", col=col, gid=gid)


def _picker():
    files = sorted(M13.glob("*.rsd"))
    if not files:
        pytest.skip("example_data/M13 has no .rsd files")
    doc = load_trace(files[0])
    return g.PeakPicker(doc, files[0])


def test_save_table_empty_does_not_crash(tmp_path):
    g.save_table(tmp_path / "empty.csv", [])
    g.save_table(tmp_path / "empty.json", [])
    assert (tmp_path / "empty.csv").exists()
    assert (tmp_path / "empty.json").exists()


def test_save_table_roundtrip(tmp_path):
    rows = [{"file": "x", "well": "A01", "scan": 10, "channel": 1}]
    g.save_table(tmp_path / "peaks.csv", rows)
    text = (tmp_path / "peaks.csv").read_text(encoding="utf-8-sig")
    assert "A01" in text and "scan" in text


def test_clusters_same_channel_window():
    pk = _picker()
    pk.records = [_rec(100, 0, area=100), _rec(102, 0, area=100),
                  _rec(100, 1, area=50), _rec(500, 0, area=100)]
    clusters = pk.clusters()
    assert sorted(len(c) for c in clusters) == [1, 1, 2]


def test_clust_frac_and_duplex():
    pk = _picker()
    a, b = _rec(100, 0, area=100, gid=1), _rec(102, 0, area=20, gid=2)
    pk.records = [a, b]
    assert pk.clust_frac(b) == pytest.approx(20 / 120)


def test_mark_std_labels_standard_peaks():
    pk = _picker()
    pk.records = [_rec(s, 0, gid=i) for i, s in enumerate((10, 30, 50, 70))]
    msg = pk.mark_std(length_bp=120)
    assert "HOM1" in msg and "HET2" in msg
    assert pk.duplex_of(pk.records[0]) == "HOM1"
    assert pk.duplex_of(_rec(999, 0)) == ""


def test_export_rows_has_expected_columns():
    pk = _picker()
    pk.records = [_rec(100, 0, area=100, gid=1)]
    row = pk.export_rows()[0]
    for key in ("file", "well", "scan", "channel", "base", "kind",
                "height_V", "area_Vscan", "duplex", "fraction"):
        assert key in row
