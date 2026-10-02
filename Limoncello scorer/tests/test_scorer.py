"""Scorer tests: allele grouping, genotype calls, QC report and I/O.

Pure logic -- no tkinter, no trace files needed."""
import json

import pytest

import scorer as s


def main_row(well, scan, channel, area, duplex="", file="x.rsd", kind="main"):
    return {"file": file, "well": well, "scan": scan, "channel": channel,
            "base": "A", "kind": kind, "area_Vscan": area, "height_V": 1.0,
            "duplex": duplex, "length_bp": 300, "fraction": ""}


def test_median_peak_spacing():
    rows = [main_row("A01", 100, 1, 1), main_row("A01", 110, 1, 1),
            main_row("A01", 125, 1, 1)]
    assert s.median_peak_spacing(rows) == pytest.approx(12.5)


def test_group_alleles_same_channel_window():
    rows = [main_row("A01", 100, 1, 100),
            main_row("A01", 103, 1, 100),   # same channel, inside window
            main_row("A01", 400, 1, 100),   # far away -> own cluster
            main_row("A01", 100, 2, 50)]    # different channel -> own cluster
    clusters = s.group_alleles(rows, spacing=10.0)
    sizes = sorted(len(c) for c in clusters)
    assert sizes == [1, 1, 2]


def test_heterozygote_fraction_and_confidence():
    cl = [main_row("A01", 100, 1, 100), main_row("A01", 103, 1, 100)]
    assert s.genotype_call(cl) == "het"
    assert s.minor_fraction(cl) == pytest.approx(0.5)
    assert s.cluster_confidence(cl, "het") == 100.0


def test_homozygote_from_duplex_label():
    cl = [main_row("A01", 100, 1, 500, duplex="HOM1"),
          main_row("A01", 101, 1, 5)]      # second peak is noise
    assert s.genotype_call(cl) == "hom-major"
    assert s.cluster_confidence(cl, "hom-major") == 80.0


def test_single_peak_is_no_call():
    cl = [main_row("A01", 100, 1, 500)]
    assert s.genotype_call(cl) == "no-call"
    assert s.cluster_confidence(cl, "no-call") == 0.0


def test_boundary_fraction_is_het():
    cl = [main_row("A01", 100, 1, 100), main_row("A01", 103, 1, 20)]
    assert s.minor_fraction(cl) == pytest.approx(20 / 120)
    assert s.genotype_call(cl) == "het"          # 0.167 is inside 0.15..0.85


def test_score_rows_handles_plain_column_names():
    rows = [{"file": "x", "well": "A", "scan": 1, "channel": 1, "base": "A",
             "kind": "main", "area": 100, "height": 2.0},
            {"file": "x", "well": "A", "scan": 4, "channel": 1, "base": "A",
             "kind": "main", "area": 100, "height": 2.0}]
    scored = s.score_rows(rows)
    assert all(r["call"] == "het" for r in scored)


def test_satellite_inherits_call():
    rows = [main_row("A01", 100, 1, 100), main_row("A01", 103, 1, 100),
            main_row("A01", 106, 1, 4, kind="plusA")]
    scored = {r["kind"]: r for r in s.score_rows(rows)}
    assert scored["plusA"]["call"] == "satellite-het"


def test_well_summary():
    rows = [main_row("A01", 100, 1, 100), main_row("A01", 103, 1, 100),
            main_row("A02", 100, 1, 500, duplex="HOM1")]
    summary = {w["well"]: w for w in s.well_summary(rows)}
    assert summary["A01"]["genotype"] == "het"
    assert summary["A02"]["genotype"] == "hom-major"


def test_html_report_escapes_and_renders():
    qc = [{"well": "<script>", "source": "rsd", "length": 10, "q_mean": 1.0,
           "q_min": 0.1, "q_10pct": 0.2, "low_signal_frac": 0.0, "n_bases": 0,
           "mean_spacing": 12.0, "spacing_cv": 0.1, "preset": "p", "pass": True}]
    html = s.html_report(qc_rows=qc, title="t")
    assert "<script>" not in html and "&lt;script&gt;" in html
    assert "PASS" in html


def test_read_write_rows_roundtrip(tmp_path):
    rows = [main_row("A01", 100, 1, 100, duplex="HET1")]
    csv_path = tmp_path / "rows.csv"
    s.write_rows(str(csv_path), rows)
    assert s.read_rows(str(csv_path))[0]["well"] == "A01"

    json_path = tmp_path / "rows.json"
    s.write_rows(str(json_path), rows)
    assert json.loads(json_path.read_text())[0]["scan"] == 100


def test_write_rows_empty(tmp_path):
    csv_path = tmp_path / "empty.csv"
    s.write_rows(str(csv_path), [])
    assert csv_path.exists()
    json_path = tmp_path / "empty.json"
    s.write_rows(str(json_path), [])
    assert json.loads(json_path.read_text()) == []
