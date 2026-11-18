"""Dye-overlap correction, and sizing with it: pull-up is removed, the peak
table carries standards, quality and flags; and on the GP 2.2 demo plate the
result is checked against Genetic Profiler 2.2's own exported peak list."""
import csv
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest

import fragment_sizing as fs
import spectral as sp

# a MegaBACE-like overlap: HEX (dye 4) reads at 109 % in channel 3
M_TRUE = np.array([[1.0, 0.164, 0.306, 0.237],
                   [0.077, 1.0, 0.001, 0.006],
                   [0.080, 0.020, 1.0, 1.092],
                   [0.002, 0.020, 0.189, 1.0]])
N = 6000
SCAN = lambda bp: 1500.0 + 9.0 * bp                       # noqa: E731
ET400 = list(fs.load_ladder("et400r").lengths)


def gauss(centres, h=1000.0, w=4.0):
    x = np.arange(N)
    y = np.zeros(N)
    for c in centres:
        y += h * np.exp(-0.5 * ((x - c) / w) ** 2)
    return y


def measured(true_dyes, seed=0):
    """(scans, 4) as the detector sees per-dye signals mixed by M_TRUE."""
    rng = np.random.default_rng(seed)
    obs = np.column_stack(true_dyes) @ M_TRUE.T
    return obs + 30.0 + rng.normal(0, 2.0, obs.shape)


class Doc:
    def __init__(self, channels, well="B01"):
        self.channels = channels
        self.acgt = channels
        self.well, self.path, self.base_order = well, f"{well}.rsd", "ACGT"
        self.dyes = ["ET-ROX", "FAM", "NED", "HEX"]


def test_reads_a_genetic_profiler_xtalk_tdb(tmp_path):
    p = tmp_path / "xtalk.tdb"
    p.write_text("auto\r\n1 0 0 0 \r\n0 1 0 0 \r\n0 0 1 0 \r\n0 0 0 1 \r\nDemo\r\n"
                 + "".join(" ".join(f"{v:.4f}" for v in row) + " \r\n" for row in M_TRUE),
                 encoding="latin-1")
    mats = sp.read_xtalk_tdb(p)
    assert set(mats) == {"auto", "Demo"}
    assert mats["auto"].is_identity()
    assert np.allclose(mats["Demo"].m, M_TRUE, atol=1e-4)


def test_matrix_wells_are_found_and_assigned_even_when_hex_beats_ned(tmp_path):
    z = np.zeros(N)
    peaks = gauss([2000, 3000, 4000])
    wells = {"A02": measured([peaks, z, z, z]), "A03": measured([peaks, z, z, z], 1),
             "A04": measured([z, peaks, z, z]), "A05": measured([z, peaks, z, z], 1),
             "A06": measured([z, z, peaks, z]), "A07": measured([z, z, peaks, z], 1),
             "A08": measured([z, z, z, peaks]), "A09": measured([z, z, z, peaks], 1),
             # a sample well: ladder + independent fragments in three dyes
             "B01": measured([gauss([SCAN(b) for b in ET400]), gauss([2600]),
                              gauss([3100]), gauss([3700])])}
    found = sp.find_matrix_wells(wells)
    assert found == {0: ["A02", "A03"], 1: ["A04", "A05"], 2: ["A06", "A07"],
                     3: ["A08", "A09"]}
    m = sp.plate_matrix(tmp_path, wells, prefer="plate")
    assert np.allclose(m.m, M_TRUE, atol=0.02)
    assert "A08" in m.source


def test_correction_removes_pullup_and_sizing_ignores_it():
    hex_peak = gauss([SCAN(150.0)], h=4000.0)               # strong HEX fragment
    fam = gauss([SCAN(203.0)])
    obs = measured([gauss([SCAN(b) for b in ET400]), fam, np.zeros(N), hex_peak])
    lad = fs.load_ladder("et400r")
    raw = fs.size_trace(Doc(obs), lad, ladder_channel=1, sample_channels=[2, 3, 4],
                        base_order="ACGT")
    ned_raw = [r for r in raw.rows if r["channel"] == 3 and r["height"] > 1000]
    assert ned_raw and abs(ned_raw[0]["length_bp"] - 150.0) < 0.5   # pull-up sized
    assert "possible pull-up" not in ned_raw[0]["flags"]            # tallest there
    m = sp.SpectralMatrix("t", M_TRUE)
    res = fs.size_trace(Doc(obs), lad, ladder_channel=1, sample_channels=[2, 3, 4],
                        base_order="ACGT", matrix=m)
    assert not [r for r in res.rows if r["channel"] == 3 and r["height"] > 200]
    sizes = {r["channel"]: r["length_bp"] for r in res.rows if r["quality"] >= 50}
    assert sizes[2] == pytest.approx(203.0, abs=0.2)
    assert sizes[4] == pytest.approx(150.0, abs=0.2)
    assert res.matrix == "t"


def test_table_has_standards_and_scored_peaks():
    obs = measured([gauss([SCAN(b) for b in ET400]), gauss([SCAN(118.0), SCAN(250.0)]),
                    np.zeros(N), np.zeros(N)])
    res = fs.size_trace(Doc(obs), "et400r", ladder_channel=1, sample_channels=[2],
                        base_order="ACGT", matrix=sp.SpectralMatrix("t", M_TRUE))
    rows = fs.table_rows([res])
    std = [r for r in rows if r["kind"] == "standard"]
    smp = [r for r in rows if r["kind"] == "sample"]
    assert [r["band_bp"] for r in std] == ET400
    assert all(r["scan"] == pytest.approx(SCAN(r["band_bp"]), abs=1) for r in std)
    assert all(abs(r["residual_bp"]) < 0.3 for r in std[1:-1])
    assert sorted(r["length_bp"] for r in smp) == pytest.approx([118.0, 250.0], abs=0.2)
    assert all(r["quality"] >= 80 and r["flags"] == "" for r in smp)
    assert {"width", "area", "snr", "quality", "flags", "dye"} <= set(smp[0])


def test_quality_score_reasons():
    q, f = fs.peak_quality(snr=200, width_ratio=1.0, pullup_ratio=0.1, in_range=True,
                           flat=False, ladder_rel=0.01)
    assert q >= 90 and f == []
    q2, f2 = fs.peak_quality(snr=200, width_ratio=3.0, pullup_ratio=5.0, in_range=False,
                             flat=True, ladder_rel=0.01)
    assert q2 < 20
    assert {"broad", "possible pull-up", "outside ladder", "flat top (saturated)"} <= set(f2)


# --------------------------------------------------------------------------- #
# against Genetic Profiler 2.2's own analysis (GP20 demo plate, run under Wine)
# --------------------------------------------------------------------------- #
DEMO = Path(os.environ.get("LIMONCELLO_GP22_DEMO",
                           "/run/media/per/NVMe Disk/limoncello-data/"
                           "fragment_sizing/GP22_demo"))
GP_REF = DEMO / "gp_reference" / "GP20_genericpeaks_ET400R_Demo-matrix_Local3rd.txt"


@pytest.mark.skipif(not GP_REF.is_file(),
                    reason="Genetic Profiler reference export not present")
@pytest.mark.parametrize("well", ["A05", "C07", "E10", "H03"])
def test_matches_genetic_profiler_on_gp20(well):
    from analyzer_core import load_trace
    gp = defaultdict(list)
    with open(GP_REF, encoding="latin-1") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["Well ID"] == well:
                gp[int(r["Channel"])].append((float(r["Scan No."]), float(r["Size"])))
    m = sp.run_matrices(DEMO / "GP20_Demo_Data")["Demo"]
    doc = load_trace(DEMO / "GP20_Demo_Data" / f"{well}.rsd", base_order="TGCA")
    res = fs.size_trace(doc, "et400r", ladder_channel=1, sample_channels=[2, 3, 4],
                        base_order="TGCA", matrix=m)
    assert not res.warnings
    ours = {r["band_bp"]: r["scan"] for r in res.standards}
    assert all(abs(ours[b] - s) <= 1.5 for s, b in gp[1] if b in ours)
    curve = fs.SizeCurve([s for s, _ in res.anchors], [b for _, b in res.anchors])
    for ch in (2, 3, 4):
        for s, size in gp[ch]:
            if 60 <= size <= 400:
                assert float(curve.bp([s])[0]) == pytest.approx(size, abs=0.5)
    # GP's sample peaks are found: at least 85 % of them, by position
    hit = tot = 0
    for ch in (2, 3, 4):
        mine = [r["scan"] for r in res.rows if r["channel"] == ch]
        for s, size in gp[ch]:
            if size >= 30:
                tot += 1
                hit += any(abs(x - s) <= 2 for x in mine)
    assert hit >= 0.85 * tot
