"""Five-dye ABIF (.fsa) runs: loading every channel and sizing a LIZ ladder on
channel 5.  The files are built by tests/abif_fixture.py to the ABIF layout."""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import fragment_sizing as fs
from analyzer_core import (SUPPORTED_EXTS, channel_dye, load_trace, n_channels,
                           physical_channel)
from abif_fixture import five_dye_fsa, short_array, write_abif

GS500 = fs.load_ladder("genescan500_liz")
SCAN = lambda bp: 1200.0 + 8.0 * bp + 0.004 * bp ** 2      # noqa: E731
TRUTH = [118.0, 207.5, 312.0]


@pytest.fixture
def fsa(tmp_path):
    path = tmp_path / "B07_run.fsa"
    five_dye_fsa(path, [SCAN(b) for b in TRUTH],
                 [SCAN(b) for b in GS500.lengths])
    return path


def test_fsa_is_a_supported_extension():
    assert ".fsa" in SUPPORTED_EXTS


def test_loads_all_five_channels_with_their_dyes(fsa):
    doc = load_trace(fsa)
    assert n_channels(doc) == doc.n_channels == 5
    assert doc.dyes == ["6-FAM", "VIC", "NED", "PET", "LIZ"]
    assert channel_dye(doc, 5) == "LIZ"
    assert doc.acgt.shape[1] == 4                    # the base-caller view
    assert doc.well == "B07"                         # TUBE1 "B7", normalised
    assert "5 dye channels" in doc.meta and "analyzed" in doc.meta
    liz = physical_channel(doc, 5)
    assert liz.max() > 500 and physical_channel(doc, 2).max() < 200
    with pytest.raises(ValueError):
        physical_channel(doc, 6)


def test_a_file_with_only_raw_data_still_loads(tmp_path):
    path = tmp_path / "raw.fsa"
    five_dye_fsa(path, [2000], [3000], raw_too=False)     # analyzed only
    assert load_trace(path).n_channels == 5
    path2 = tmp_path / "raw2.fsa"
    write_abif(path2, [short_array("DATA", k, np.full(500, 10.0 * k))
                       for k in (1, 2, 3, 4, 105)])
    doc = load_trace(path2)
    assert doc.n_channels == 5 and "raw data" in doc.meta
    assert physical_channel(doc, 5)[0] == pytest.approx(1050.0)


def test_sizes_against_a_liz_ladder_on_channel_5(fsa):
    res = fs.size_trace(load_trace(fsa), GS500, ladder_channel=5, sample_channel=1)
    assert not res.warnings
    assert 250.0 not in [b for _, b in res.anchors]      # vendor: not for sizing
    assert len(res.anchors) == len(GS500.lengths) - 1
    got = sorted(r["length_bp"] for r in res.rows if r["height"] > 300)
    assert got == pytest.approx(TRUTH, abs=0.3)
    assert all(r["dye"] == "6-FAM" and r["ladder_ok"] for r in res.rows)


def test_liz_ladders_default_to_channel_5():
    for key in ("genescan120_liz", "genescan500_liz", "genescan600_liz",
                "genescan1200_liz"):
        assert fs.load_ladder(key).default_channel == 5


def test_cli_sizes_a_fsa_with_the_ladders_default_channel(fsa, tmp_path):
    out = tmp_path / "sizes.csv"
    root = Path(__file__).resolve().parent.parent
    r = subprocess.run([sys.executable, str(root / "scorer.py"), "size", str(fsa),
                        "--ladder", "genescan500_liz", "--sample-channel", "1",
                        "--out", str(out)], capture_output=True, text=True, cwd=root)
    assert r.returncode == 0, r.stderr
    assert "anchors 15" in r.stdout and "!" not in r.stdout
    assert "6-FAM" in out.read_text(encoding="utf-8-sig")


# --------------------------------------------------------------------------- #
# a real ABIF file (Biopython's test data; not shipped here)
# --------------------------------------------------------------------------- #
import os                                                    # noqa: E402

PUBLIC = Path(os.environ.get("LIMONCELLO_ABIF_PUBLIC",
                             "/run/media/per/NVMe Disk/limoncello-data/abif_public"))


@pytest.mark.skipif(not (PUBLIC / "test.fsa").is_file(),
                    reason="public ABIF test files not present (LIMONCELLO_ABIF_PUBLIC)")
def test_a_real_fsa_loads_with_its_dyes_and_well():
    doc = load_trace(PUBLIC / "test.fsa")
    assert doc.n_channels == 4
    assert doc.dyes == ["5-FAM", "JOE", "NED", "ROX"]
    assert doc.well == "E01"
    # the ROX size standard: evenly sized bands across the run
    rox = physical_channel(doc, 4)
    peaks = fs.detect_peaks(rox, min_prominence_frac=0.015, min_height_frac=0.02)
    assert len([p for p in peaks if 150 < p.height < 300]) >= 15
