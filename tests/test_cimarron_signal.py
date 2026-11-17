"""Contract of the clean-room Cimarron signal stage (cimarron_signal)."""
import os
from pathlib import Path

import numpy as np
import pytest

from cimarron_basecaller import cimarron_signal as cs
from cimarron_basecaller.rsd_io import read_rsd

ROOT = Path(__file__).resolve().parent.parent
M13 = ROOT / "example_data" / "M13"
CLEANROOM = Path(os.environ.get(
    "CIMARRON_CLEANROOM_DATA",
    "/tmp/claude-1000/-home-tv-Limoncello-scorer/1433b067-08ec-4601-aba0-99feee20fecf/scratchpad/cleanroom/data"))


def _rsd(name="A01.rsd"):
    path = M13 / name
    if not path.exists():
        pytest.skip(f"example_data/M13/{name} missing")
    return read_rsd(str(path))


def synthetic_raw(n=6000, seed=0):
    """Pre-run noise, then Gaussian bands every ~9 scans in random dyes, mixed
    through the PLAYER1 dye spectra, on a constant offset."""
    rng = np.random.default_rng(seed)
    mix = np.linalg.inv(cs.PLAYER1_SEPARATION)
    sig = np.zeros((n, 4))
    t = np.arange(n)
    p = 2100.0
    while p < n - 300:
        d = rng.integers(4)
        sig[:, d] += 1000 * np.exp(-0.5 * ((t - p) / 2.2) ** 2)
        p += 9 + rng.normal(0, 0.8)
    raw = sig @ mix.T + 300 + rng.normal(0, 3, size=(n, 4))
    return raw


def test_footer_metadata_and_dye_order():
    r = _rsd()
    meta = cs.read_meta(r.footer_raw)
    assert meta.get("MACHINE ID") == "PLAYER1"
    assert cs.dye_order_from_footer(r.footer_raw) == "TGCA"
    assert cs.read_meta(b"") == {}
    assert cs.dye_order_from_footer(b"") is None


def test_segments_cover_every_row_once():
    for L in (500, 2048, 2049, 3900, 7235, 9001):
        starts = cs.segment_starts(L)
        own = cs.segment_ownership(starts, L)
        assert own[0][0] == 0 and own[-1][1] == L
        for (a, b), (c, d) in zip(own, own[1:]):
            assert b == c and a < b
        for s, (a, b) in zip(starts, own):
            assert s <= a and b <= s + cs.SEG_LEN
        if L > cs.SEG_LEN:
            assert starts[-1] == L - cs.SEG_LEN


def test_homomorphic_deconvolution_keeps_band_positions():
    N = cs.SEG_LEN
    t = np.arange(N)
    x = np.exp(-0.5 * ((t - 700) / 3.0) ** 2) + 0.6 * np.exp(-0.5 * ((t - 1300) / 3.0) ** 2) + 0.01
    y = cs.homomorphic_deconvolve(x[:, None], 140.0)[:, 0]
    assert abs(int(np.argmax(y)) - 700) <= 1
    assert abs(int(np.argmax(y[1000:])) + 1000 - 1300) <= 1


def test_process_frame_and_flattening_relations():
    raw = synthetic_raw()
    res = cs.process(raw, "TGCA", {"MACHINE ID": "PLAYER1"})
    S, E = res.start_tp, res.end_tp
    assert 1 <= S < E <= len(raw)
    n = E - S
    assert res.trace.shape == (n, 4) and res.internal.shape == (n, 4) and res.divisor.shape == (n,)
    assert np.all(res.trace >= 0)
    assert np.all(res.divisor > 0)
    # exact relation between the internal and the flattened trace
    again = cs.flatten(res.internal, cs.flatten_divisor(res.internal), res.gain, res.offset)
    assert np.allclose(again, res.trace)
    # the caller trace is the flattened trace with the divisor multiplied back
    assert np.allclose(res.caller_trace, res.trace * res.divisor[:, None] / cs.FLATTEN_SCALE)
    # unflatten inverts flatten wherever neither clip was active
    back = cs.unflatten(res.trace, res.divisor, res.gain, res.offset)
    pos = (res.internal > 0) & (res.trace > 0)
    assert np.allclose(back[pos], res.internal[pos])
    assert res.estimated is False


def test_unknown_instrument_is_flagged_estimated():
    raw = synthetic_raw(seed=1)
    res = cs.process(raw, "TGCA", {"MACHINE ID": "SOMETHING ELSE"})
    assert res.estimated is True
    assert res.info["separation_source"] == "estimated"
    assert np.all(np.isfinite(res.trace))


def test_estimated_separation_recovers_dye_spectra():
    raw = synthetic_raw(seed=2)
    S, E = cs.read_span(raw)
    bs = cs.remove_baseline(raw[S - 1:E - 1])
    M = cs.estimate_separation(bs)
    # separating the mixed spectra should give a near-diagonal product
    prod = M @ np.linalg.inv(cs.PLAYER1_SEPARATION)
    prod = prod / np.diag(prod)[:, None]
    off = prod - np.diag(np.diag(prod))
    assert np.max(np.abs(off)) < 0.15


def test_block_delays_shift_whole_rows():
    x = np.zeros((1000, 4))
    x[500, :] = 1.0
    table = np.array([[1, 2, 3, 4]] * 20)
    y = cs.apply_block_delays(x, table)
    assert [int(np.argmax(y[:, c])) for c in range(4)] == [501, 502, 503, 504]


def test_real_trace_matches_cimarron_processed_trace():
    """On an m13 well with Cimarron's result file, the flattened trace
    correlates with Cimarron's stored processed trace (needs the clean-room data)."""
    rsd_path = CLEANROOM / "m13" / "rsd" / "A02.rsd"
    esd_path = CLEANROOM / "m13" / "esd" / "A02.esd"
    if not (rsd_path.exists() and esd_path.exists()):
        pytest.skip("clean-room m13 data not available")
    b = esd_path.read_bytes()
    # START/END TRACE PROCESSING of A02 (from its .esd tags)
    S, E = 2164, 9399
    proc = np.frombuffer(b[:(E - S) * 20], dtype=np.dtype([("cur", "<i4"), ("ch", "<f4", 4)]))["ch"]
    r = read_rsd(str(rsd_path))
    res = cs.process(r.trace, "TGCA", cs.read_meta(r.footer_raw), span=(S, E))
    best = []
    for c in range(4):
        best.append(max(np.corrcoef(np.roll(res.trace[:, c], k)[10:-10], proc[10:-10, c])[0, 1] for k in range(-5, 6)))
    assert min(best) > 0.85, best
