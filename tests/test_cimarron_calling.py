"""Contract of the clean-room Cimarron calling stage (cimarron_calling)."""
import os
import re
import struct
from pathlib import Path

import numpy as np
import pytest

from cimarron_basecaller import cimarron_calling as cc

CLEANROOM = Path(os.environ.get(
    "CIMARRON_CLEANROOM_DATA",
    "/tmp/claude-1000/-home-tv-Limoncello-scorer/1433b067-08ec-4601-aba0-99feee20fecf/scratchpad/cleanroom/data"))


def band_trace(letters, spacing=9.0, width=2.0, start=100, tail=200, order="TGCA", heights=None):
    n = int(start + spacing * len(letters) + tail)
    t = np.arange(n)
    tr = np.zeros((n, 4))
    centres = []
    for i, b in enumerate(letters):
        c = start + spacing * i
        h = 1.0 if heights is None else heights[i]
        tr[:, order.index(b)] += h * np.exp(-0.5 * ((t - c) / width) ** 2)
        centres.append(int(round(c)))
    return tr, np.array(centres)


def test_band_detection_plateau_and_edges():
    E = np.array([0, 1, 2, 2, 1, 0, 1, 3, 3, 3, 2, 1, 0, 0.5, 1.0])
    b = cc.detect_bands(E)
    # the peak is the LAST scan of a plateau; the first peak (3) has no trough
    # on its left and the trailing rise has no peak, so one band remains
    assert list(b.p) == [9] and list(b.l) == [5] and list(b.r) == [12]
    E2 = np.array([1.0, 2, 1, 0, 2, 0, 3, 1, 0])
    b2 = cc.detect_bands(E2)
    assert list(b2.p) == [4] and list(b2.l) == [3] and list(b2.r) == [5]


def test_median_gap_rule():
    assert cc.median_int([3, 9, 5]) == 5
    assert cc.median_int([4, 7]) == 5          # (4 + 7) // 2
    assert cc.median_int([8, 9, 10, 11]) == 9


def test_expected_curve_rounds_and_holds_vertex():
    p = np.arange(0, 2000, 10.0)
    y = np.round(8.5 - 1e-3 * p + 4e-7 * p ** 2)   # convex, vertex beyond the last band
    e = cc.expected_curve(p, y)
    assert np.all(e == np.floor(e))
    assert len(set(e.tolist())) == 1               # held to the last band's value


def test_quality_single_rule_centroids():
    """SPEC_calling §13.3: one rule alone gives these exact qualities."""
    def q(rule, other=None, s=1.0):
        st = np.zeros((1, 9)); st[0, rule - 1] = s
        if other:
            st[0, other - 1] = s
        return 100 * cc.QUAL_DEFUZZ(st)[0]
    assert q(9) == pytest.approx(99.0, abs=0.02)
    assert q(8) == pytest.approx(87.167, abs=0.02)
    assert q(6) == pytest.approx(62.5, abs=0.02)
    assert q(5) == pytest.approx(50.0, abs=0.02)
    assert q(5, 6) == pytest.approx(56.25, abs=0.02)


def test_stitch_fallback_never_discards_the_read():
    w1 = cc.WindowResult(np.arange(100, 2000, 10), np.zeros(190, int), np.full(190, 99.0), 9.0, np.zeros(190, bool))
    # second window's bands are all 5 scans off the first window's: no join within 4 scans
    p2 = np.arange(1905, 3900, 10)
    w2 = cc.WindowResult(p2, np.ones(len(p2), int), np.full(len(p2), 99.0), 9.0, np.zeros(len(p2), bool))
    p, let, q, ins, nfb = cc.stitch([w1, None, w2])
    assert nfb == 1
    assert p[-1] == p2[-1]
    assert np.all(np.diff(p) > 3)
    # a join at distance 0 keeps earlier bands and continues with the new window
    p3 = np.arange(1900, 3900, 10)
    w3 = cc.WindowResult(p3, np.ones(len(p3), int), np.full(len(p3), 99.0), 9.0, np.zeros(len(p3), bool))
    p, let, q, ins, nfb = cc.stitch([w1, w3])
    assert nfb == 0 and np.array_equal(p, np.arange(100, 3900, 10))
    assert let[np.searchsorted(p, 1940)] == 1


def test_call_trace_reads_clean_bands():
    rng = np.random.default_rng(0)
    letters = "".join(rng.choice(list("ACGT"), 400))
    tr, centres = band_trace(letters)
    res = cc.call_trace(tr, "TGCA", start_tp=1000)
    # the first band has no left trough and the end margin drops the last
    # bands; everything in between must be exact
    k = len(res.sequence)
    assert k > 350
    assert res.sequence == letters[1:k + 1]
    assert np.array_equal(res.positions, centres[1:k + 1] + 1000 + 1)
    assert np.all((res.qualities >= 0) & (res.qualities <= 100))
    assert 0 <= res.good_start <= res.good_end < k


def test_two_channels_at_80_percent_give_n():
    tr, centres = band_trace("ACGTACGTAC" * 30)
    j = 150
    other = "TGCA".index("C") if "ACGTACGTAC"[j % 10] != "C" else "TGCA".index("A")
    t = np.arange(len(tr))
    tr[:, other] += 0.9 * np.exp(-0.5 * ((t - centres[j]) / 2.0) ** 2)
    res = cc.call_trace(tr, "TGCA", start_tp=0)
    i = int(np.argmin(np.abs(res.positions - 1 - centres[j])))
    assert res.sequence[i] == "N"
    assert res.sequence.count("N") == 1


def test_good_limits_whole_read():
    a, b = cc.good_limits(np.full(300, 99.0))
    assert (a, b) == (0, 298)
    a, b = cc.good_limits(np.r_[np.full(50, 10.0), np.full(200, 99.0), np.full(50, 10.0)])
    assert 40 <= a <= 52 and 245 <= b <= 252


def _esd(path):
    b = path.read_bytes()
    tags = {}
    for m in re.finditer(rb"\x05(.)([A-Z][A-Z0-9 /_\-\(\)\.#]{1,40})\x00", b, re.S):
        if m.group(1)[0] != len(m.group(2)) + 1:
            continue
        j = m.end()
        t = b[j]
        if t in (0xb1, 0xb2):
            n, k = (b[j + 1], j + 2) if t & 1 else (struct.unpack_from("<H", b, j + 1)[0], j + 3)
            v = np.frombuffer(b[k:k + n], "<i4")
            tags.setdefault(m.group(2).decode(), v[0] if v.size == 1 else v)
        elif t in (0x05, 0x06):
            n, k = (b[j + 1], j + 2) if t == 0x05 else (struct.unpack_from("<H", b, j + 1)[0], j + 3)
            tags.setdefault(m.group(2).decode(), b[k:k + n].rstrip(b"\0").decode("latin-1"))
    return b, tags


def test_calls_on_cimarrons_own_trace_agree():
    """On Cimarron's own stored trace (flattening inverted) the calls reproduce
    Cimarron's sequence closely (needs the clean-room data)."""
    pytest.importorskip("Bio")
    from Bio import Align
    path = CLEANROOM / "m13" / "esd" / "A02.esd"
    if not path.exists():
        pytest.skip("clean-room m13 data not available")
    b, tags = _esd(path)
    S, E = int(tags["START TRACE PROCESSING"]), int(tags["END TRACE PROCESSING"])
    proc = np.frombuffer(b[:(E - S) * 20], dtype=np.dtype([("cur", "<i4"), ("ch", "<f4", 4)]))["ch"].astype(float)
    res = cc.call_trace(proc, "TGCA", S, flattened=True)
    ref = tags["SEQUENCE"]
    al = Align.PairwiseAligner(mode="global", match_score=1, mismatch_score=-1, open_gap_score=-2,
                               extend_gap_score=-1, end_deletion_score=0, end_insertion_score=0)
    a = al.align(ref, res.sequence)[0]
    agree = sum(ref[t0 + k] == res.sequence[q0 + k]
                for (t0, t1), (q0, q1) in zip(*a.aligned) for k in range(t1 - t0))
    assert agree / len(ref) > 0.97
