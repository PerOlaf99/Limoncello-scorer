"""Interface contract of the clean-room Cimarron 3.12 clone (cimarron_clone)."""
from pathlib import Path

import numpy as np
import pytest

from cimarron_basecaller.rsd_io import read_rsd, to_acgt_trace

M13 = Path(__file__).resolve().parent.parent / "example_data" / "M13"


def _path(name="A01.rsd"):
    path = M13 / name
    if not path.exists():
        pytest.skip(f"example_data/M13/{name} missing")
    return path


def _trace(name="A01.rsd"):
    return to_acgt_trace(read_rsd(str(_path(name))), base_order="TGCA")


def test_call_returns_sequence_and_int_positions():
    from cimarron_basecaller.cimarron_clone import call
    trace, order = _trace()
    seq, pos = call(trace, order)
    assert isinstance(seq, str) and len(seq) > 300
    # Cimarron calls N where two channels compete; everything else is a base
    assert set(seq) <= set(order) | {"N"}
    assert seq.count("N") < 0.05 * len(seq)
    pos = np.asarray(pos)
    assert pos.shape == (len(seq),)
    assert np.issubdtype(pos.dtype, np.integer)
    assert np.all(np.diff(pos) > 0)


def test_call_respects_base_order_alphabet():
    from cimarron_basecaller.cimarron_clone import call
    trace, _ = _trace()
    seq_acgt, pos_acgt = call(trace, "ACGT")
    seq_lower, pos_lower = call(trace, "acgt")
    assert seq_lower == seq_acgt.lower()
    assert np.array_equal(pos_acgt, pos_lower)


def test_call_is_independent_of_column_order():
    """The same raw data given in another column order (with the matching
    base_order) gives the same calls."""
    from cimarron_basecaller.cimarron_clone import call
    trace, order = _trace()
    perm = [3, 1, 0, 2]
    seq1, pos1 = call(trace, order)
    seq2, pos2 = call(trace[:, perm], "".join(order[i] for i in perm))
    assert seq1 == seq2 and np.array_equal(pos1, pos2)


def test_call_is_deterministic():
    from cimarron_basecaller.cimarron_clone import call
    trace, order = _trace()
    a = call(trace, order)
    b = call(trace, order)
    assert a[0] == b[0] and np.array_equal(a[1], b[1])


def test_call_on_flat_trace_returns_no_bases():
    """A capillary without signal: Cimarron's calls there are noise, so the clone
    returns nothing rather than thousands of noise bases."""
    from cimarron_basecaller.cimarron_clone import call
    rng = np.random.default_rng(0)
    trace = 200.0 + rng.normal(0, 5, size=(9000, 4))
    seq, pos = call(trace, "ACGT")
    assert seq == "" and len(pos) == 0


def test_call_rsd_matches_call_and_reports_everything():
    from cimarron_basecaller.cimarron_clone import call, call_rsd
    res = call_rsd(str(_path()))
    trace, order = _trace()
    seq, pos = call(trace, order)
    assert res.dye_order == "TGCA"
    assert res.sequence == seq and np.array_equal(res.positions, pos)
    assert len(res.qualities) == len(seq)
    assert 0 <= res.good_start <= res.good_end < len(seq)
    sig = res.signal
    assert sig is not None and sig.estimated is False
    assert sig.trace.shape == sig.internal.shape == (sig.end_tp - sig.start_tp, 4)
    # positions are in Cimarron's PEAK POSITIONS convention: raw scan + 2
    assert sig.start_tp + 1 <= pos[0] and pos[-1] <= sig.end_tp + 1
