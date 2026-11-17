"""Interface contract of the clean-room Cimarron 3.12 clone (cimarron_clone.call)."""
from pathlib import Path

import numpy as np
import pytest

from cimarron_basecaller.rsd_io import read_rsd, to_acgt_trace

M13 = Path(__file__).resolve().parent.parent / "example_data" / "M13"


def _trace(name="A01.rsd"):
    path = M13 / name
    if not path.exists():
        pytest.skip(f"example_data/M13/{name} missing")
    return to_acgt_trace(read_rsd(str(path)), base_order="TGCA")


def test_call_returns_sequence_and_int_positions():
    from cimarron_basecaller.cimarron_clone import call
    trace, order = _trace()
    seq, pos = call(trace, order)
    assert isinstance(seq, str) and len(seq) > 300
    assert set(seq) <= set(order)
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


def test_shared_signal_path_matches_track_bases_golden():
    """The clone re-implements the GOLDEN signal path; pin it to track_bases."""
    import configs
    from cimarron_basecaller import track_bases
    from cimarron_basecaller.cimarron_clone import golden_equivalent
    trace, order = _trace()
    seq, _q, bands = track_bases(trace, base_order=order, **configs.GOLDEN)
    seq2, bands2 = golden_equivalent(trace, order)
    assert seq == seq2
    assert [b.position for b in bands] == [b.position for b in bands2]
