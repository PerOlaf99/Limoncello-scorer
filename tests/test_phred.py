"""Calibrated FASTQ qualities from cimarron_basecaller.phred.

The calibration itself was measured on the 96-well M13 plate (see the module
docstring); these tests pin its contract and catch a regression back to the
old behaviour of writing truncated peak heights (Q0-Q3) as Phred.
"""
import pathlib
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import basecall  # noqa: E402
from cimarron_basecaller.phred import Q_MAX, error_probabilities, phred_qualities  # noqa: E402
from cimarron_basecaller.spacing_caller import TrackedBase  # noqa: E402

M13 = ROOT / "example_data" / "M13" / "A01.rsd"


def _bands(heights, jitter=0.0):
    out, pos = [], 1000.0
    for h in heights:
        pos += 12.0
        out.append(TrackedBase(position=int(pos + jitter), channel=0, height=float(h),
                               expected_position=pos, spacing_used=12.0))
    return out


def test_empty_read_has_no_qualities():
    assert phred_qualities([]) == []


def test_qualities_are_integers_in_range():
    q = phred_qualities(_bands(np.linspace(0.05, 3.0, 200)))
    assert len(q) == 200
    assert all(isinstance(v, int) and 0 <= v <= Q_MAX for v in q)


def test_a_base_weaker_than_its_neighbours_is_trusted_less():
    # Heights are locally normalised, so the model prices a base against its
    # neighbourhood rather than on an absolute scale.
    even = error_probabilities(_bands([1.5] * 40))[20]
    dip = [1.5] * 40
    dip[20] = 0.2
    assert error_probabilities(_bands(dip))[20] > even


def test_a_misplaced_peak_is_trusted_less():
    on_time = error_probabilities(_bands([1.5] * 40))[20]
    late = _bands([1.5] * 40)
    late[20] = TrackedBase(position=late[20].position + 5, channel=0, height=1.5,
                           expected_position=late[20].expected_position, spacing_used=12.0)
    assert error_probabilities(late)[20] > on_time


@pytest.mark.skipif(not M13.exists(), reason="example M13 trace not present")
@pytest.mark.parametrize("config", ["golden", "precision"])
def test_fastq_qualities_are_phred_not_peak_heights(config):
    cfg = basecall.WIN_CONFIG if config == "golden" else basecall.WIN_CONFIG_PRECISION
    seq, quals = basecall.basecall_well(str(M13), cfg)
    assert len(seq) == len(quals) > 300
    q = np.asarray(quals)
    # The old output truncated heights (~0-3) to ints: nearly every base Q0-Q3.
    assert np.median(q) >= 10
    assert q.max() <= Q_MAX
