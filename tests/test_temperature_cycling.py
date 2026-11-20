"""CTCE runs are recognised from their electrophoresis current: a cycled
capillary temperature makes the current ripple; one temperature does not."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import analyzer_core as ac  # noqa: E402


def _decline(n=6000):
    x = np.arange(n, dtype=float)
    return 120.0 - 40.0 * (1 - np.exp(-x / 1500.0))


def test_a_cycled_run_ripples():
    n = 6000
    c = _decline(n)
    x = np.arange(n)
    on = x < 2400                                   # 20 cycles, then held
    c = c + np.where(on, 3.0 * np.sin(2 * np.pi * x / 120.0), 0.0)
    c = np.round(c * 10) / 10                       # the instrument's 0.1 steps
    verdict, score, per = ac.temperature_cycling(c)
    assert verdict == "cycling" and abs(per - 120) <= 6


def test_one_temperature_declines_in_steps():
    c = np.round(_decline() / 2) * 2               # coarse staircase
    verdict, score, _ = ac.temperature_cycling(c)
    assert verdict in ("constant", "check") and score < ac.CYCLING_YES


def test_no_current_is_unknown():
    assert ac.temperature_cycling(None)[0] == "unknown"
