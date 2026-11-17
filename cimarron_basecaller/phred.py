"""Phred base qualities for ``track_bases`` calls.

``TrackedBase.height`` is a normalised peak height, not an error probability:
written straight into a FASTQ it read as Q0-Q3 for almost every base.  This
module turns each call into a calibrated Phred score, Q = -10*log10(p_error),
with a small logistic model of p_error.

Calibration (2026-11-17): every called base of the 96-well MB1000_M13_DT plate
was labelled against M13mp18 (M77815.1; BLAST to place the read, then a local
alignment -- a mismatch, an inserted base or an unaligned base counts as an
error).  Fitted on the even-numbered columns, checked on the odd ones:

==========  ===========  ==========  ======================================
config      error rate   AUC         held-out calibration
==========  ===========  ==========  ======================================
golden      13.5 %       0.949       predicted vs observed Q within ~1.5 up
                                     to Q25; conservative above (Q26 -> 36)
precision   1.4 %        0.896       within ~1 Q up to Q25; conservative
                                     at Q25-30 (27 -> 31)
==========  ===========  ==========  ======================================

Height alone gives AUC 0.77 / 0.74; the neighbourhood, spacing and
position-in-read terms carry the rest.  The coefficients are specific to the
config the read was called with (the precision config's deep quality trim
removes most of the tail the golden model has to price), and to MegaBACE 1000
M13-like reads; re-fit them for a different instrument or read type.
"""
from __future__ import annotations

from typing import Iterable, List, Sequence

import numpy as np
from scipy.ndimage import minimum_filter1d, uniform_filter1d

# intercept, then one weight per feature in ``_features`` order
MODELS = {
    "golden": (-2.78059, -0.814353, -0.463298, 1.47614, 0.549588, -0.549768,
               2.65294, -17.9316, 22.9528),
    "precision": (-1.67444, -0.87409, -0.475065, -1.87207, 2.09169, -1.11694,
                  2.25432, -6.00894, 3.05199),
}
Q_MAX = 40          # the plate holds too few errors to vouch for more


def _features(bands: Sequence) -> np.ndarray:
    """log height, its 5-base minimum and 11-base mean, |position - predicted|,
    |gap to previous - 1| and |gap to next - 1| (in spacings), and relative
    position in the read with its square."""
    h = np.maximum(np.array([b.height for b in bands], float), 1e-3)
    lh = np.log(h)
    pos = np.array([b.position for b in bands], float)
    exp = np.array([b.expected_position for b in bands], float)
    sp = np.maximum(np.array([b.spacing_used for b in bands], float), 1e-6)
    n = len(bands)
    gap = np.diff(pos)
    gp = np.r_[1.0, gap / sp[1:]] if n > 1 else np.ones(1)
    gn = np.r_[gap / sp[:-1], 1.0] if n > 1 else np.ones(1)
    rel = np.arange(n) / max(n - 1, 1)
    return np.column_stack([
        lh,
        minimum_filter1d(lh, 5, mode="nearest"),
        uniform_filter1d(lh, 11, mode="nearest"),
        np.abs(pos - exp) / sp,
        np.abs(gp - 1.0),
        np.abs(gn - 1.0),
        rel,
        rel ** 2,
    ])


def error_probabilities(bands: Sequence, model: str = "golden") -> np.ndarray:
    """Estimated probability that each called base is wrong."""
    if not len(bands):
        return np.zeros(0)
    w = np.asarray(MODELS[model], float)
    z = w[0] + _features(bands) @ w[1:]
    return 1.0 / (1.0 + np.exp(-z))


def phred_qualities(bands: Iterable, model: str = "golden") -> List[int]:
    """Integer Phred qualities (0..Q_MAX) for a list of ``TrackedBase``."""
    bands = list(bands)
    p = np.clip(error_probabilities(bands, model), 10 ** (-Q_MAX / 10), 1.0)
    return [int(q) for q in np.clip(np.rint(-10 * np.log10(p)), 0, Q_MAX)]
