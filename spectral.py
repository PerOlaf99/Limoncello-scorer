"""Dye overlap (spectral) correction for multi-dye fragment runs.

Every dye shows up a little -- or a lot -- in its neighbours' channels: on the
MegaBACE genotyping plates HEX reads at 56-109 % in the NED channel and the
ET-ROX ladder at 8-14 % in FAM and NED.  Uncorrected, each strong peak leaves
"pull-up" copies in the other channels, and those copies get sized as
fragments.  Genetic Profiler removes them with a spectral overlap matrix; so
does this module.

The matrix ``M`` is ``channels x dyes``: ``M[c, d]`` is how much of dye ``d``
the detector sees in channel ``c`` (the dye's own channel is 1).  A measured
scan is ``observed = M @ true``, so the correction is
``true = inverse(M) @ observed``.

Where a matrix comes from:

* a Genetic Profiler ``xtalk.tdb`` (named 4x4 matrices -- "auto" is the
  identity, a run's own calibration sits in ``<run>/xtm/xtalk.tdb``), or
* single-dye matrix-standard wells of the plate (each well holds one dye):
  :func:`matrix_from_wells` measures how each dye's peaks read in every
  channel, as GP's "Create New" does.

Headless (numpy only); the CLI, the GUI and the sizing table share it.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np


@dataclass
class SpectralMatrix:
    name: str
    m: np.ndarray                     # channels x dyes, own-channel 1
    source: str = ""
    dyes: List[str] = field(default_factory=list)

    @property
    def n(self) -> int:
        return int(self.m.shape[0])

    def correct(self, channels: np.ndarray) -> np.ndarray:
        """Per-dye signal from (scans, n) measured channels."""
        x = np.asarray(channels, float)
        k = self.n
        if x.ndim != 2 or x.shape[1] < k:
            raise ValueError(f"need {k} channels, got {x.shape}")
        out = x.copy()
        out[:, :k] = x[:, :k] @ np.linalg.inv(self.m).T
        return out

    def is_identity(self) -> bool:
        return bool(np.allclose(self.m, np.eye(self.n)))


# --------------------------------------------------------------------------- #
# Genetic Profiler xtalk.tdb
# --------------------------------------------------------------------------- #
def read_xtalk_tdb(path) -> Dict[str, SpectralMatrix]:
    """The named matrices of a Genetic Profiler ``xtalk.tdb``: a name line,
    then one line per channel of whitespace-separated values."""
    path = Path(path)
    lines = [ln.strip() for ln in path.read_text(encoding="latin-1").splitlines()]
    out: Dict[str, SpectralMatrix] = {}
    i = 0
    while i < len(lines):
        name = lines[i]
        i += 1
        if not name:
            continue
        rows = []
        while i < len(lines):
            try:
                vals = [float(t) for t in lines[i].split()]
            except ValueError:
                break
            if not vals:
                break
            rows.append(vals)
            i += 1
        if rows and all(len(r) == len(rows) for r in rows):
            out[name] = SpectralMatrix(name, np.asarray(rows, float),
                                       source=f"{path.name}: {name}")
    return out


def run_matrices(run_folder) -> Dict[str, SpectralMatrix]:
    """Matrices stored with a run (``<run>/xtm/xtalk.tdb``), if any."""
    p = Path(run_folder) / "xtm" / "xtalk.tdb"
    return read_xtalk_tdb(p) if p.is_file() else {}


# --------------------------------------------------------------------------- #
# a matrix from matrix-standard wells
# --------------------------------------------------------------------------- #
def _dye_ratios(x: np.ndarray, own: int, top: int = 40) -> Optional[np.ndarray]:
    """How a single-dye trace reads in every channel, relative to its own:
    the median ratio over its tallest peaks (baseline removed), so one
    saturated or odd peak does not decide it."""
    from scipy.signal import find_peaks
    x = np.asarray(x, float)
    x = x - np.median(x, axis=0)
    y = x[:, own]
    idx, props = find_peaks(y, prominence=0.1 * max(y.max(), 1e-9), distance=5)
    if len(idx) == 0:
        return None
    idx = idx[np.argsort(props["prominences"])[::-1][:top]]
    sat = y.max() * 0.98                     # skip clipped tops
    idx = [i for i in idx if y[i] < sat] or list(idx)
    ratios = x[idx, :] / y[idx, None]
    col = np.median(ratios, axis=0)
    col[own] = 1.0
    return np.clip(col, 0.0, None)


def matrix_from_wells(traces_by_dye: Dict[int, Sequence[np.ndarray]],
                      name: str = "from plate", dyes: Sequence[str] = ()) -> SpectralMatrix:
    """Build ``M`` from single-dye wells.

    *traces_by_dye* maps a dye's own channel index (0-based) to the measured
    (scans, n) traces of the wells that hold only that dye; each dye column
    is the median over its wells."""
    n = max(traces_by_dye) + 1
    for tr in (t for ts in traces_by_dye.values() for t in ts):
        n = max(n, np.shape(tr)[1])
    m = np.eye(n)
    missing = []
    for d in range(n):
        cols = [c for c in (_dye_ratios(t, d) for t in traces_by_dye.get(d, []))
                if c is not None]
        if cols:
            m[:, d] = np.median(np.asarray(cols)[:, :n], axis=0)
        else:
            missing.append(d + 1)
    src = "matrix-standard wells" + (f" (no wells for channel(s) {missing})"
                                     if missing else "")
    return SpectralMatrix(name, m, source=src, dyes=list(dyes))


def identity(n: int = 4) -> SpectralMatrix:
    return SpectralMatrix("none", np.eye(n), source="no correction")


# --------------------------------------------------------------------------- #
# finding the matrix-standard wells on a plate
# --------------------------------------------------------------------------- #
def dye_profile(channels: np.ndarray, max_leftover: float = 0.1):
    """(profile, leftover) of a well: how every channel reads relative to its
    strongest one, and how much signal a scaled copy of that one leaves
    unexplained in any channel (relative to the strongest's amplitude).  A
    single-dye well leaves almost nothing; a sample well with independent
    peaks in several dyes leaves a lot."""
    x = np.asarray(channels, float)
    x = x - np.median(x, axis=0)
    amp = np.percentile(x, 99.9, axis=0)
    d = int(np.argmax(amp))
    y = x[:, d]
    yy = float(np.dot(y, y))
    if yy <= 0 or amp[d] <= 0:
        return None, 1.0
    a = x.T @ y / yy
    res = x - np.outer(y, a)
    leftover = float(np.max(np.percentile(np.abs(res), 99.9, axis=0)) / amp[d])
    return np.clip(a, 0, None) / max(a.max(), 1e-12), leftover


def find_matrix_wells(channels_by_well: Dict[str, np.ndarray],
                      max_leftover: float = 0.1) -> Dict[int, List[str]]:
    """The plate's single-dye wells, grouped by dye: {channel index: wells}.

    Wells are grouped by their channel profile, then every group gets its
    own channel by an optimal one-to-one assignment -- a dye can read
    stronger in a neighbour's channel than in its own (HEX in the NED channel
    on MegaBACE plates), so the strongest channel alone cannot name it."""
    from scipy.optimize import linear_sum_assignment
    prof = {}
    for well, ch in channels_by_well.items():
        p, left = dye_profile(ch)
        if p is not None and left <= max_leftover:
            prof[well] = p
    groups: List[List[str]] = []
    for well in sorted(prof):
        for g in groups:
            if np.abs(prof[g[0]] - prof[well]).sum() < 0.25:
                g.append(well)
                break
        else:
            groups.append([well])
    if not groups:
        return {}
    n = len(next(iter(prof.values())))
    score = np.array([[np.mean([prof[w][c] for w in g]) for c in range(n)]
                      for g in groups])
    rows, cols = linear_sum_assignment(-score)
    out: Dict[int, List[str]] = {}
    for r, c in zip(rows, cols):
        if score[r, c] >= 0.5:            # the group really reads in that channel
            out[int(c)] = sorted(groups[r])
    return out


def plate_matrix(run_folder, channels_by_well: Optional[Dict[str, np.ndarray]] = None,
                 prefer: str = "auto") -> SpectralMatrix:
    """The matrix to use for a run.

    *prefer*: "run" -- the run's own Genetic Profiler calibration
    (``xtm/xtalk.tdb``, any matrix but the identity); "plate" -- measured
    from the plate's single-dye wells; "none"; "auto" (default) -- the run's
    calibration if it has one, else the plate's wells, else none."""
    if prefer in ("auto", "run"):
        for name, m in run_matrices(run_folder).items():
            if not m.is_identity():
                return m
        if prefer == "run":
            return identity()
    if prefer in ("auto", "plate") and channels_by_well:
        wells = find_matrix_wells(channels_by_well)
        if len(wells) >= 2:
            m = matrix_from_wells({d: [channels_by_well[w] for w in ws]
                                   for d, ws in wells.items()},
                                  name=f"{Path(run_folder).name} (plate wells)")
            m.source = "matrix-standard wells " + "; ".join(
                f"Ch{d + 1}: {', '.join(ws)}" for d, ws in sorted(wells.items()))
            return m
    return identity()
