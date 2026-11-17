"""
Clean-room clone of the Cimarron 3.12 base caller (MegaBACE).

Chain: raw trace -> ``cimarron_signal.process`` (span, baseline, colour
separation, homomorphic deconvolution, dye-mobility delays, flattening) ->
``cimarron_calling.call_trace`` (fuzzy band editing, letters, qualities,
stitching, good-sequence limits).  Both stages were written from functional
specifications and black-box measurements only; see CLEANROOM_rewrite.md (and
CLEANROOM.md for phase 1) for the provenance of every rule.

Entry points
------------
``call(trace, base_order)``
    The evaluator interface: *trace* is an (n, 4) raw trace whose column i
    holds base ``base_order[i]`` (e.g. ``to_acgt_trace`` output with "ACGT");
    returns ``(sequence, positions)``.  Letters use *base_order*'s case; N is
    returned as "N" / "n".  Positions follow Cimarron's PEAK POSITIONS
    convention (raw scan index + 2).
``call_rsd(path)``
    Raw .rsd file -> ``CloneResult`` (sequence, positions, qualities,
    good-sequence limits, the signal result), using the file's own channel
    order and machine id.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import cimarron_calling as calling
from . import cimarron_signal as signal

# Channel order the PLAYER1 profile (the default dye chemistry) is written in;
# a trace given by base letters is put into this order before processing.
DEFAULT_DYE_ORDER = "TGCA"

# Scans 200..1200 precede the injected sample's arrival in every capillary of
# both plates (phase 1), so they give the per-channel pre-run level and noise.
PRE_RUN = (200, 1200)


def prerun_z(trace: np.ndarray) -> np.ndarray:
    """Per-channel z-score of the raw trace against the pre-run level/noise."""
    pre = trace[PRE_RUN[0]:PRE_RUN[1]]
    lvl = np.median(pre, axis=0)
    sd = 1.4826 * np.median(np.abs(pre - lvl), axis=0) + 1e-9
    return (trace - lvl) / sd


def has_signal(trace: np.ndarray, threshold: float = 12.0) -> bool:
    """False for a capillary that never rises above its pre-run noise (phase 1:
    the 90th percentile of the max-channel z-score separates the six empty
    mb4000 wells, 0.3-7.7, from every well with signal, >= 20).  Cimarron emits
    ~100 noise calls on such wells, which no caller can reproduce."""
    trace = np.asarray(trace, float)
    if len(trace) < PRE_RUN[1] + 500:
        return True
    return float(np.percentile(prerun_z(trace).max(axis=1), 90)) >= threshold


@dataclass
class CloneResult:
    sequence: str
    positions: np.ndarray
    qualities: np.ndarray
    good_start: int
    good_end: int
    signal: signal.SignalResult | None
    dye_order: str


def call_instrument(raw: np.ndarray, dye_order: str, meta: dict | None = None) -> CloneResult:
    """Raw trace in instrument channel order -> CloneResult (letters in upper case)."""
    raw = np.asarray(raw, float)
    if not has_signal(raw):
        return CloneResult("", np.zeros(0, int), np.zeros(0), 0, 0, None, dye_order)
    sig = signal.process(raw, dye_order, meta or {})
    r = calling.call_trace(sig.caller_trace, dye_order, sig.start_tp)
    return CloneResult(r.sequence, r.positions, r.qualities, r.good_start, r.good_end, sig, dye_order)


def call_rsd(path: str) -> CloneResult:
    """Raw .rsd file -> CloneResult, with the file's own channel order and metadata."""
    from .rsd_io import read_rsd
    rsd = read_rsd(path)
    meta = signal.read_meta(rsd.footer_raw)
    order = signal.dye_order_from_footer(rsd.footer_raw) or DEFAULT_DYE_ORDER
    return call_instrument(rsd.trace, order, meta)


def call(trace: np.ndarray, base_order: str = "ACGT"):
    """Cimarron-3.12-like calls for a trace whose column i holds base_order[i]:
    returns (sequence, positions in Cimarron's coordinates)."""
    trace = np.asarray(trace, float)
    up = base_order.upper()
    if sorted(up) != list("ACGT"):
        raise ValueError(f"base_order must be a permutation of ACGT, got {base_order!r}")
    raw = np.column_stack([trace[:, up.index(b)] for b in DEFAULT_DYE_ORDER])
    res = call_instrument(raw, DEFAULT_DYE_ORDER)
    lower = base_order.islower()
    table = {b: base_order[up.index(b)] for b in "ACGT"}
    table["N"] = "n" if lower else "N"
    seq = "".join(table[c] for c in res.sequence)
    return seq, np.asarray(res.positions, dtype=int)
