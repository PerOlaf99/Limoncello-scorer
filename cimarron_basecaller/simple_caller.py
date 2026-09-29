"""
A simpler, empirically-tuned base-calling pipeline: per-channel peak
detection with cross-channel amplitude normalization, rather than the
patent's under-specified blind-deconvolution + fuzzy-logic pipeline
(see pipeline.py / peaks.py / quality.py for that attempt).

Rationale: the patent (EP0944739A1) leaves several core pieces
under-specified (the exact deconvolution kernel, exact OmitOkN/GapCheck
membership shapes), and a faithful-as-possible reconstruction of it
scored very poorly against real ground truth (~2% query coverage). This
module instead uses standard, well-understood Sanger-trace base-calling
techniques -- per-channel local-maxima peak detection, robust baseline
subtraction, and cross-channel normalization -- tuned and validated
directly against real MegaBACE data (A01, matched to its actual Cimarron
3.12 SCF output) rather than against an incomplete textual spec.
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import minimum_filter1d, uniform_filter1d


@dataclass
class CalledBase:
    position: int
    channel: int
    height: float
    quality: float


def robust_baseline_subtract(trace: np.ndarray, window: int = 151) -> np.ndarray:
    """Rolling-minimum baseline subtraction per channel: a simple, robust
    alternative to the patent's rising-exponential-threshold algorithm.
    A rolling minimum tracks the noise floor; smoothing it avoids chasing
    individual low-noise samples."""
    out = np.empty_like(trace)
    for c in range(trace.shape[1]):
        col = trace[:, c]
        baseline = minimum_filter1d(col, size=window, mode="nearest")
        baseline = uniform_filter1d(baseline, size=window // 3 + 1, mode="nearest")
        out[:, c] = np.clip(col - baseline, 0, None)
    return out


def normalize_channels(trace: np.ndarray, percentile: float = 99.0) -> np.ndarray:
    """Scale each channel so its high-percentile amplitude is 1.0, so no
    single channel's raw dynamic range dominates cross-channel comparison
    (channels were observed to have very different raw amplitude ranges)."""
    out = np.empty_like(trace)
    for c in range(trace.shape[1]):
        scale = np.percentile(trace[:, c], percentile)
        out[:, c] = trace[:, c] / scale if scale > 1e-9 else trace[:, c]
    return out


def normalize_channels_local(trace: np.ndarray, window: int = 400, percentile: float = 95.0) -> np.ndarray:
    """Like normalize_channels, but the scale factor is a ROLLING local
    high-percentile rather than one global value -- compensates for signal
    amplitude naturally rising/falling over the length of a read (e.g.
    weaker signal near the very start/end), which a single global
    normalization factor doesn't account for.
    """
    from scipy.ndimage import maximum_filter1d, percentile_filter

    out = np.empty_like(trace)
    for c in range(trace.shape[1]):
        local_scale = percentile_filter(trace[:, c], percentile=percentile, size=window, mode="nearest")
        local_scale = np.maximum(local_scale, np.percentile(trace[:, c], 90) * 0.1 + 1e-6)
        out[:, c] = trace[:, c] / local_scale
    return out


def detect_candidate_peaks(norm_trace: np.ndarray, min_distance: int, prominence: float) -> list[CalledBase]:
    """Find local-maxima peaks independently on each (normalized) channel."""
    candidates = []
    for c in range(norm_trace.shape[1]):
        idx, props = find_peaks(
            norm_trace[:, c], distance=min_distance, prominence=prominence
        )
        for i, p in zip(idx, props["prominences"]):
            candidates.append(CalledBase(position=int(i), channel=c, height=float(norm_trace[i, c]), quality=float(p)))
    candidates.sort(key=lambda b: b.position)
    return candidates


def resolve_conflicts(candidates: list[CalledBase], min_distance: int) -> list[CalledBase]:
    """Where candidate peaks from different channels fall within
    `min_distance` of each other, keep only the tallest (winner-take-all) --
    this is how a single base identity gets assigned at each trace
    position despite detecting peaks per-channel."""
    if not candidates:
        return []
    resolved = [candidates[0]]
    for cand in candidates[1:]:
        if cand.position - resolved[-1].position < min_distance:
            if cand.height > resolved[-1].height:
                resolved[-1] = cand
        else:
            resolved.append(cand)
    return resolved


def filter_by_dominance(
    resolved: list[CalledBase], norm_trace: np.ndarray, min_ratio: float = 1.3
) -> list[CalledBase]:
    """Drop any called base whose winning channel doesn't clearly dominate
    the runner-up channel at that trace position -- an ambiguous/noisy
    call where two channels are nearly tied is more likely a spurious
    peak (cross-talk, baseline noise) than a real base."""
    out = []
    for b in resolved:
        row = norm_trace[b.position]
        sorted_vals = np.sort(row)[::-1]
        top, second = sorted_vals[0], sorted_vals[1]
        ratio = top / second if second > 1e-9 else np.inf
        if ratio >= min_ratio:
            out.append(b)
    return out


def estimate_spacing(trace: np.ndarray, norm_trace: np.ndarray) -> float:
    """Rough initial base-spacing estimate via autocorrelation of the
    combined channel envelope, used to set peak-detection distance."""
    envelope = norm_trace.max(axis=1)
    envelope = envelope - envelope.mean()
    n = len(envelope)
    fft = np.fft.rfft(envelope, n=2 * n)
    acf = np.fft.irfft(fft * np.conj(fft))[:n]
    # search for the first strong peak beyond a minimum plausible spacing
    search_lo, search_hi = 4, 40
    local = acf[search_lo:search_hi]
    if len(local) == 0 or local.max() <= 0:
        return 12.0
    best = search_lo + int(np.argmax(local))
    return float(best)


def apply_mobility_correction(trace: np.ndarray, shifts: list[int] | tuple[int, ...]) -> np.ndarray:
    """Shift each channel's samples by a per-channel integer offset to
    correct for dye mobility differences (each dye migrates through the
    capillary at a slightly different rate, so the same physical base
    shows up at a slightly different scan index per channel). `shifts[c]`
    is the amount channel c's peaks are observed to be delayed by --
    positive means the channel's real peak appears `shifts[c]` samples
    LATER than it should, so we shift it earlier to correct.
    """
    out = np.zeros_like(trace)
    n = trace.shape[0]
    for c, s in enumerate(shifts):
        s = int(round(s))
        if s == 0:
            out[:, c] = trace[:, c]
        elif s > 0:
            out[:n - s, c] = trace[s:, c]
        else:
            out[-s:, c] = trace[:n + s, c]
    return out


def call_bases(
    trace: np.ndarray,
    base_order: str = "ACGT",
    min_distance: int | None = None,
    prominence: float = 0.03,
    baseline_window: int = 151,
    mobility_shifts: list[int] | None = None,
    use_local_norm: bool = False,
    local_norm_window: int = 400,
    dominance_ratio: float | None = None,
) -> tuple[str, list[float], list[CalledBase]]:
    """Full simplified pipeline: baseline subtract -> normalize channels ->
    (optional mobility correction) -> detect per-channel peaks -> resolve
    cross-channel conflicts -> (optional) dominance filter -> call bases.
    """
    baseline_subtracted = robust_baseline_subtract(trace, window=baseline_window)
    if use_local_norm:
        norm_trace = normalize_channels_local(baseline_subtracted, window=local_norm_window)
    else:
        norm_trace = normalize_channels(baseline_subtracted)
    if mobility_shifts is not None:
        norm_trace = apply_mobility_correction(norm_trace, mobility_shifts)

    if min_distance is None:
        spacing = estimate_spacing(trace, norm_trace)
        min_distance = max(3, int(round(spacing * 0.6)))

    candidates = detect_candidate_peaks(norm_trace, min_distance=min_distance, prominence=prominence)
    resolved = resolve_conflicts(candidates, min_distance=min_distance)
    if dominance_ratio is not None:
        resolved = filter_by_dominance(resolved, norm_trace, min_ratio=dominance_ratio)

    sequence = "".join(base_order[b.channel] for b in resolved)
    qualities = [b.quality for b in resolved]
    return sequence, qualities, resolved
