"""
Dynamic-programming base caller.

Both simple_caller.py (independent per-channel peak picking) and
spacing_caller.py (greedy left-to-right spacing tracker) hit a similar
accuracy ceiling (~70-80%) with a persistent tension between coverage and
identity -- fundamentally because greedy left-to-right decisions can't
recover from one bad local call, and independent peak picking ignores the
periodicity prior entirely. This module fixes that properly: generate a
generous candidate-peak pool (favor not missing real peaks over
precision), then use dynamic programming to find the single path through
those candidates that jointly maximizes peak strength while penalizing
deviation from a smoothly-varying expected local spacing -- so the whole
called sequence is chosen together, not one greedy step at a time.
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.signal import find_peaks

from .simple_caller import robust_baseline_subtract, normalize_channels_local, apply_mobility_correction
from .spacing_caller import find_start_position, _windowed_area


@dataclass
class Candidate:
    position: int
    channel: int
    area: float


def generate_candidates(norm_trace: np.ndarray, min_prominence: float, half_width: float) -> list[Candidate]:
    """Generous per-channel candidate peaks using local-area strength
    (not just point height), erring toward over-generating candidates --
    the DP step is responsible for picking the right subset, so precision
    isn't needed here, only recall."""
    candidates = []
    n = norm_trace.shape[0]
    for ch in range(norm_trace.shape[1]):
        idx, _ = find_peaks(norm_trace[:, ch], prominence=min_prominence, distance=2)
        if len(idx) == 0:
            continue
        areas = _windowed_area(norm_trace[:, ch], 0, n - 1, half_width)
        for i in idx:
            candidates.append(Candidate(position=int(i), channel=ch, area=float(areas[i])))
    candidates.sort(key=lambda c: c.position)
    return candidates


def fit_spacing_curve(positions: np.ndarray, spacings: np.ndarray, trace_len: int, degree: int = 2) -> np.ndarray:
    """Smooth (polynomial) fit of local spacing vs. position, evaluated at
    every position 0..trace_len-1, used as the DP's expected-spacing
    reference (spacing genuinely drifts over a read -- broader peaks and
    wider spacing later in the run, per migration-time diffusion)."""
    if len(positions) < degree + 1:
        return np.full(trace_len, float(np.median(spacings)) if len(spacings) else 12.0)
    coeffs = np.polyfit(positions, spacings, deg=degree)
    xs = np.arange(trace_len)
    fitted = np.polyval(coeffs, xs)
    return np.clip(fitted, 2.0, None)


def dp_call_bases(
    trace: np.ndarray,
    base_order: str = "ACGT",
    min_prominence: float = 0.05,
    area_half_width_frac: float = 0.35,
    initial_spacing: float = 11.5,
    spacing_penalty_weight: float = 0.15,
    max_lookback_factor: float = 2.3,
    min_lookback_factor: float = 0.35,
    baseline_window: int = 151,
    local_norm_window: int = 300,
    mobility_shifts: list[int] | None = None,
) -> tuple[str, list[float]]:
    baseline_subtracted = robust_baseline_subtract(trace, window=baseline_window)

    from .spacing_caller import detect_signal_region
    sig_start, sig_end = detect_signal_region(baseline_subtracted)

    norm_trace_full = normalize_channels_local(baseline_subtracted, window=local_norm_window)
    if mobility_shifts is not None:
        norm_trace_full = apply_mobility_correction(norm_trace_full, mobility_shifts)

    # Restrict everything downstream to the real signal region -- local
    # normalization blows pre-injection/post-run noise up to full scale,
    # so start/end must be fixed using the global detector BEFORE any
    # candidate generation happens in that region.
    norm_trace = norm_trace_full[sig_start:sig_end]

    n = norm_trace.shape[0]
    half_width = max(1.0, initial_spacing * area_half_width_frac)
    candidates = generate_candidates(norm_trace, min_prominence, half_width)
    if not candidates:
        return "", []

    # First pass: cheap greedy spacing tracker just to get an empirical
    # spacing-vs-position profile to fit a smooth curve against.
    from .spacing_caller import track_bases
    _, _, tracked = track_bases(
        trace[sig_start:sig_end], base_order=base_order, initial_spacing=initial_spacing,
        mobility_shifts=mobility_shifts, baseline_window=baseline_window,
        local_norm_window=local_norm_window,
    )
    if len(tracked) >= 5:
        positions = np.array([b.position for b in tracked[1:]])
        spacings = np.array([b.spacing_used for b in tracked[1:]])
        spacing_curve = fit_spacing_curve(positions, spacings, n)
    else:
        spacing_curve = np.full(n, initial_spacing)

    areas = np.array([c.area for c in candidates])
    area_scale = np.percentile(areas, 90) if len(areas) else 1.0
    norm_areas = areas / (area_scale + 1e-9)

    positions_arr = np.array([c.position for c in candidates])
    m = len(candidates)

    dp = np.full(m, -np.inf)
    backptr = np.full(m, -1, dtype=int)

    # Candidates are position-sorted; for each i, only look back at
    # candidates within a plausible spacing window of i (bounded lookback
    # keeps this roughly linear rather than O(m^2)).
    start_idx = 0  # detect_signal_region already trimmed to the real signal region

    for i in range(m):
        pos_i = positions_arr[i]
        expected = spacing_curve[pos_i]
        lo_gap = expected * min_lookback_factor
        hi_gap = expected * max_lookback_factor

        # Binary-search-ish bounded window over j via direct scan back
        # (candidate density is low enough this stays cheap).
        j = i - 1
        best_score, best_j = -np.inf, -1
        while j >= 0 and (pos_i - positions_arr[j]) <= hi_gap:
            gap = pos_i - positions_arr[j]
            if gap >= lo_gap:
                local_expected = spacing_curve[positions_arr[j]]
                penalty = spacing_penalty_weight * abs(gap - local_expected) / max(local_expected, 1.0)
                score = dp[j] - penalty
                if score > best_score:
                    best_score, best_j = score, j
            j -= 1

        base_score = norm_areas[i]
        if best_j == -1:
            # No valid predecessor: can still start here (e.g. the very
            # first real base), with a small bonus if it's near the
            # detected start-of-signal region.
            start_bonus = 0.0 if pos_i - start_idx < expected * 3 else -1.0
            dp[i] = base_score + start_bonus
        else:
            dp[i] = best_score + base_score
        backptr[i] = best_j

    # Traceback from the highest-scoring endpoint NEAR THE END of the
    # trace (not the global argmax over all candidates) -- otherwise the
    # DP can "win" by picking a short, locally strong subsequence instead
    # of calling the whole read, since summed reward doesn't inherently
    # favor longer paths.
    last_pos = positions_arr[-1]
    near_end = np.where(positions_arr >= last_pos - initial_spacing * 5)[0]
    if len(near_end) == 0:
        near_end = np.array([m - 1])
    end_i = int(near_end[np.argmax(dp[near_end])])
    path = []
    cur = end_i
    while cur != -1:
        path.append(cur)
        cur = backptr[cur]
    path.reverse()

    chosen = [candidates[i] for i in path]
    sequence = "".join(base_order[c.channel] for c in chosen)
    qualities = [norm_areas[path[k]] for k in range(len(path))]
    return sequence, qualities
