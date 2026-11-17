"""
Black-box clone of the Cimarron 3.12 base caller (MegaBACE), phase 1.

Built clean-room: every rule here was learned from Cimarron 3.12's own
outputs (.esd result files) on raw traces, or from the public patent
EP0944739A1.  See CLEANROOM.md at the repository root for the provenance of
each decision.

Entry point: ``call(trace, base_order) -> (sequence, positions)``.

The signal path is the same as ``track_bases`` with ``configs.GOLDEN``
(re-implemented here, lean, so the read region can be chosen freely without
touching track_bases); what differs is *where* calling starts and stops and the
coordinate convention of the returned positions.
"""

from __future__ import annotations

import numpy as np

from .simple_caller import robust_baseline_subtract, normalize_channels_local, apply_mobility_correction
from .spacing_caller import (
    DEFAULT_CHM, DEFAULT_MOBILITY_SHIFTS, TrackedBase, apply_gaussian_reconstruction_windowed,
    apply_spectral_separation, detect_signal_region, estimate_global_spacing,
    precompute_channel_peak_masks, smooth_trace, trim_to_quality_block, _ramp_value,
)

# GOLDEN preset values (configs.GOLDEN), frozen here so later preset tuning
# cannot silently change the clone.
_G = dict(
    smoothing_window=2, baseline_window=151, local_norm_window=1800,
    gaussian_recon_segment_size=384, gaussian_recon_noise_reg=0.06, gaussian_recon_sigma_scale=1.05,
    channel_peak_bonus=1.2, window_frac=(0.75, 1.25), min_prominence=0.05, ema_alpha=0.08,
    pullback_weight=(0.008, 0.001), profile_fracs=(0.33, 1.0), max_misses=6,
    trim_quality_percentile=10.0,
)


def _process(trace: np.ndarray):
    """GOLDEN signal path. Returns (baseline_subtracted, envelope, best_channel, norm_trace)."""
    g = _G
    bs = robust_baseline_subtract(smooth_trace(trace, g['smoothing_window']), window=g['baseline_window'])
    bs = apply_spectral_separation(bs, DEFAULT_CHM)
    bs = apply_mobility_correction(bs, DEFAULT_MOBILITY_SHIFTS)
    norm = normalize_channels_local(bs, window=g['local_norm_window'])
    s0, s1 = detect_signal_region(bs)
    sp = estimate_global_spacing(norm.max(axis=1), s0, s1)
    norm = apply_gaussian_reconstruction_windowed(
        norm, np.full(norm.shape[0], sp), segment_size=g['gaussian_recon_segment_size'],
        noise_reg=g['gaussian_recon_noise_reg'], sigma_scale=g['gaussian_recon_sigma_scale'])
    norm = np.clip(norm, 0, None)
    mask = precompute_channel_peak_masks(norm)
    boosted = norm * (1.0 + mask.astype(float) * g['channel_peak_bonus'])
    return bs, boosted.max(axis=1), boosted.argmax(axis=1), norm


def _track(envelope, best_channel, sig_start, sig_end, global_spacing, prof_start=None, prof_end=None,
           clamp=(0.5, 2.0)):
    """The track_bases forward walk (GOLDEN parameters), starting at sig_start.

    prof_start/prof_end define the scan-fraction axis of the parameter
    profiles (track_bases uses its detected signal region).  `clamp` bounds the
    running spacing relative to `global_spacing` (GOLDEN: 0.5-2.0)."""
    g = _G
    n = len(envelope)
    ps = sig_start if prof_start is None else prof_start
    pe = sig_end if prof_end is None else prof_end
    span = max(1, pe - ps)
    f0, f1 = g['profile_fracs']
    pos = sig_start
    spacing = global_spacing
    misses = 0
    out: list[TrackedBase] = []
    while pos < n:
        frac = (pos - ps) / span
        mp = g['min_prominence']
        ea = g['ema_alpha']
        pw = _ramp_value(g['pullback_weight'], frac, f0, f1)
        expected_next = pos + spacing
        lo = int(expected_next - spacing * (1 - g['window_frac'][0]))
        hi = int(expected_next + spacing * (g['window_frac'][1] - 1))
        lo = max(lo, pos + 2)
        hi = min(hi, n - 1)
        if lo > hi:
            break
        window = envelope[lo:hi + 1]
        if len(window) == 0:
            break
        k = int(np.argmax(window))
        pk, val = lo + k, window[k]
        if val < mp:
            misses += 1
            pos = int(expected_next)
            if misses > g['max_misses']:
                break
            continue
        misses = 0
        ref = global_spacing
        spacing = (1 - ea) * spacing + ea * (pk - pos)
        spacing = (1 - pw) * spacing + pw * ref
        spacing = float(np.clip(spacing, ref * clamp[0], ref * clamp[1]))
        out.append(TrackedBase(position=pk, channel=int(best_channel[pk]), height=float(val),
                               expected_position=expected_next, spacing_used=spacing))
        pos = pk
    return out


def golden_equivalent(trace: np.ndarray, base_order: str = "ACGT"):
    """Re-implementation of track_bases(**configs.GOLDEN); used by the tests to
    pin the shared signal path to the original."""
    bs, env, ch, _ = _process(trace)
    s0, s1 = detect_signal_region(bs)
    sp = estimate_global_spacing(env, s0, s1)
    bands = _track(env, ch, s0, s1, sp)
    if bands:
        a, b = trim_to_quality_block(bands, min_quality_percentile=_G['trim_quality_percentile'])
        bands = bands[a:b]
    return "".join(base_order[b.channel] for b in bands), bands


# ---------------------------------------------------------------------------
# Read span (where Cimarron starts and stops calling)
# ---------------------------------------------------------------------------

# Scans 200..1200 precede the injected sample's arrival in every capillary of
# both plates (signal never starts before scan ~1500), so they give the
# per-channel pre-run level and noise.
PRE_RUN = (200, 1200)


def prerun_z(trace: np.ndarray) -> np.ndarray:
    """Per-channel z-score of the raw trace against the pre-run level/noise."""
    pre = trace[PRE_RUN[0]:PRE_RUN[1]]
    lvl = np.median(pre, axis=0)
    sd = 1.4826 * np.median(np.abs(pre - lvl), axis=0) + 1e-9
    return (trace - lvl) / sd


def has_signal(z: np.ndarray, threshold: float = 12.0) -> bool:
    """False for a capillary that never rises above its pre-run noise.

    On such wells Cimarron still emits ~100 noise calls (in a coordinate frame
    past the end of the trace), which no caller can reproduce; returning nothing
    avoids thousands of spurious calls.  The 90th percentile of the max-channel
    z-score separates the six empty mb4000 wells (0.3-7.7) from every well with
    signal (>= 20)."""
    if len(z) < PRE_RUN[1] + 500:
        return True
    return float(np.percentile(z.max(axis=1), 90)) >= threshold


def signal_onset(z: np.ndarray, level: float = 8.0, sustain: int = 20, smooth: int = 11) -> int | None:
    """First scan after the pre-run window where the (smoothed) max-channel
    z-score stays above `level` for `sustain` scans: the arrival of the first
    labelled fragments, where Cimarron begins calling (it calls the faint
    early bands too, ~120 scans before the 10%-of-max region track_bases uses)."""
    from scipy.ndimage import uniform_filter1d
    s = uniform_filter1d(z.max(axis=1), smooth)
    above = s > level
    run = 0
    for i in range(PRE_RUN[1], len(s)):
        run = run + 1 if above[i] else 0
        if run >= sustain:
            return i - sustain + 1
    return None


def end_blob(trace: np.ndarray, search: int = 1500, ratio: float = 3.0):
    """Peak scan of the end-of-run blob (a large, broad rise in all channels
    followed by a fall to the pre-run floor, seen on most m13 capillaries), or
    None.  Cimarron's END RAW DATA sits next to it."""
    from scipy.ndimage import uniform_filter1d
    n = len(trace)
    if n < 2 * search:
        return None
    s = uniform_filter1d(trace.sum(axis=1), 51)
    seg = s[n - search:]
    k = int(np.argmax(seg))
    pk = n - search + k
    base = np.median(s[n - 2 * search:n - search])
    if s[pk] < ratio * base:
        return None
    # must come back down after the peak (not just a rising trace end)
    if s[pk:].min() > 0.5 * s[pk] and pk < n - 60:
        return None
    if pk >= n - 60:
        return None
    return pk


def autocorr_spacing(norm: np.ndarray, s0: int, s1: int, gsp: float, win: int = 1000,
                     step: int = 500, frac: float = 0.6, lo: float = 0.8, hi: float = 2.2):
    """Band spacing from the autocorrelation of the channel-summed processed
    signal (one pulse per base whatever its colour): the strongest local
    maximum of the autocorrelation within [lo, hi] x the peak-picking estimate
    `gsp`, median over 1000-scan windows in the first `frac` of the signal
    region.  None when no window shows a periodicity in that range."""
    from scipy.signal import find_peaks
    S = norm.sum(axis=1)
    a, b = int(np.floor(lo * gsp)), int(np.ceil(hi * gsp))
    est = []
    for w0 in range(s0, s0 + int(frac * (s1 - s0)) - win + 1, step):
        x = S[w0:w0 + win] - S[w0:w0 + win].mean()
        f = np.fft.rfft(x, 2 * win)
        ac = np.fft.irfft(f * np.conj(f))[:b + 2]
        if ac[0] <= 0:
            continue
        ac = ac / ac[0]
        pk, _ = find_peaks(ac[a:b + 1])
        if len(pk):
            k = a + pk[np.argmax(ac[a + pk])]
            y0, y1, y2 = ac[k - 1], ac[k], ac[k + 1]
            den = y0 - 2 * y1 + y2
            est.append(k + (0.5 * (y0 - y2) / den if den < 0 else 0.0))
    if not est:
        return None
    return float(np.median(est))


def reference_spacing(gsp: float, acr: float | None, ac_ratio: float = 1.2, ac_fallback: float = 1.4) -> float:
    """The spacing the tracker is anchored and clamped to.

    track_bases' peak-picking estimate `gsp` (median distance between envelope
    peaks) matches Cimarron's band spacing on clean reads, but on broad-peak
    reads it picks up sub-peaks and comes out 15-50% short, which makes the
    tracker call ~2 bases per band (mb4000: D02, M10, O02, O14 dropped from
    0.65-0.96 agreement to ~0 that way).  The autocorrelation estimate `acr`
    does not have that failure but reads ~5-20% long on clean reads, so the
    reference is max(gsp, acr / 1.2); when no periodicity is found at all,
    gsp is taken as 1.4x short."""
    if acr is None:
        return gsp * ac_fallback
    return max(gsp, float(acr) / ac_ratio)


def call_debug(trace: np.ndarray, base_order: str = "ACGT", start_mode: str = "onset",
               end_blob_offset: int = -60, end_margin: int = 240, pos_offset: float = 7.5,
               pos_drift: float = 0.0003, oracle_span: tuple | None = None, clamp=(0.7, 1.5),
               ac_ratio: float = 1.2, ac_fallback: float = 1.4):
    """call() with its knobs exposed (for development); returns (seq, positions, info)."""
    trace = np.asarray(trace, dtype=float)
    n = len(trace)
    info = {}
    z = prerun_z(trace)
    if not has_signal(z):
        info['empty'] = True
        return "", np.zeros(0, dtype=int), info
    bs, env, ch, norm = _process(trace)
    s0, s1 = detect_signal_region(bs)
    gsp = estimate_global_spacing(env, s0, s1)
    acr = autocorr_spacing(norm, s0, s1, gsp)
    ref = reference_spacing(gsp, acr, ac_ratio, ac_fallback)
    start = s0
    if start_mode == "onset":
        o = signal_onset(z)
        if o is not None:
            start = min(o, s0)
    end = n - end_margin
    pk = end_blob(trace)
    if pk is not None:
        end = min(end, pk + end_blob_offset)
    if oracle_span is not None:
        start, end = oracle_span
    info.update(start=start, end=end, blob=pk, s0=s0, s1=s1, gsp=gsp, acr=acr, ref=ref)
    bands = _track(env, ch, max(0, start - 2), s1, ref, prof_start=s0, prof_end=s1, clamp=tuple(clamp))
    bands = [b for b in bands if start <= b.position <= end]
    seq = "".join(base_order[b.channel] for b in bands)
    pos = np.array([b.position for b in bands], dtype=float)
    pos = pos + pos_offset + pos_drift * (pos - 2000.0)
    return seq, np.round(pos).astype(int), info


def call(trace: np.ndarray, base_order: str = "ACGT"):
    """Cimarron-3.12-like base calls: (sequence, positions in Cimarron's coordinates)."""
    seq, pos, _ = call_debug(trace, base_order)
    return seq, pos
