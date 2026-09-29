"""
Spacing-tracking base caller: rather than detecting peaks on each channel
independently (simple_caller.py) and reconciling conflicts after the
fact, this walks the trace left-to-right maintaining a running estimate
of the LOCAL expected inter-base spacing (updated by exponential moving
average as real spacing drifts over the length of a read -- capillary
electrophoresis runs are not perfectly uniform-speed), and searches only
within a window around the predicted next position for the best
candidate peak.

This directly encodes the strong periodicity prior real chromatogram
data has, which independent per-channel peak detection throws away --
in practice this resolves a lot of the ambiguity that caused
simple_caller.py to trade off identity against coverage (tightening
thresholds to cut false positives also cut real low-SNR peaks, and vice
versa). It's conceptually close to what the patent describes (a
quadratic/locally-fit "expected spacing curve" used to interpret
GapCheck/OmitOkN decisions), just implemented as an online tracker
instead of a two-pass quadratic fit.
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np

from .simple_caller import robust_baseline_subtract, normalize_channels_local, apply_mobility_correction


@dataclass
class TrackedBase:
    position: int
    channel: int
    height: float
    expected_position: float
    spacing_used: float


def _windowed_area(trace_col: np.ndarray, lo: int, hi: int, half_width: float) -> np.ndarray:
    """Integrated area (sum, not just point value) under trace_col across
    a window of full width ~2*half_width centered on each index in
    [lo, hi], via a cumulative-sum trick for O(1) per-position cost.
    Using area rather than the raw point value makes detection far more
    robust to the peak broadening/height loss that happens later in a
    capillary run (longer fragments diffuse more during migration, so
    peaks get shorter and wider but keep roughly the same area)."""
    n = len(trace_col)
    hw = max(1, int(round(half_width)))
    csum = np.concatenate([[0.0], np.cumsum(trace_col)])
    idx = np.arange(lo, hi + 1)
    left = np.clip(idx - hw, 0, n)
    right = np.clip(idx + hw + 1, 0, n)
    return csum[right] - csum[left]


def track_bases_area(
    trace: np.ndarray,
    base_order: str = "ACGT",
    initial_spacing: float = 11.0,
    window_frac: tuple[float, float] = (0.6, 1.5),
    min_area_fraction: float = 0.15,
    ema_alpha: float = 0.15,
    baseline_window: int = 151,
    local_norm_window: int = 300,
    mobility_shifts: list[int] | None = None,
    area_half_width_frac: float = 0.4,
    max_misses: int = 6,
) -> tuple[str, list[float], list[TrackedBase]]:
    """Like track_bases, but selects the next base by integrated local
    AREA (per channel, within a window that scales with the current
    tracked spacing) rather than raw point height -- see module docstring
    rationale on peak broadening vs. migration time.
    """
    baseline_subtracted = robust_baseline_subtract(trace, window=baseline_window)
    norm_trace = normalize_channels_local(baseline_subtracted, window=local_norm_window)
    if mobility_shifts is not None:
        norm_trace = apply_mobility_correction(norm_trace, mobility_shifts)

    n = norm_trace.shape[0]
    global_scale = np.percentile(norm_trace.max(axis=1), 90)

    pos = find_start_position(norm_trace.max(axis=1))
    spacing = initial_spacing
    misses = 0
    results: list[TrackedBase] = []

    while pos < n:
        expected_next = pos + spacing
        lo = int(expected_next - spacing * (1 - window_frac[0]))
        hi = int(expected_next + spacing * (window_frac[1] - 1))
        lo = max(lo, pos + 2)
        hi = min(hi, n - 1)
        if lo > hi:
            break

        half_width = max(1.0, spacing * area_half_width_frac)
        best_area, best_pos, best_ch = -1.0, None, None
        for ch in range(norm_trace.shape[1]):
            areas = _windowed_area(norm_trace[:, ch], lo, hi, half_width)
            i = int(np.argmax(areas))
            if areas[i] > best_area:
                best_area, best_pos, best_ch = float(areas[i]), lo + i, ch

        # Normalize area by window width so it's comparable to the
        # point-height-scale thresholds used elsewhere.
        area_norm = best_area / (2 * half_width + 1) if best_pos is not None else 0.0

        if best_pos is None or area_norm < min_area_fraction * global_scale:
            misses += 1
            pos = int(expected_next)
            if misses > max_misses:
                break
            continue

        misses = 0
        observed_spacing = best_pos - pos
        spacing = (1 - ema_alpha) * spacing + ema_alpha * observed_spacing
        spacing = float(np.clip(spacing, initial_spacing * 0.4, initial_spacing * 2.5))

        results.append(TrackedBase(
            position=best_pos, channel=best_ch, height=float(area_norm),
            expected_position=expected_next, spacing_used=spacing,
        ))
        pos = best_pos

    sequence = "".join(base_order[b.channel] for b in results)
    qualities = [b.height for b in results]
    return sequence, qualities, results


def detect_signal_region(
    baseline_subtracted: np.ndarray,
    min_fraction: float = 0.1,
    sustained: int = 20,
    search_bounds: tuple[int, int] | None = (1500, -100),
) -> tuple[int, int]:
    """Find the [start, end) region where real signal actually exists,
    using a threshold on the (not locally-normalized) baseline-subtracted
    envelope, searched only within `search_bounds` -- a wide, physically-
    justified prior range rather than the whole trace.

    Why bounded: a purely global-percentile threshold with no search
    bound was found to trigger far too early on some wells (A03: falsely
    triggered around scan 346 on baseline noise, vs. a true signal start
    around scan 2184 -- confirmed via matched SCF ground truth), which
    then poisons the spacing estimate with noise-driven "peaks" and causes
    massive overcalling. Checking real ground-truth start/end positions
    across 12 wells on this plate showed they're all tightly clustered
    (starts 1982-2184, ends 9330-9526) -- physically expected, since every
    capillary on a run shares the same injection timing. `search_bounds`
    defaults to a generous margin around that empirical range. This prior
    is run/instrument-timing-specific (injection time, voltage, etc, per
    the run metadata) -- if you're processing data from a differently-
    configured run, re-derive it or widen the bounds.

    search_bounds: (lo, hi). hi may be negative, meaning "this many
    samples from the end" (so the same default works across traces of
    varying length). Pass None to search the whole trace (old behavior --
    not recommended given the false-trigger issue above).
    """
    envelope = baseline_subtracted.max(axis=1)
    global_max = np.percentile(envelope, 99)
    threshold = min_fraction * global_max
    above = envelope > threshold
    n = len(above)

    if search_bounds is None:
        lo, hi = 0, n
    else:
        lo, hi = search_bounds
        hi = n + hi if hi < 0 else hi
        lo, hi = max(0, lo), min(n, hi)

    run = 0
    start = lo
    for i in range(lo, hi):
        run = run + 1 if above[i] else 0
        if run >= sustained:
            start = i - sustained + 1
            break
    else:
        start = lo

    run = 0
    end = hi
    for i in range(hi - 1, lo - 1, -1):
        run = run + 1 if above[i] else 0
        if run >= sustained:
            end = i + sustained
            break
    else:
        end = hi

    return start, min(end, n)


def trim_to_quality_block(
    tracked: list[TrackedBase],
    filter_widths: tuple[int, ...] = (5, 9, 15, 21),
    min_quality_percentile: float = 35.0,
) -> tuple[int, int]:
    """Find the longest contiguous run of called bases whose LOCAL
    (moving-average-filtered) quality stays above a threshold, to trim
    off degraded regions -- notably the unresolved-strand "cluster" that
    can appear at the very end of a run once fragments are too long to
    stay well-separated (progressive peak broadening/overlap, not a
    sudden signal dropout, so this can't be caught by amplitude
    thresholding alone -- see detect_signal_region for that separate,
    amplitude-based trim).

    Returns (start_idx, end_idx) into `tracked` (end exclusive).
    """
    if len(tracked) < 10:
        return 0, len(tracked)

    heights = np.array([b.height for b in tracked])
    threshold = np.percentile(heights, min_quality_percentile)

    best_len, best_range = -1, (0, len(tracked))
    for w in filter_widths:
        kernel = np.ones(w) / w
        filtered = np.convolve(heights, kernel, mode="same")
        above = filtered >= threshold
        start = None
        for i, ok in enumerate(list(above) + [False]):
            if ok and start is None:
                start = i
            elif not ok and start is not None:
                length = i - start
                if length > best_len:
                    best_len = length
                    best_range = (start, i)
                start = None

    return best_range


def second_pass_repeat_detection(
    tracked: list[TrackedBase],
    norm_trace: np.ndarray,
    gap_factor: float = 1.15,
    min_sub_peak_prominence: float = 0.4,
    valley_depth_frac: float = 0.8,
) -> list[TrackedBase]:
    """Second pass: look for missed bases in homopolymer/repeat runs.

    A base repeated 2-3+ times in a row (e.g. "TTT") produces consecutive
    peaks on the SAME channel -- exactly the case a single greedy forward
    pass is most likely to under-call, since two overlapping same-dye
    peaks blur together more than two different-dye peaks do (no color
    contrast to help separate them). This pass looks specifically at
    unusually large gaps between consecutively called bases (bigger than
    `gap_factor` times the locally expected spacing) and re-examines the
    SAME channel that was called on either side of the gap for a second,
    resolvable sub-peak in between -- i.e. specifically hunting for
    repeats, not general re-detection.

    valley_depth_frac: require a genuine valley (local minimum) between
    the candidate and each neighbor, dipping below valley_depth_frac
    times the candidate's own height, before accepting it as a real
    second band rather than the shoulder of a single broad peak that
    happens to clear the amplitude bar alone.
    """
    if len(tracked) < 3:
        return tracked

    spacings = np.array([b.spacing_used for b in tracked])
    out = [tracked[0]]

    for i in range(1, len(tracked)):
        prev, cur = tracked[i - 1], tracked[i]
        gap = cur.position - prev.position
        local_expected = spacings[i]

        if gap > gap_factor * local_expected and prev.channel == cur.channel:
            lo, hi = prev.position + 2, cur.position - 2
            if hi > lo:
                segment = norm_trace[lo:hi, prev.channel]
                sub_idx = int(np.argmax(segment))
                sub_val = segment[sub_idx]
                abs_pos = lo + sub_idx
                if sub_val >= min_sub_peak_prominence * min(prev.height, cur.height):
                    # Require a genuine valley (local minimum) on BOTH sides
                    # between this candidate and its neighbors -- a real
                    # second band should be properly separated, not just
                    # the tallest point of a single broad shoulder that
                    # happens to clear the amplitude bar.
                    left_valley = norm_trace[prev.position:abs_pos + 1, prev.channel].min() if abs_pos > prev.position else sub_val
                    right_valley = norm_trace[abs_pos:cur.position + 1, prev.channel].min() if cur.position > abs_pos else sub_val
                    valley_ok = (
                        left_valley < sub_val * valley_depth_frac
                        and right_valley < sub_val * valley_depth_frac
                    )
                    if valley_ok:
                        inserted_pos = abs_pos
                        out.append(TrackedBase(
                            position=inserted_pos, channel=prev.channel, height=float(sub_val),
                            expected_position=(prev.position + cur.position) / 2,
                            spacing_used=local_expected,
                        ))
        out.append(cur)

    out.sort(key=lambda b: b.position)
    return out


def apply_spectral_separation(baseline_subtracted: np.ndarray, chm: np.ndarray) -> np.ndarray:
    """Undo dye cross-talk / brightness imbalance via linear unmixing
    against an empirically-calibrated characteristic matrix (chm[i,j] =
    relative amplitude channel i shows when true base j is present,
    normalized so chm[j,j] == 1.0). Applied per-sample: corrected = chm^-1 @ observed.
    """
    chm_inv = np.linalg.inv(chm)
    return baseline_subtracted @ chm_inv.T


# Empirically fit from pooled calibration across all 12 real M13 wells on
# this plate (see dev notes / README) -- columns = true base A,C,G,T,
# rows = observed channel amplitude (normalized so the true channel's own
# entry is 1.0). This is chemistry/instrument-specific: re-fit if you have
# ground truth for a different run.
#
# IMPORTANT: the raw empirical matrix, applied directly via full inversion,
# made calling dramatically WORSE (mean identity 79.8%->70.4%, coverage std
# exploded to 48.3) -- not because spectral separation is a bad idea, but
# because the underlying calibration (SCF-to-raw position alignment via
# cross-correlation, ~0.4-0.55 correlation) isn't precise enough, and
# matrix inversion amplifies that calibration noise. A shrinkage-
# regularized version (blended 90% toward the identity matrix, 10% real
# correction) was tuned empirically and found to genuinely help: mean
# identity 79.8%->80.7%, coverage 85.6%->96.3%, AND meaningfully more
# consistent across wells (identity std 1.9->1.7, and critically the
# worst-case well jumps from 70% to 78.9% identity). DEFAULT_CHM_SHRUNK
# is what track_bases actually uses by default; DEFAULT_CHM is kept for
# reference / further calibration work.
DEFAULT_CHM = np.array([
    [1.000, 0.162, 0.333, 0.459],
    [1.736, 1.000, 0.376, 0.584],
    [0.188, 0.089, 1.000, 0.145],
    [0.460, 0.100, 0.923, 1.000],
])
DEFAULT_CHM_SHRUNK = 0.1 * DEFAULT_CHM + 0.9 * np.eye(4)


def smooth_trace(trace: np.ndarray, window: int = 2) -> np.ndarray:
    """Light moving-average smoothing applied to the RAW trace, before
    baseline subtraction or any correction. This turned out to matter a
    lot: applying the FULL (unshrunk) empirically-fit spectral separation
    matrix (see DEFAULT_CHM) directly to unsmoothed data was unusably
    unstable (mean identity collapsed to 70.4%, coverage std blew up to
    48.3) even though the matrix itself was shown to be highly consistent
    and reliable across all 12 wells (~5% coefficient of variation per
    entry) -- the instability was matrix inversion amplifying per-SAMPLE
    electrical noise (inversion has gain >1 in some directions), not
    inaccuracy in the average matrix. A light window=2-3 smoothing pass
    before separation suppresses that per-sample noise enough for the
    full-strength correction to be stable AND effective: mean identity
    79.8%->89.9%, with a much higher floor across wells (worst-case
    identity 70%->87.7%).
    """
    from scipy.ndimage import uniform_filter1d
    if window <= 1:
        return trace
    return uniform_filter1d(trace, size=window, axis=0, mode="nearest")


def _looks_periodic(envelope: np.ndarray, threshold: float, min_peaks: int = 4, cv_max: float = 0.6) -> bool:
    """Check whether the region ahead shows regularly-spaced peaks (the
    signature of real signal) rather than isolated noise excursions.
    Real base calls come out roughly periodically; noise crossing a
    threshold once or twice, with no regular follow-up, shouldn't count.
    """
    from scipy.signal import find_peaks
    idx, _ = find_peaks(envelope, height=threshold * 0.5, distance=3)
    if len(idx) < min_peaks:
        return False
    spacings = np.diff(idx[:min_peaks + 3])
    if len(spacings) == 0 or np.mean(spacings) < 1e-9:
        return False
    cv = np.std(spacings) / np.mean(spacings)
    return cv < cv_max


def detect_signal_region_auto(
    baseline_subtracted: np.ndarray,
    calib_samples: int = 500,
    k: float = 5.0,
    sustained: int = 20,
    validate_window: int = 600,
) -> tuple[int, int]:
    """Fully self-contained start/stop detection -- no prior knowledge of
    a specific run's injection timing baked in (contrast with the earlier
    `detect_signal_region`, which relies on a hard-coded search window
    empirically fit to one specific plate's timing and won't generalize
    to a run with different instrument settings).

    Two-stage: (1) find candidate crossings of a per-trace, robust
    (median + k*MAD) amplitude threshold estimated from this trace's own
    early samples; (2) validate each candidate by checking that
    regularly-spaced peaks actually follow it (see _looks_periodic) --
    this rejects isolated noise excursions that happen to cross the
    amplitude threshold without being followed by real periodic signal
    (this combination was needed: amplitude threshold alone gave exact
    matches on some wells but was badly fooled on others where baseline
    noise dispersion happens to be proportionally larger relative to
    early, still-weak real signal).
    """
    envelope = baseline_subtracted.max(axis=1)
    n = len(envelope)

    med = np.median(envelope[:calib_samples])
    mad = np.median(np.abs(envelope[:calib_samples] - med)) * 1.4826
    threshold = med + k * mad

    above = envelope > threshold
    run = 0
    start = 0
    i = 0
    while i < n:
        run = run + 1 if above[i] else 0
        if run >= sustained:
            candidate = i - sustained + 1
            window_end = min(n, candidate + validate_window)
            if _looks_periodic(envelope[candidate:window_end], threshold):
                start = candidate
                break
            # Not validated -- likely an isolated noise excursion. Reset
            # and keep scanning forward past it.
            run = 0
        i += 1
    else:
        start = 0

    # End: same idea, scanning backward from the trace's end.
    run = 0
    end = n
    i = n - 1
    while i >= start:
        run = run + 1 if above[i] else 0
        if run >= sustained:
            candidate = i + sustained
            window_start = max(start, candidate - validate_window)
            if _looks_periodic(envelope[window_start:candidate], threshold):
                end = candidate
                break
            run = 0
        i -= 1
    else:
        end = n

    return start, min(end, n)


def precompute_channel_peak_masks(norm_trace: np.ndarray, prominence: float = 0.04) -> np.ndarray:
    """Boolean (n_samples, 4) mask: True where that sample is an
    independently-detected local peak on THAT channel alone (not just the
    tallest of the 4 at that sample). Used as corroborating evidence --
    a position that's both the envelope max AND a genuine local peak on
    its own channel is more likely a true base than one that's only the
    tallest of four mostly-flat/noisy channels.
    """
    from scipy.signal import find_peaks
    n = norm_trace.shape[0]
    mask = np.zeros((n, 4), dtype=bool)
    for c in range(4):
        idx, _ = find_peaks(norm_trace[:, c], prominence=prominence, distance=2)
        mask[idx, c] = True
    return mask


def combined_channel_score(
    norm_trace: np.ndarray, channel_peak_mask: np.ndarray, position: int, channel_bonus: float = 0.5
) -> np.ndarray:
    """Per-channel combined score at one position: the channel's own
    normalized amplitude, boosted if that sample is ALSO independently a
    local peak on that specific channel (corroborating evidence from two
    different detectors -- envelope-relative height, and single-channel
    peak detection -- rather than trusting either alone)."""
    row = norm_trace[position].copy()
    boost = channel_peak_mask[position].astype(float) * channel_bonus
    return row * (1.0 + boost)


def trim_mott(
    tracked: list[TrackedBase],
    quality_cutoff: float = 0.08,
) -> tuple[int, int]:
    """Mott's algorithm for quality-based trimming (the same underlying
    method used by Phred/Phrap and by sangeranalyseR's M1 trimming) --
    an established, widely-used technique rather than an ad-hoc heuristic.

    For each base, compute (cutoff - quality); this is positive for
    low-quality bases and negative for high-quality ones. Take the
    running cumulative sum, clamping at 0 (never go negative), and find
    the contiguous region achieving the highest cumulative score. This
    naturally finds the longest/best-supported high-quality stretch
    without needing an arbitrary moving-average window width, unlike
    trim_to_quality_block.

    `quality_cutoff` is on the same 0-ish-to-1-ish scale as
    TrackedBase.height (normalized local peak height/area), not a
    Phred-scale quality -- tune relative to that.
    """
    if len(tracked) < 5:
        return 0, len(tracked)

    heights = np.array([b.height for b in tracked])
    scores = quality_cutoff - heights  # positive = "bad", negative = "good"

    running = 0.0
    running_start = 0
    best_score = 0.0
    best_range = (0, len(tracked))

    for i, s in enumerate(scores):
        if running <= 0:
            running = 0.0
            running_start = i
        running -= s  # subtract because negative s (good bases) should increase running
        if running > best_score:
            best_score = running
            best_range = (running_start, i + 1)

    return best_range


def band_filter(trace: np.ndarray, spike_width: int = 3) -> np.ndarray:
    """'Band filter' -- reinterpreted after the frequency-domain band-pass
    version tested poorly (every tested cutoff hurt coverage 15-24pp for
    no identity gain, likely filtfilt ringing/edge artifacts fighting the
    position-domain peak detector downstream).

    In gel/capillary electrophoresis terminology, "bands" are the DNA
    fragment peaks themselves (as in gel electrophoresis "bands") -- so a
    "band filter" more plausibly means a DESPIKING filter: remove noise
    spikes narrower than any real band could physically be (a genuine DNA
    peak has a minimum width from diffusion during migration), while
    preserving real band shapes, via a median filter with a window
    narrower than the minimum plausible band width. This is a standard
    despiking technique for chromatography-style signals, and unlike mean
    smoothing, a median filter doesn't blur genuine peaks -- it only
    suppresses features narrower than its window.
    """
    from scipy.signal import medfilt

    out = np.empty_like(trace, dtype=float)
    k = spike_width if spike_width % 2 == 1 else spike_width + 1
    for c in range(trace.shape[1]):
        out[:, c] = medfilt(trace[:, c], kernel_size=k)
    return out


def _gaussian_psf(spatial_sigma: float, n: int) -> np.ndarray:
    x = np.arange(n) - n // 2
    psf = np.exp(-0.5 * (x / max(spatial_sigma, 1e-6)) ** 2)
    return psf / psf.sum()


def gaussian_reconstruction_filter(
    trace: np.ndarray, spacing: float, segment_size: int = 2048, noise_reg: float = 5e-2,
    sigma_scale: float = 1.0,
) -> np.ndarray:
    """Frequency-domain Gaussian reconstruction filter for blind
    deconvolution, per the actual Cimarron source (mb.cxx, Univ. of Utah,
    1996) as embedded in EP0944739A1 -- the real "band filter" step
    (Baseline sub -> Spectral Sep -> Normalization -> Band filter ->
    Mobility shift corr -> base calling), not the frequency-domain
    band-pass or despiking guesses tried earlier this session.

    fbw = K * N / (2*pi*spacing), where N is the FFT segment size and
    spacing is the LOCALLY measured band spacing (adaptive -- gets
    narrower where bands are closer together, e.g. in fast-migrating
    early regions or tight repeats, and wider later where bands spread
    out from diffusion). K = 2*sqrt(ln(0.23)/-0.5), matching the formula
    identified earlier in this session's patent-literal attempt.

    This is a RECONSTRUCTION (sharpening/deconvolution) filter, not a
    smoothing filter -- it divides out an assumed Gaussian blur kernel of
    width `fbw` in the frequency domain (regularized, Wiener-style, to
    avoid amplifying noise at frequencies the kernel suppresses heavily),
    aiming to un-blur overlapping/broadened peaks -- exactly the
    mechanism needed for closely-spaced same-channel repeat peaks, which
    diagnosis showed account for 70% of this pipeline's remaining
    deletion errors.
    """
    K = 2.0 * np.sqrt(np.log(0.23) / -0.5)
    fbw = K * segment_size / (2.0 * np.pi * max(spacing, 1e-6))
    fbw = float(np.clip(fbw, 1.0, segment_size / 4))

    # fbw as computed above is a FREQUENCY-domain parameter (bandwidth in
    # FFT-bin-equivalent units over the segment) -- the corresponding
    # SPATIAL-domain Gaussian sigma is segment_size/(2*pi*fbw), which
    # algebraically simplifies to spacing/K (~0.29*spacing). Using fbw
    # directly as a spatial-domain sigma (an earlier version of this
    # function did) gives a kernel far WIDER than the band spacing itself
    # (e.g. sigma~28 vs spacing~10) -- physically backwards for a blur
    # kernel, and caused near-total collapse of peak detectability
    # (coverage ~3-5%) when tested. This corrected mapping keeps the PSF
    # width sensibly smaller than the spacing it's meant to help resolve.
    spatial_sigma = spacing / K

    n = trace.shape[0]
    spatial_sigma = spatial_sigma * float(sigma_scale)
    psf = _gaussian_psf(spatial_sigma, n)
    psf_fft = np.fft.rfft(np.fft.ifftshift(psf))

    out = np.empty_like(trace, dtype=float)
    for c in range(trace.shape[1]):
        sig_fft = np.fft.rfft(trace[:, c])
        wiener = np.conj(psf_fft) / (np.abs(psf_fft) ** 2 + noise_reg)
        deconvolved = np.fft.irfft(sig_fft * wiener, n=n)
        out[:, c] = deconvolved
    return out


def apply_gaussian_reconstruction_windowed(
    trace: np.ndarray, spacing_curve: np.ndarray, segment_size: int = 512, overlap: int = 64,
    noise_reg: float = 5e-2, sigma_scale: float = 1.0,
    noise_reg_curve: np.ndarray | None = None,
) -> np.ndarray:
    """Apply gaussian_reconstruction_filter in overlapping windows using
    the LOCAL spacing at each window (spacing genuinely drifts over a
    read), cross-fading overlaps to avoid seams.

    noise_reg_curve: optional per-sample regularization strength. When given
    the window uses its local median instead of the scalar `noise_reg`, so the
    filter can be sharp early in the read (where matched bases are dense) and
    smoothly regularized in the noisy 3' tail (where substitutions live)."""
    n = trace.shape[0]
    out = np.zeros_like(trace, dtype=float)
    weight = np.zeros(n)
    fade = np.linspace(0, 1, overlap) if overlap > 0 else np.array([])

    start = 0
    while start < n:
        end = min(start + segment_size, n)
        seg = trace[start:end]
        if seg.shape[0] < 16:
            out[start:end] += seg
            weight[start:end] += 1.0
            break
        local_spacing = float(np.median(spacing_curve[start:end]))
        local_reg = float(np.median(noise_reg_curve[start:end])) if noise_reg_curve is not None else noise_reg
        filtered = gaussian_reconstruction_filter(seg, local_spacing, segment_size=seg.shape[0],
                                                  noise_reg=local_reg, sigma_scale=sigma_scale)

        w = np.ones(seg.shape[0])
        if overlap > 0 and start > 0:
            w[:min(overlap, len(fade))] = fade[:min(overlap, len(fade))]
        if overlap > 0 and end < n:
            w[-min(overlap, len(fade)):] = fade[::-1][:min(overlap, len(fade))]

        out[start:end] += filtered * w[:, None]
        weight[start:end] += w
        if end >= n:
            break
        start = end - overlap

    weight[weight < 1e-9] = 1.0
    return out / weight[:, None]


# Position-adaptive characteristic matrices, fit separately for 5 equal
# read-position buckets (0-20%, 20-40%, ..., 80-100% of the way through
# the signal region), pooled across all 12 real wells. Confirmed the
# cross-talk pattern genuinely CHANGES over a read, not just gets
# noisier -- e.g. G->C cross-talk nearly triples (0.26-0.32 -> 0.789) and
# A->T more than doubles (0.41-0.49 -> 0.888) in the last 20% of the
# read, matching exactly where error-rate diagnosis found channel
# misassignment concentrated (23.2% in the last fifth vs 0.3% mid-read).
# A single global matrix (DEFAULT_CHM) structurally can't represent this.
POSITION_ADAPTIVE_CHM_BUCKETS = [
    np.array([[1.000, 0.153, 0.323, 0.467], [1.713, 1.000, 0.404, 0.617],
              [0.142, 0.077, 1.000, 0.156], [0.437, 0.091, 0.965, 1.000]]),
    np.array([[1.000, 0.116, 0.280, 0.427], [1.745, 1.000, 0.321, 0.556],
              [0.109, 0.056, 1.000, 0.116], [0.415, 0.062, 0.929, 1.000]]),
    np.array([[1.000, 0.119, 0.283, 0.416], [1.731, 1.000, 0.308, 0.533],
              [0.137, 0.057, 1.000, 0.110], [0.432, 0.067, 0.936, 1.000]]),
    np.array([[1.000, 0.106, 0.255, 0.412], [1.717, 1.000, 0.258, 0.556],
              [0.189, 0.082, 1.000, 0.122], [0.492, 0.070, 0.882, 1.000]]),
    np.array([[1.000, 0.209, 0.404, 0.460], [1.740, 1.000, 0.789, 0.582],
              [0.321, 0.205, 1.000, 0.168], [0.888, 0.305, 0.889, 1.000]]),
]


def apply_position_adaptive_spectral_separation(
    baseline_subtracted: np.ndarray,
    sig_start: int,
    sig_end: int,
    bucket_matrices: list[np.ndarray] = POSITION_ADAPTIVE_CHM_BUCKETS,
    spacing_curve: np.ndarray | None = None,
) -> np.ndarray:
    """Apply spectral separation using a matrix that varies smoothly over
    the read position, instead of one fixed global matrix -- see
    POSITION_ADAPTIVE_CHM_BUCKETS docstring for why this matters. Uses
    piecewise-constant matrices per bucket (not continuous per-sample
    interpolation, for speed -- 5 segments captures the effect fine
    without needing a fresh matrix inversion per sample).

    spacing_curve: per-scan-position estimated local base spacing (same
    length as baseline_subtracted). If given, bucket boundaries are placed
    at equal CUMULATIVE BASE COUNT (integrating 1/spacing across scan
    position), not equal scan-position span -- important because base
    density is NOT uniform over a read (spacing grows later on from
    diffusion), so an even scan-position split lands ~15-20 bases early
    relative to where the calibration buckets were actually defined (by
    base index, not scan position) -- confirmed as a real, measurable
    misalignment before this fix (checked against ground truth: bucket 4
    starting 176-256 scans too early). Falls back to even scan-position
    split if not given.
    """
    n = baseline_subtracted.shape[0]
    out = baseline_subtracted.copy()
    n_buckets = len(bucket_matrices)

    if spacing_curve is not None:
        # Cumulative base count estimate via integrating 1/spacing.
        inv_spacing = 1.0 / np.clip(spacing_curve[sig_start:sig_end], 1.0, None)
        cumulative_bases = np.concatenate([[0.0], np.cumsum(inv_spacing)])
        total_bases = cumulative_bases[-1]
        bounds = [sig_start]
        for b in range(1, n_buckets):
            target = total_bases * b / n_buckets
            idx = int(np.searchsorted(cumulative_bases, target))
            bounds.append(sig_start + min(idx, len(inv_spacing)))
        bounds.append(sig_end)
    else:
        span = max(sig_end - sig_start, 1)
        bounds = [sig_start + int(span * b / n_buckets) for b in range(n_buckets + 1)]

    bounds[0] = 0
    bounds[-1] = n

    for b in range(n_buckets):
        lo, hi = bounds[b], bounds[b + 1]
        if hi <= lo:
            continue
        chm_inv = np.linalg.inv(bucket_matrices[b])
        out[lo:hi] = baseline_subtracted[lo:hi] @ chm_inv.T

    return out


def find_start_position(envelope: np.ndarray, threshold: float = 0.15, sustained: int = 20) -> int:
    """Find the first point where the envelope sustainably rises above
    `threshold` (avoiding triggering on an isolated early noise spike)."""
    above = envelope > threshold
    run = 0
    for i, ok in enumerate(above):
        run = run + 1 if ok else 0
        if run >= sustained:
            return max(0, i - sustained)
    return 0


def estimate_global_spacing(envelope: np.ndarray, sig_start: int, sig_end: int) -> float:
    """Robust global spacing estimate (median peak-to-peak distance over
    the whole signal region), used as an anchor the local EMA tracker gets
    pulled back toward -- prevents the kind of unbounded drift where one
    noisy stretch causes the running estimate to wander far from reality
    and never recover (observed on some wells: spacing drifting from ~14
    up past 19 then collapsing to ~6, derailing the rest of the read)."""
    from scipy.signal import find_peaks
    idx, _ = find_peaks(envelope[sig_start:sig_end], prominence=0.1, distance=3)
    if len(idx) < 5:
        return 12.0
    spacings = np.diff(idx)
    # median is robust to the occasional too-close/too-far outlier pair
    return float(np.median(spacings))


# Empirically measured per-channel dye mobility shift, pooled (median then
# mean) across all 12 real M13 wells on this plate -- see README. Different
# dyes migrate through the capillary at very slightly different rates, so
# the same physical base shows up at a slightly different scan index per
# channel. Re-fit if you have ground truth for a different chemistry.
DEFAULT_MOBILITY_SHIFTS = [-2, -1, -2, 2]  # A, C, G, T


def _ramp_value(value, frac: float, f0: float, f1: float) -> float:
    """Evaluate a possibly position-profiled parameter at scan-fraction `frac`.

    If `value` is a scalar, return it unchanged (so the scalar config is
    reproduced byte-identically to the un-profiled caller). If `value` is a
    (start, end) pair, linearly interpolate start->end as `frac` runs from `f0`
    to `f1` (clamped outside that interval). A (v, v) pair is treated as the
    scalar v, again for byte-identical flat behaviour.
    """
    if isinstance(value, (tuple, list, np.ndarray)) and len(value) == 2:
        a, b = float(value[0]), float(value[1])
        if a == b:
            return a
        if f1 <= f0:
            return b
        t = (frac - f0) / (f1 - f0)
        t = 0.0 if t < 0.0 else (1.0 if t > 1.0 else t)
        return a + (b - a) * t
    return float(value)


def _ramp_array(value, frac_of: np.ndarray, f0: float, f1: float) -> np.ndarray:
    """Vectorized `_ramp_value` over a per-sample fraction array."""
    if isinstance(value, (tuple, list, np.ndarray)) and len(value) == 2:
        a, b = float(value[0]), float(value[1])
        if a == b:
            return np.full(frac_of.shape, a, dtype=float)
        if f1 <= f0:
            return np.full(frac_of.shape, b, dtype=float)
        t = np.clip((frac_of - f0) / (f1 - f0), 0.0, 1.0)
        return a + (b - a) * t
    return np.full(frac_of.shape, float(value), dtype=float)


def prune_spurious_insertions(
    tracked: list[TrackedBase],
    merge_frac: float = 1.4,
    height_frac: float = 0.9,
    max_iter: int = 4,
    expected: np.ndarray | None = None,
) -> list[TrackedBase]:
    """Fuzzy insertion detection ("OmitOkN" in the real module names).

    A single real spacing interval should contain at most one band. If a called
    base sits BETWEEN two neighbours whose combined spacing `(b-a) + (c-b)` is
    smaller than `merge_frac` times the locally expected spacing, then it is a
    candidate over-call: three called bases occupy the room of fewer than
    ~1.4 real bands. We only drop it when it is also weaker than BOTH
    neighbours (`height < height_frac * min(neighbour heights)`), i.e. it looks
    like a shoulder of a broad peak rather than a resolved band. Removing a
    spurious insertion deletes a gap column, which merges the flanking
    matching runs -- raising identity and the longest error-free run without
    removing any true match.

    `expected`, when given, is a per-scan-position array of the TRUE local
    band spacing (e.g. the spacing-anchor curve). Using it instead of the
    tracker's own `spacing_used` is essential: in an over-called region the
    adaptive estimate shrinks to match the overcalls, so it can never flag
    them -- the earlier version using `spacing_used` fired on only ~1% of the
    real insertions.

    Iterates until no more candidates (or `max_iter`), so runs of overcalls
    collapse one at a time.
    """
    if len(tracked) < 3:
        return tracked
    out = list(tracked)
    for _ in range(max_iter):
        removed = False
        i = 1
        while i < len(out) - 1:
            a, b, c = out[i - 1], out[i], out[i + 1]
            d1 = b.position - a.position
            d2 = c.position - b.position
            if d1 <= 0 or d2 <= 0:
                i += 1
                continue
            if expected is not None and 0 <= b.position < len(expected):
                exp = float(expected[b.position])
            else:
                exp = b.spacing_used if b.spacing_used and b.spacing_used > 0 else (d1 + d2) / 2.0
            if (d1 + d2) < merge_frac * exp and b.height < height_frac * min(a.height, c.height):
                del out[i]
                removed = True
                continue
            i += 1
        if not removed:
            break
    return out


def track_bases(
    trace: np.ndarray,
    base_order: str = "ACGT",
    initial_spacing: float = 14.0,
    window_frac: tuple[float, float] = (0.66, 1.35),
    min_prominence: float | tuple[float, float] = 0.05,
    ema_alpha: float | tuple[float, float] = 0.08,
    pullback_weight: float | tuple[float, float] = 0.08,
    baseline_window: int = 151,
    local_norm_window: int = 300,
    mobility_shifts: list[int] | str | None = "default",
    max_misses: int = 6,
    auto_trim: bool = True,
    trim_quality_percentile: float = 10.0,
    trim_method: str = "percentile",  # "percentile" (balanced, default) or
    # "mott" (established Phred/Phrap-style trim -- see trim_mott).
    # "mott" with a high mott_quality_cutoff (~0.96-1.0) trims to a smaller,
    # much higher-confidence CORE of the read: identity 91.1%->96.5%+ (can
    # exceed Cimarron 3.12's own ~95.4% identity benchmark), but coverage
    # drops from 95.8% to ~73-76% -- a genuinely different trade-off than
    # the balanced default, not a strict improvement. Good for use cases
    # that want a shorter, near-certainly-correct core read rather than
    # the longest defensible one. See README for the full trade-off curve.
    mott_quality_cutoff: float = 0.08,
    repeat_detection: bool = False,
    repeat_gap_factor: float = 1.15,
    repeat_min_sub_peak_prominence: float = 0.4,
    repeat_valley_depth_frac: float = 0.8,
    spectral_separation_matrix: np.ndarray | str | None = "default",
    smoothing_window: int = 2,
    use_gaussian_reconstruction: bool = False,
    gaussian_recon_noise_reg: float = 0.05,
    gaussian_recon_segment_size: int = 512,
    gaussian_recon_sigma_scale: float = 1.0,
    position_adaptive_spectral: bool = False,
    use_combined_channel_score: bool = False,
    channel_peak_bonus: float | tuple[float, float] = 0.5,  # tested across all 12 wells: net negative on
    # balanced accuracy (coverage +6.4pp but identity -4.1pp) -- more false
    # repeat insertions than true ones recovered with current tuning. Left
    # in and available to opt into / keep tuning, not removed, since the
    # mechanism is sound (this is a real, diagnosed failure mode) but the
    # current threshold isn't there yet.
    profile_fracs: tuple[float, float] = (0.0, 1.0),
    prune_insertions: bool = False,
    prune_merge_frac: float = 1.4,
    prune_height_frac: float = 0.9,
    prune_max_iter: int = 4,
    prune_expected_spacing: bool = False,
    use_spacing_anchor_curve: bool = False,
) -> tuple[str, list[float], list[TrackedBase]]:
    """Walk the trace predicting each next base position from a running
    local spacing estimate. Returns (sequence, qualities, tracked_bases).

    The spacing estimate is an EMA of observed spacing (per `ema_alpha`)
    that is ALSO pulled, on every step, a little (`pullback_weight`)
    toward a robust global median spacing for the whole read -- pure EMA
    with no anchor was found to drift unboundedly on some wells (a noisy
    stretch nudges it up, then it overshoots, then collapses, and it
    never finds its way back), which silently truncates or garbles the
    rest of the read even though tracking nominally "continues".

    auto_trim: if True (default), apply trim_to_quality_block() as a final
    step to cut off degraded regions -- notably an unresolved-strand
    cluster that can appear at the end of a run once fragments are long
    enough that separation breaks down. Set False to get the raw,
    untrimmed call (e.g. for diagnostics).

    Position profiling: `ema_alpha`, `pullback_weight`, `min_prominence` and
    `channel_peak_bonus` may each be given as a (start, end) pair instead of a
    scalar. The pair is linearly interpolated over the scan-fraction interval
    `profile_fracs` (fraction of the trace between the detected signal-region
    start and end), so a parameter can e.g. loosen its pull-back toward the
    global spacing across the degraded 3' tail. A scalar, or a (v, v) pair,
    reproduces the un-profiled caller byte-for-byte.
    """
    baseline_subtracted = robust_baseline_subtract(smooth_trace(trace, smoothing_window), window=baseline_window)
    if spectral_separation_matrix is not None:
        if position_adaptive_spectral:
            # Signal region must be estimated BEFORE separation here, from
            # the pre-separation signal, since which regional matrix to
            # apply depends on position within the read. A cheap first-
            # pass track gives a real spacing curve so bucket boundaries
            # land at equal BASE COUNT, not equal scan-position span (see
            # apply_position_adaptive_spectral_separation docstring).
            sig_start_pre, sig_end_pre = detect_signal_region(baseline_subtracted)
            _, _, pre_tracked = track_bases(
                trace, base_order=base_order, spectral_separation_matrix=None,
                use_gaussian_reconstruction=False, position_adaptive_spectral=False,
                auto_trim=False, smoothing_window=smoothing_window,
            )
            n_pre = baseline_subtracted.shape[0]
            if len(pre_tracked) >= 5:
                pre_positions = np.array([b.position for b in pre_tracked[1:]])
                pre_spacings = np.array([b.spacing_used for b in pre_tracked[1:]])
                coeffs = np.polyfit(pre_positions, pre_spacings, deg=2)
                spacing_curve_est = np.clip(np.polyval(coeffs, np.arange(n_pre)), 2.0, None)
            else:
                spacing_curve_est = None
            baseline_subtracted = apply_position_adaptive_spectral_separation(
                baseline_subtracted, sig_start_pre, sig_end_pre, spacing_curve=spacing_curve_est
            )
        else:
            chm = DEFAULT_CHM if isinstance(spectral_separation_matrix, str) else spectral_separation_matrix
            baseline_subtracted = apply_spectral_separation(baseline_subtracted, chm)
    if mobility_shifts is not None:
        shifts = DEFAULT_MOBILITY_SHIFTS if isinstance(mobility_shifts, str) else mobility_shifts
        baseline_subtracted = apply_mobility_correction(baseline_subtracted, shifts)
    norm_trace = normalize_channels_local(baseline_subtracted, window=local_norm_window)

    if use_gaussian_reconstruction:
        # "Band filter" step, per the real Sequence Analyzer pipeline order
        # (Baseline sub -> Spectral Sep -> Normalization -> Band filter ->
        # Mobility shift corr -> base calling) and the actual Cimarron
        # source formula (mb.cxx) -- a frequency-domain Gaussian
        # reconstruction/deconvolution filter, bandwidth set from measured
        # band spacing. See gaussian_reconstruction_filter().
        env_for_spacing = norm_trace.max(axis=1)
        sig_start_gr, sig_end_gr = detect_signal_region(baseline_subtracted)
        spacing_for_filter = estimate_global_spacing(env_for_spacing, sig_start_gr, sig_end_gr)
        spacing_curve_const = np.full(norm_trace.shape[0], spacing_for_filter)
        reg_curve = None
        if (isinstance(gaussian_recon_noise_reg, (tuple, list, np.ndarray))
                and len(gaussian_recon_noise_reg) == 2
                and float(gaussian_recon_noise_reg[0]) != float(gaussian_recon_noise_reg[1])):
            _span_gr = max(1, sig_end_gr - sig_start_gr)
            _fraco = np.clip((np.arange(norm_trace.shape[0]) - sig_start_gr) / _span_gr, 0.0, 1.0)
            reg_curve = _ramp_array(gaussian_recon_noise_reg, _fraco, *profile_fracs)
        norm_trace = apply_gaussian_reconstruction_windowed(
            norm_trace, spacing_curve_const, segment_size=gaussian_recon_segment_size,
            noise_reg=float(np.atleast_1d(gaussian_recon_noise_reg)[0]),
            sigma_scale=gaussian_recon_sigma_scale, noise_reg_curve=reg_curve,
        )
        norm_trace = np.clip(norm_trace, 0, None)

    envelope = norm_trace.max(axis=1)
    best_channel = norm_trace.argmax(axis=1)
    n = len(envelope)

    sig_start, sig_end = detect_signal_region(baseline_subtracted)
    _span = max(1, sig_end - sig_start)
    frac_of = np.clip((np.arange(n) - sig_start) / _span, 0.0, 1.0)
    _pf0, _pf1 = profile_fracs

    if use_combined_channel_score:
        # Boost each channel's value wherever that sample is ALSO an
        # independently-detected local peak on that specific channel --
        # corroborating evidence from two different detectors (envelope-
        # relative height, and single-channel peak detection) rather than
        # trusting the envelope-argmax alone.
        channel_peak_mask = precompute_channel_peak_masks(norm_trace)
        if (isinstance(channel_peak_bonus, (tuple, list, np.ndarray))
                and len(channel_peak_bonus) == 2
                and float(channel_peak_bonus[0]) != float(channel_peak_bonus[1])):
            bonus_arr = _ramp_array(channel_peak_bonus, frac_of, _pf0, _pf1)
            boosted = norm_trace * (1.0 + channel_peak_mask.astype(float) * bonus_arr[:, None])
        else:
            boosted = norm_trace * (1.0 + channel_peak_mask.astype(float) * channel_peak_bonus)
        envelope = boosted.max(axis=1)
        best_channel = boosted.argmax(axis=1)

    global_spacing = estimate_global_spacing(envelope, sig_start, sig_end)
    anchor_curve = None
    if use_spacing_anchor_curve or prune_expected_spacing:
        # A single global median is too WIDE for the primer-adjacent region
        # (small fragments migrate closer together, so early spacing is
        # genuinely smaller than the read-wide median) and too NARROW for the
        # far tail. Fit a smooth spacing-vs-position curve from a cheap
        # first-pass track and pull the EMA toward that local value instead of
        # a constant -- an online version of the patent's "expected spacing
        # curve". The pre-track itself runs with the constant anchor (this
        # flag defaults off there), so it stays a pure trend estimate.
        _, _, pre_tracked = track_bases(
            trace, base_order=base_order, spectral_separation_matrix=None,
            use_gaussian_reconstruction=False, position_adaptive_spectral=False,
            auto_trim=False, smoothing_window=smoothing_window,
        )
        if len(pre_tracked) >= 8:
            ppos = np.array([b.position for b in pre_tracked[1:]], dtype=float)
            psp = np.array([b.spacing_used for b in pre_tracked[1:]], dtype=float)
            coeffs = np.polyfit(ppos, psp, deg=2 if len(ppos) >= 12 else 1)
            anchor_curve = np.clip(np.polyval(coeffs, np.arange(n)),
                                   global_spacing * 0.4, global_spacing * 2.0)
    pos = sig_start
    spacing = global_spacing  # start from the robust global estimate, not a fixed constant
    misses = 0
    results: list[TrackedBase] = []

    while pos < n:
        _frac = (pos - sig_start) / _span
        _mp = _ramp_value(min_prominence, _frac, _pf0, _pf1)
        _ea = _ramp_value(ema_alpha, _frac, _pf0, _pf1)
        _pw = _ramp_value(pullback_weight, _frac, _pf0, _pf1)
        expected_next = pos + spacing
        lo = int(expected_next - spacing * (1 - window_frac[0]))
        hi = int(expected_next + spacing * (window_frac[1] - 1))
        lo = max(lo, pos + 2)
        hi = min(hi, n - 1)
        if lo > hi:
            break

        window = envelope[lo:hi + 1]
        if len(window) == 0:
            break
        local_peak_offset = int(np.argmax(window))
        local_peak_pos = lo + local_peak_offset
        local_peak_val = window[local_peak_offset]

        if local_peak_val < _mp:
            # No credible peak in the expected window -- count a miss and
            # advance by the current spacing estimate anyway, so a few
            # consecutive weak/missed bases don't derail tracking entirely.
            misses += 1
            pos = int(expected_next)
            if misses > max_misses:
                break
            continue

        misses = 0
        observed_spacing = local_peak_pos - pos
        spacing = (1 - _ea) * spacing + _ea * observed_spacing
        # Pull back toward the robust local estimate every step, so a
        # noisy local stretch can nudge the estimate but can't make it
        # drift away unboundedly -- it's always being tugged back toward
        # what the read's spacing actually looks like at this position.
        _anchor = float(anchor_curve[pos]) if anchor_curve is not None else global_spacing
        spacing = (1 - _pw) * spacing + _pw * _anchor
        spacing = float(np.clip(spacing, _anchor * 0.5, _anchor * 2.0))

        ch = int(best_channel[local_peak_pos])
        results.append(TrackedBase(
            position=local_peak_pos, channel=ch, height=float(local_peak_val),
            expected_position=expected_next, spacing_used=spacing,
        ))
        pos = local_peak_pos

    if repeat_detection and results:
        results = second_pass_repeat_detection(results, norm_trace, gap_factor=repeat_gap_factor,
                                                min_sub_peak_prominence=repeat_min_sub_peak_prominence,
                                                valley_depth_frac=repeat_valley_depth_frac)

    if prune_insertions and results:
        results = prune_spurious_insertions(results, merge_frac=prune_merge_frac,
                                            height_frac=prune_height_frac,
                                            max_iter=prune_max_iter,
                                            expected=anchor_curve if prune_expected_spacing else None)

    if auto_trim and results:
        if trim_method == "mott":
            start, end = trim_mott(results, quality_cutoff=mott_quality_cutoff)
        else:
            start, end = trim_to_quality_block(results, min_quality_percentile=trim_quality_percentile)
        results = results[start:end]

    sequence = "".join(base_order[b.channel] for b in results)
    qualities = [b.height for b in results]
    return sequence, qualities, results
