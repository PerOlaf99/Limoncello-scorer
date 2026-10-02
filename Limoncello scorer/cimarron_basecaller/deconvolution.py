"""
Iterative blind deconvolution with adaptive filter bandwidth (FBW).

The patent describes deblurring each 2048-sample segment against an
"unknown Laurentian blurring function", re-estimating the filter
bandwidth between an initial narrow-band guess and a refined value based
on measured band spacing, per:

    K = 2 * sqrt(ln(0.23) / -0.5)
    FBW = 0.5 + (K * median_band_spacing) * (2048 / (2*pi))

That formula is given explicitly in the patent and is used verbatim here
(it's a numeric relationship/fact, not a copyrightable expression). The
underlying deconvolution kernel itself ("Laurentian" = Lorentzian-shaped
blur, per the referenced Stockham & Ives blind-deconvolution patent
5,273,632) is implemented here as a standard frequency-domain Wiener-style
deconvolution against a Lorentzian PSF of the given FBW -- a reasonable,
independent implementation of "blind deconvolution with a Lorentzian
kernel", since the patent does not specify implementation-level details
beyond naming the technique and giving the FBW-selection formula.
"""

from __future__ import annotations
import numpy as np

SEGMENT_SIZE = 2048
SEGMENT_OVERLAP = 148


def fbw_from_band_spacing(median_band_spacing: float) -> float:
    """Map an estimated median band spacing (in samples) to a filter
    bandwidth value.

    CAVEAT: the patent's OCR'd text gives "FBW = 0.5 + (K * spacing) *
    (2048/(2*pi))" -- but the patent's own prose separately states "Band
    spacing and filter band width are inversely related." Those two
    statements are inconsistent (the OCR'd formula only multiplies), which
    strongly suggests the OCR flattened a fraction (division) into a
    multiplication when the source PDF's layout was parsed as plain text.
    This implementation follows the *inverse* relationship stated in the
    prose (dividing by spacing) since that's the more clearly and directly
    stated fact, but you may want to check the original PDF's rendered
    equation directly if exact numeric fidelity to the patent matters --
    Google Patents' OCR is not fully reliable for embedded equations.
    """
    k = 2.0 * np.sqrt(np.log(0.23) / -0.5)
    spacing = max(median_band_spacing, 1e-6)
    fbw = 0.5 + k * (SEGMENT_SIZE / (2.0 * np.pi)) / spacing
    # Clamp to a sane range relative to typical segment/band sizes.
    return float(np.clip(fbw, 1.0, SEGMENT_SIZE / 4))


def _lorentzian_psf(n: int, fbw: float) -> np.ndarray:
    """Lorentzian point-spread function of width `fbw`, centered, length n."""
    x = np.arange(n) - n // 2
    gamma = max(fbw, 1e-3)
    psf = 1.0 / (1.0 + (x / gamma) ** 2)
    psf /= psf.sum()
    return psf


def blind_deconvolve_segment(
    segment: np.ndarray, fbw: float, noise_reg: float = 1e-2
) -> np.ndarray:
    """Deconvolve one (n, 4) segment against a Lorentzian PSF of the given
    filter bandwidth, channel by channel, using regularized (Wiener-style)
    frequency-domain division. Also min-max normalizes each channel's
    amplitude into [0, 1] per the patent's stated normalization step.
    """
    n = segment.shape[0]
    psf = _lorentzian_psf(n, fbw)
    psf_fft = np.fft.rfft(np.fft.ifftshift(psf))

    out = np.empty_like(segment, dtype=float)
    for c in range(segment.shape[1]):
        sig_fft = np.fft.rfft(segment[:, c])
        wiener = np.conj(psf_fft) / (np.abs(psf_fft) ** 2 + noise_reg)
        deconvolved = np.fft.irfft(sig_fft * wiener, n=n)
        deconvolved = np.clip(deconvolved, 0, None)
        rng = deconvolved.max() - deconvolved.min()
        if rng > 1e-9:
            deconvolved = (deconvolved - deconvolved.min()) / rng
        out[:, c] = deconvolved
    return out


def estimate_median_band_spacing(segment: np.ndarray) -> float:
    """Rough initial band-spacing estimate from envelope peak spacing,
    used to bootstrap the first FBW guess for a segment."""
    envelope = segment.max(axis=1)
    peaks = []
    for i in range(1, len(envelope) - 1):
        if envelope[i] > envelope[i - 1] and envelope[i] >= envelope[i + 1]:
            peaks.append(i)
    if len(peaks) < 2:
        return 20.0  # fallback default spacing guess
    spacings = np.diff(peaks)
    return float(np.median(spacings))


def iterative_deconvolve(segment: np.ndarray, initial_fbw: float | None = None):
    """Two-pass deconvolution of a segment: an initial narrow-band pass to
    get a conservative band-density estimate, then a refined pass with the
    FBW re-derived from that estimate -- matching the patent's stated
    "read twice" procedure per segment.

    Returns (deconvolved_segment, final_fbw, median_band_spacing).
    """
    if initial_fbw is None:
        # Narrow-band guess: assume tight spacing so the first pass is
        # conservative and doesn't overestimate band density.
        initial_fbw = fbw_from_band_spacing(6.0)

    first_pass = blind_deconvolve_segment(segment, initial_fbw)
    spacing_estimate = estimate_median_band_spacing(first_pass)
    refined_fbw = fbw_from_band_spacing(spacing_estimate)
    second_pass = blind_deconvolve_segment(segment, refined_fbw)

    return second_pass, refined_fbw, spacing_estimate


def iter_segments(trace: np.ndarray, size: int = SEGMENT_SIZE, overlap: int = SEGMENT_OVERLAP):
    """Yield (start, end, segment) tuples covering `trace` in overlapping
    windows, per the patent's segment-reading scheme. Short final segments
    are padded with pseudo-random noise up to `size`, as the patent
    specifies (to avoid non-random artifacts from zero-padding)."""
    n = trace.shape[0]
    start = 0
    rng = np.random.default_rng(0)
    while start < n:
        end = min(start + size, n)
        seg = trace[start:end]
        if seg.shape[0] < size:
            pad_len = size - seg.shape[0]
            noise = rng.random((pad_len, trace.shape[1])) * seg.max(initial=1e-6) * 0.01
            seg = np.vstack([seg, noise])
        yield start, end, seg
        if end >= n:
            break
        start = end - overlap
