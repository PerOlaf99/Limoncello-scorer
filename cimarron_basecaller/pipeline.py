"""
Top-level base-calling pipeline: preprocess -> per-segment (iterative
blind deconvolution -> Monte Carlo alignment -> putative peak detection
-> OmitOkN -> GapCheck) -> reassembly -> BaseQual scoring -> best
contiguous high-quality block selection.

This mirrors the overall flow of EP0944739A1 Figure 1. It is an
independent, from-scratch implementation guided by the patent's public
disclosure (application withdrawn; no patent was ever granted on this
method), not a port of the original compiled DLLs or their C++ source.
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np

from . import preprocessing, deconvolution, alignment, peaks, quality

BASES = "ACGT"  # caller supplies the channel->base mapping; default order


@dataclass
class BaseCallResult:
    sequence: str
    qualities: list[float]
    left_cutoff: int
    right_cutoff: int
    bands: list[peaks.Band]


def call_segment(segment: np.ndarray, prev_fbw: float | None = None):
    """Run deconvolution -> alignment -> peak detection -> fuzzy filtering
    on a single 2048-sample segment. Returns (bands, final_fbw)."""
    deconvolved, fbw, _spacing = deconvolution.iterative_deconvolve(segment, prev_fbw)

    best_shift = alignment.monte_carlo_align(deconvolved)
    aligned = alignment.apply_alignment(deconvolved, best_shift)

    band_list = peaks.detect_putative_peaks(aligned)
    peaks.measure_spacing_and_crossbanding(aligned, band_list)

    labels = peaks.omit_ok_n(band_list)
    band_list = [b for b, lbl in zip(band_list, labels) if lbl != "OMIT"]
    peaks.measure_spacing_and_crossbanding(aligned, band_list)

    split_flags = peaks.gap_check(band_list)
    band_list = peaks.apply_gap_splits(band_list, split_flags)
    peaks.measure_spacing_and_crossbanding(aligned, band_list)

    # Second OmitOkN pass to clean up any spurious insertions introduced
    # by GapCheck, per the patent's described two-pass refinement.
    labels = peaks.omit_ok_n(band_list)
    band_list = [b for b, lbl in zip(band_list, labels) if lbl != "OMIT"]

    return band_list, fbw


def basecall(trace: np.ndarray, base_order: str = BASES) -> BaseCallResult:
    """Run the full pipeline on a raw (n_samples, 4) trace and return a
    called sequence with per-base quality values and suggested cutoffs.
    """
    begin, end = preprocessing.detect_begin_end(trace)
    trimmed = trace[begin:end]

    baseline_subtracted = preprocessing.subtract_baseline(trimmed)

    all_bands: list[peaks.Band] = []
    fbw = None
    offset = 0
    for seg_start, seg_end, segment in deconvolution.iter_segments(baseline_subtracted):
        band_list, fbw = call_segment(segment, prev_fbw=fbw)
        for b in band_list:
            b.position += seg_start
        # Drop bands in the overlap region already covered by the previous
        # segment, to avoid double-calling.
        band_list = [b for b in band_list if b.position >= offset]
        all_bands.extend(band_list)
        offset = seg_end - deconvolution.SEGMENT_OVERLAP

    all_bands.sort(key=lambda b: b.position)

    quality.score_bands(all_bands)

    left_cut, right_cut = select_high_quality_block(all_bands)

    called = all_bands[left_cut:right_cut]
    sequence = "".join(base_order[b.channel] for b in called)
    qualities = [b.quality if b.quality is not None else 0.0 for b in called]

    return BaseCallResult(
        sequence=sequence,
        qualities=qualities,
        left_cutoff=left_cut,
        right_cutoff=right_cut,
        bands=all_bands,
    )


def select_high_quality_block(
    bands: list[peaks.Band],
    filter_widths: tuple[int, ...] = (1, 3, 5, 7, 9, 11),
    thresholds: tuple[float, ...] = (30, 40, 50, 60, 70, 80, 85, 90, 95),
) -> tuple[int, int]:
    """Find the longest contiguous block of above-threshold, moving-
    average-filtered quality values, scanning several (filter width,
    threshold) pairs and favoring narrow-filter/high-threshold results,
    per the patent's described selection surface (Section IV.D).
    """
    if not bands:
        return 0, 0

    qualities = np.array([b.quality if b.quality is not None else 0.0 for b in bands])

    best_len = -1
    best_range = (0, len(bands))
    best_score = -np.inf

    for w in filter_widths:
        kernel = np.ones(w) / w
        filtered = np.convolve(qualities, kernel, mode="same")
        for t in thresholds:
            above = filtered >= t
            start = None
            for i, ok in enumerate(list(above) + [False]):
                if ok and start is None:
                    start = i
                elif not ok and start is not None:
                    length = i - start
                    # Favor narrow filter / high threshold, per the patent.
                    score = length * (1.0 / w) * (t / 100.0)
                    if score > best_score:
                        best_score = score
                        best_len = length
                        best_range = (start, i)
                    start = None

    if best_len <= 0:
        return 0, len(bands)
    return best_range
