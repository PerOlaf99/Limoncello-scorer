"""
Monte Carlo lane-alignment search.

Per EP0944739A1 (Fig. 14): rather than aligning lanes using band position
information, search a 3D space of shift offsets between channel pairs
(e.g. A-vs-G, {A,G}-vs-C, {A,G,C}-vs-T) to MAXIMIZE the integral of the
envelope (max across the 4 shifted channels) at each sample. Good
alignment => bands sit "shoulder to shoulder" => a taller, larger-area
envelope. Uses a lattice-based Monte Carlo refinement (attributed in the
patent to W.L. Price, The Computer Journal Vol. 20 No. 4): keep a
population of candidate shift-triples, repeatedly replace the worst
candidate with a random perturbation of the best, until convergence.

This is an independent implementation of that search strategy.
"""

from __future__ import annotations
import numpy as np


def envelope_integral(segment: np.ndarray, shifts: tuple[int, int, int, int]) -> float:
    """Sum of the per-sample max across 4 channels, each shifted by the
    given integer offsets (one of which is conventionally 0)."""
    n = segment.shape[0]
    shifted = []
    for c, s in enumerate(shifts):
        if s == 0:
            shifted.append(segment[:, c])
        elif s > 0:
            col = np.zeros(n)
            col[s:] = segment[:n - s, c]
            shifted.append(col)
        else:
            col = np.zeros(n)
            col[:n + s] = segment[-s:, c]
            shifted.append(col)
    envelope = np.max(np.stack(shifted, axis=1), axis=1)
    return float(envelope.sum())


def _triple_to_shifts(triple: tuple[int, int, int]) -> tuple[int, int, int, int]:
    """Convert an (x, y, z) alignment triple into 4 per-channel shifts,
    with channel 0 (A) held fixed at 0, matching the patent's axis
    definitions (x: A-G, y: AG-C, z: AGC-T)."""
    x, y, z = triple
    return (0, x, x + y, x + y + z)


def monte_carlo_align(
    segment: np.ndarray,
    center: tuple[int, int, int] = (0, 0, 0),
    n_iterations: int = 200,
    lattice_radius: int = 6,
    rng: np.random.Generator | None = None,
) -> tuple[int, int, int]:
    """Search for the (x, y, z) shift triple that maximizes envelope
    integral, via lattice-seeded Monte Carlo refinement.

    Returns the best (x, y, z) triple found.
    """
    if rng is None:
        rng = np.random.default_rng(0)

    # Seed initial population from lattice points of concentric cubes
    # around `center`, as described.
    candidates = []
    for r in range(1, lattice_radius + 1):
        for dx in (-r, 0, r):
            for dy in (-r, 0, r):
                for dz in (-r, 0, r):
                    if dx == dy == dz == 0:
                        continue
                    candidates.append((
                        center[0] + dx, center[1] + dy, center[2] + dz
                    ))
    if not candidates:
        candidates = [center]

    scored = [
        (envelope_integral(segment, _triple_to_shifts(t)), t) for t in candidates
    ]
    scored.sort(key=lambda p: p[0])  # ascending: worst first

    best_score, best_triple = scored[-1]

    for _ in range(n_iterations):
        # Replace worst with a random perturbation of best.
        perturb = tuple(
            int(best_triple[i] + rng.integers(-1, 2)) for i in range(3)
        )
        new_score = envelope_integral(segment, _triple_to_shifts(perturb))
        scored[0] = (new_score, perturb)
        scored.sort(key=lambda p: p[0])
        if scored[-1][0] > best_score:
            best_score, best_triple = scored[-1]

    return best_triple


def apply_alignment(segment: np.ndarray, triple: tuple[int, int, int]) -> np.ndarray:
    """Apply the found shift triple to a segment, returning aligned traces
    (channels shifted and zero-padded at their leading/trailing edge)."""
    shifts = _triple_to_shifts(triple)
    n = segment.shape[0]
    out = np.zeros_like(segment, dtype=float)
    for c, s in enumerate(shifts):
        if s == 0:
            out[:, c] = segment[:, c]
        elif s > 0:
            out[s:, c] = segment[:n - s, c]
        else:
            out[:n + s, c] = segment[-s:, c]
    return out
