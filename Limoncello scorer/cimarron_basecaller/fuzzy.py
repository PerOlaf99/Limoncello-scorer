"""
Fuzzy set / fuzzy logic primitives.

Implements a piecewise-linear ("trapezoidal") fuzzy set with standard
AND (min), OR (max), NOT (1-x) operators and centroid defuzzification.
This mirrors the general fuzzy-logic approach described in EP0944739A1
(University of Utah Research Foundation / Cimarron base-caller patent,
application withdrawn, no patent ever granted) but is an independent
implementation, not a transcription of the original C++.

A fuzzy set here is just a piecewise-linear membership function defined
by a list of (x, y) breakpoints, with y in [0, 1]. membership(x)
interpolates linearly between breakpoints and clamps outside the range.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence
import numpy as np


@dataclass
class FuzzySet:
    """A piecewise-linear fuzzy membership function."""
    xs: np.ndarray
    ys: np.ndarray
    name: str = ""

    def __post_init__(self):
        self.xs = np.asarray(self.xs, dtype=float)
        self.ys = np.asarray(self.ys, dtype=float)
        if self.xs.shape != self.ys.shape:
            raise ValueError("xs and ys must have the same shape")
        if len(self.xs) < 2:
            raise ValueError("A fuzzy set needs at least two breakpoints")

    def membership(self, x: float) -> float:
        """Degree of membership (0..1) of value x in this fuzzy set."""
        if x <= self.xs[0]:
            return float(self.ys[0])
        if x >= self.xs[-1]:
            return float(self.ys[-1])
        return float(np.interp(x, self.xs, self.ys))

    def scaled(self, amplitude: float, n_points: int = 50) -> "FuzzySet":
        """Return a new set whose membership values are capped at `amplitude`.

        Used to combine a rule's firing strength with its output/conclusion
        fuzzy set before defuzzification.
        """
        xs = np.linspace(self.xs[0], self.xs[-1], n_points)
        ys = np.minimum([self.membership(x) for x in xs], amplitude)
        return FuzzySet(xs, ys, name=self.name)

    def union_with(self, other: "FuzzySet", n_points: int = 200) -> "FuzzySet":
        """Fuzzy OR (max) of this set and another, resampled onto a shared grid."""
        lo = min(self.xs[0], other.xs[0])
        hi = max(self.xs[-1], other.xs[-1])
        xs = np.linspace(lo, hi, n_points)
        ys = np.maximum(
            [self.membership(x) for x in xs],
            [other.membership(x) for x in xs],
        )
        return FuzzySet(xs, ys)

    def centroid(self) -> float:
        """Defuzzify by centroid (center of mass) of the membership curve."""
        area = np.trapezoid(self.ys, self.xs)
        if area <= 1e-12:
            return float((self.xs[0] + self.xs[-1]) / 2.0)
        moment = np.trapezoid(self.xs * self.ys, self.xs)
        return float(moment / area)


def AND(*vals: float) -> float:
    """Fuzzy AND = min of truth values."""
    return float(min(vals))


def OR(*vals: float) -> float:
    """Fuzzy OR = max of truth values."""
    return float(max(vals))


def NOT(val: float) -> float:
    """Fuzzy NOT = complement."""
    return float(1.0 - val)


def make_zero_set(lo: float = 0.0, hi: float = 1.0) -> FuzzySet:
    """A fuzzy set that is identically zero everywhere -- the starting
    'conclusion' accumulator before rules are OR'd into it."""
    return FuzzySet(np.array([lo, hi]), np.array([0.0, 0.0]))
