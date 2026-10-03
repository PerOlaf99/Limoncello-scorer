#!/usr/bin/env python3
"""Regular fragment-length sizing for capillary-electrophoresis traces.

This is the *sizing* half of genotyping, deliberately separate from the
SNP/allele caller in :mod:`scorer`.  A sample is amplified (PCR) so every
fragment is flanked by the same primer pair; one primer carries a fluorophore.
A size standard ("ladder") -- a set of fragments of known length, labelled with
a *different* fluorophore -- is run in the same lane.  Sizing means:

1. find the ladder's peaks in its channel,
2. match them, in order, to the kit's known fragment lengths,
3. fit a migration curve ``scan -> bp`` through those anchors, and
4. read every sample peak back through the curve to get its length in bp.

Nothing here is MegaBACE-specific except the channel convention: ``acgt`` is
always A,C,G,T while the physical channels follow the plate's dye order
("ACTG" on a MegaBACE), so a channel is always resolved through
:func:`analyzer_core.acgt_index_for_channel`.  The module is headless
(numpy/scipy only) so :mod:`scorer` can use it without Tk or matplotlib.

Ladders
-------
A ladder is just a name, a fluorophore and an ascending list of lengths.  Two
well-known ROX standards are bundled (see :data:`BUILTIN_LADDERS`); any other
kit is accepted either as a JSON file or as a comma-separated list on the
command line, so a lab can use whatever ladder it actually has -- including a
home-made one with a fluorophore of its choice.  MegaBACE-compatible standards
whose fragment tables are *not* bundled are listed in
:data:`KNOWN_MEGABACE_STANDARDS` for reference.

CLI::

    python scorer.py ladders
    python scorer.py size well.rsd --ladder geneflo1000_rox \
        --ladder-channel 4 --sample-channel 2 --out sizes.csv
    python scorer.py size well.rsd --lengths 35,50,75,100,139,150,200,250
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from analyzer_core import acgt_index_for_channel

# numpy dropped ``trapz`` in 2.0 (renamed ``trapezoid``); accept both.
_trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz")

# Physical channel order on a MegaBACE plate.  A channel number is meaningless
# without it: Ch2 is C here, but would be another base on a different dye set.
DEFAULT_BASE_ORDER = "ACTG"

# Sizing tolerances.  A sample peak whose length lands outside the ladder's
# span is reported but flagged: the curve extrapolates there and its accuracy
# degrades quickly past the outermost anchor.
DEFAULT_LADDER_PROMINENCE_FRAC = 0.015
DEFAULT_LADDER_HEIGHT_FRAC = 0.02
DEFAULT_SAMPLE_PROMINENCE_FRAC = 0.02
DEFAULT_SAMPLE_HEIGHT_FRAC = 0.03
# A ladder needs at least this many matched anchors before the sizing curve
# (and its leave-one-out accuracy estimate) mean anything.
MIN_LADDER_ANCHORS = 3


@dataclass(frozen=True)
class Ladder:
    """A size standard: a name, the dye it is labelled with, and its lengths."""

    name: str
    lengths: Tuple[float, ...]
    dye: str = "ROX"
    source: str = ""
    notes: str = ""

    def __post_init__(self):
        if not self.lengths:
            raise ValueError(f"Ladder {self.name!r} has no fragment lengths")
        if any(b <= 0 for b in self.lengths):
            raise ValueError(f"Ladder {self.name!r} lengths must be positive")
        if list(self.lengths) != sorted(self.lengths):
            object.__setattr__(self, "lengths", tuple(sorted(self.lengths)))

    @property
    def span(self) -> Tuple[float, float]:
        return self.lengths[0], self.lengths[-1]


@dataclass(frozen=True)
class Peak:
    scan: int
    height: float
    prominence: float


@dataclass
class SizingResult:
    ladder: Ladder
    ladder_channel: int
    sample_channel: int
    base_order: str
    anchors: List[Tuple[int, float]]
    quality: Dict[str, object]
    rows: List[dict]
    well: str = ""
    file: str = ""
    warnings: List[str] = field(default_factory=list)

    def summary(self) -> dict:
        return {
            "well": self.well,
            "file": self.file,
            "ladder": self.ladder.name,
            "ladder_dye": self.ladder.dye,
            "ladder_channel": self.ladder_channel,
            "sample_channel": self.sample_channel,
            "base_order": self.base_order,
            "n_anchors": len(self.anchors),
            "n_sample_peaks": len(self.rows),
            "in_range_peaks": sum(1 for r in self.rows if r["in_range"]),
            **{f"ladder_{k}": v for k, v in self.quality.items()},
        }


# --------------------------------------------------------------------------- #
# built-in ladders
# --------------------------------------------------------------------------- #
# Only standards whose fragment table can be sourced are bundled.  Everything
# else is supplied by the user (``--lengths`` or a JSON file); see
# KNOWN_MEGABACE_STANDARDS for the MegaBACE-compatible kits and why their
# tables are not hard-coded here.
BUILTIN_LADDERS: Dict[str, dict] = {
    "genescan500_rox": {
        "name": "GeneScan 500 ROX",
        "dye": "ROX",
        "lengths": [35, 50, 75, 100, 139, 150, 160, 200, 250, 300,
                    340, 350, 400, 450, 490, 500],
        "source": "Applied Biosystems GeneScan 500 ROX Size Standard",
        "notes": "16 fragments; 100-mers (100,200,300,400,500) are spiked.",
    },
    "geneflo1000_rox": {
        "name": "Geneflo 1000 ROX",
        "dye": "ROX",
        "lengths": list(range(400, 1001, 25)),
        "source": "ChimerX Geneflo 1000 Size Standard (cat. 3127), insert",
        "notes": "25 equally spaced bands, 400-1000 bp; 100-mers spiked.",
    },
}

# MegaBACE 1000-compatible size standards.  These are the kits a MegaBACE user
# is most likely to have; their exact fragment tables are printed on the kit
# inserts and are not reproduced from a primary source here (a half-remembered
# ladder is worse than none -- it silently mis-sizes every peak).  Supply the
# lengths from your own insert with ``--lengths`` or a ladder JSON file.
# Resolution figures come from the MegaBACE Long Read Matrix specification.
KNOWN_MEGABACE_STANDARDS = [
    {"key": "et400r", "name": "MegaBACE ET400-R", "dye": "ET-ROX",
     "range_bp": "~20-400", "resolution_bp": 0.3,
     "catalog": "GE/Cytiva MegaBACE ET400-R Size Standard"},
    {"key": "et550r", "name": "MegaBACE ET550-R", "dye": "ET-ROX",
     "range_bp": "~20-550", "resolution_bp": 0.5,
     "catalog": "GE/Cytiva MegaBACE ET550-R Size Standard"},
    {"key": "et900r", "name": "MegaBACE ET900-R", "dye": "ET-ROX",
     "range_bp": "~20-900", "resolution_bp": 2.0,
     "catalog": "GE/Cytiva MegaBACE ET900-R, cat. 25-6900-01"},
    {"key": "geneflo1000_rox", "name": "Geneflo 1000 ROX (bundled)",
     "dye": "ROX", "range_bp": "400-1000", "resolution_bp": None,
     "catalog": "ChimerX 3127 (usable as a ROX-channel standard)"},
]


def list_ladders() -> List[str]:
    return sorted(BUILTIN_LADDERS)


def _coerce_lengths(lengths: Sequence) -> Tuple[float, ...]:
    out = []
    for value in lengths:
        if isinstance(value, str):
            value = value.strip()
            if not value:
                continue
        out.append(float(value))
    return tuple(out)


def ladder_from_lengths(lengths: Sequence, name: str = "custom",
                        dye: str = "ROX", source: str = "",
                        notes: str = "") -> Ladder:
    return Ladder(name=name, lengths=_coerce_lengths(lengths), dye=dye,
                  source=source, notes=notes)


def parse_lengths(spec: str) -> Ladder:
    """A comma/space separated list of fragment lengths -> a custom ladder."""
    if not spec or not spec.strip():
        raise ValueError("Empty ladder length list")
    parts = [p for p in spec.replace(";", ",").replace(" ", ",").split(",") if p]
    return ladder_from_lengths(parts, name="custom", dye="custom")


def load_ladder(spec) -> Ladder:
    """Resolve *spec* to a :class:`Ladder`.

    *spec* may be a built-in name (case-insensitive), a comma-separated list of
    lengths, or a path to a JSON file::

        {"name": "My ROX ladder", "dye": "ROX",
         "lengths": [50, 100, 150, 200], "notes": "..."}
    """
    if isinstance(spec, Ladder):
        return spec
    if not isinstance(spec, str) or not spec.strip():
        raise ValueError("No ladder given")

    text = spec.strip()
    path = Path(text)
    if path.is_file():
        data = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(data, list):
            data = {"lengths": data}
        return Ladder(
            name=str(data.get("name") or path.stem),
            lengths=_coerce_lengths(data.get("lengths", [])),
            dye=str(data.get("dye") or "custom"),
            source=str(data.get("source") or str(path)),
            notes=str(data.get("notes") or ""),
        )

    key = text.lower()
    if key in BUILTIN_LADDERS:
        data = BUILTIN_LADDERS[key]
        return Ladder(name=data["name"], lengths=_coerce_lengths(data["lengths"]),
                      dye=data.get("dye", "ROX"),
                      source=data.get("source", ""), notes=data.get("notes", ""))

    if any(ch.isdigit() for ch in text):
        return parse_lengths(text)

    raise ValueError(
        f"Unknown ladder {spec!r}. Known: {', '.join(list_ladders())} "
        "(or pass lengths / a JSON path)")


# --------------------------------------------------------------------------- #
# peak detection (numpy/scipy only -- mirrors genotyping's picker, not imports)
# --------------------------------------------------------------------------- #
def detect_peaks(trace, *, min_prominence_frac: float = 0.02,
                 min_height_frac: float = 0.03, cut: int = 0,
                 smooth: bool = True) -> List[Peak]:
    """Prominent local maxima of one channel, above a robust local baseline.

    The baseline is the median of the first ~400 scans (the injection region),
    which is where the trace is flattest; prominence/height thresholds are
    fractions of the largest baseline-corrected sample so they scale with the
    run's gain.
    """
    from scipy.signal import find_peaks

    y = np.asarray(trace, dtype=float).ravel()
    n = y.size
    if n == 0:
        return []
    if smooth and n >= 9:
        from scipy.signal import medfilt, savgol_filter
        y = savgol_filter(medfilt(y, 5), 9, 2, mode="interp")
    base = float(np.median(y[:min(400, n)]))
    z = y - base
    mx = float(z.max())
    if not np.isfinite(mx) or mx <= 0:
        return []
    idx, props = find_peaks(z, prominence=mx * min_prominence_frac,
                            height=mx * min_height_frac)
    peaks = []
    for k, i in enumerate(idx):
        if i < cut:
            continue
        peaks.append(Peak(int(i), float(z[i]), float(props["prominences"][k])))
    return peaks


def _peak_areas(trace, peaks: Sequence[Peak]) -> List[float]:
    """Trapezoid area of each peak down to the local valley between neighbours.

    The integration window runs from the midpoint toward each neighbour to the
    next, and the baseline is the floor of that window -- so two partially
    overlapping peaks split their shared valley instead of double-counting it.
    """
    y = np.asarray(trace, dtype=float).ravel()
    scans = [p.scan for p in peaks]
    areas = []
    for i, scan in enumerate(scans):
        left = int((scans[i - 1] + scan) / 2) if i > 0 else max(0, scan - 20)
        right = (int((scan + scans[i + 1]) / 2) if i + 1 < len(scans)
                 else min(len(y), scan + 20))
        seg = y[left:max(right, left + 1)]
        if seg.size < 2:
            areas.append(0.0)
            continue
        floor = float(np.min(seg))
        areas.append(max(0.0, float(_trapz(seg - floor, dx=1.0))))
    return areas


# --------------------------------------------------------------------------- #
# ladder <-> length matching and the sizing curve
# --------------------------------------------------------------------------- #
def _log_linear_fit(u: np.ndarray, v: np.ndarray) -> Tuple[float, float]:
    """Least-squares ``v = a + b*u`` with a guaranteed positive slope.

    In capillary electrophoresis scan number grows roughly with log(bp), so
    this is a better initial guess than a raw linear fit; ``b`` is forced
    positive because longer fragments must migrate later, never earlier.
    """
    if u.size < 2 or float(np.ptp(u)) == 0.0:
        return 0.0, 1.0
    b, a = np.polyfit(u, v, 1)
    if not np.isfinite(b) or b <= 0:
        b = max((v[-1] - v[0]) / (u[-1] - u[0]), 1e-9)
        a = float(v[0] - b * u[0])
    return float(a), float(b)


def _align_dp(scans: Sequence[int], lengths: Sequence[float],
              pred: np.ndarray) -> List[Tuple[int, int]]:
    """Monotonic alignment of detected peaks to expected lengths.

    A Needleman-Wunsch alignment: match, skip a detected peak, or skip an
    expected length.  The match cost is the scan distance to that length's
    predicted position, and every gap costs one average peak spacing, so a
    real peak is always cheaper to match than to skip.
    """
    n, k = len(scans), len(lengths)
    if n == 0 or k == 0:
        return []
    span = max(float(scans[-1]) - float(scans[0]), 1.0)
    gap = span / max(n, k)
    inf = float("inf")
    dp = np.full((n + 1, k + 1), inf)
    back = np.zeros((n + 1, k + 1), dtype=np.int8)  # 0 match, 1 skip scan, 2 skip length
    dp[0, 0] = 0.0
    for i in range(n + 1):
        for j in range(k + 1):
            if i == 0 and j == 0:
                continue
            best, choice = inf, 0
            if i and j:
                cost = dp[i - 1, j - 1] + abs(pred[j - 1] - scans[i - 1])
                if cost < best:
                    best, choice = cost, 0
            if i and dp[i - 1, j] + gap < best:
                best, choice = dp[i - 1, j] + gap, 1
            if j and dp[i, j - 1] + gap < best:
                best, choice = dp[i, j - 1] + gap, 2
            dp[i, j], back[i, j] = best, choice
    pairs = []
    i, j = n, k
    while i or j:
        choice = back[i, j]
        if choice == 0 and i and j:
            pairs.append((i - 1, j - 1))
            i, j = i - 1, j - 1
        elif choice == 1:
            i -= 1
        else:
            j -= 1
    pairs.reverse()
    return pairs


def align_peaks(scans: Sequence[int], lengths: Sequence[float],
                iterations: int = 3) -> List[Tuple[int, float]]:
    """Match ladder peak scans to the kit's known lengths.

    Bootstraps a log-linear migration model from the endpoints (or a direct
    pairing when the counts already agree), aligns with DP, refits the model on
    the matches, and repeats -- the iteration is what absorbs the curvature a
    single global fit cannot.  Returns ``(scan, length)`` anchors.
    """
    scans = [int(s) for s in scans]
    lengths = [float(b) for b in lengths]
    if not scans or not lengths:
        return []
    if len(scans) == len(lengths):
        a, b = _log_linear_fit(np.log(lengths), np.asarray(scans, float))
    else:
        a, b = _log_linear_fit(np.log([lengths[0], lengths[-1]]),
                               np.asarray([scans[0], scans[-1]], float))
    best: List[Tuple[int, int]] = []
    for _ in range(max(1, iterations)):
        pred = a + b * np.log(lengths)
        pairs = _align_dp(scans, lengths, pred)
        if len(pairs) < 2:
            best = pairs or best
            break
        u = np.log([lengths[j] for _, j in pairs])
        v = np.asarray([scans[i] for i, _ in pairs], float)
        a, b = _log_linear_fit(u, v)
        best = pairs
    return [(scans[i], lengths[j]) for i, j in best]


class SizeCurve:
    """Piecewise-cubic ``scan -> bp`` curve through the ladder anchors.

    PCHIP is used because it is shape-preserving: the curve stays monotone
    between anchors and never overshoots into a non-physical length.  Outside
    the anchor span it extrapolates linearly from the outermost segment, and
    callers are expected to flag those points (they are the least trustworthy).
    """

    def __init__(self, scans: Sequence[int], lengths: Sequence[float]):
        from scipy.interpolate import PchipInterpolator

        pairs = sorted(zip((int(s) for s in scans), (float(b) for b in lengths)))
        xs, ys = [], []
        for x, y in pairs:                       # dedupe scans, keep first
            if xs and x <= xs[-1]:
                continue
            xs.append(x)
            ys.append(y)
        if len(xs) < 2:
            raise ValueError("A sizing curve needs at least two ladder anchors")
        self.scan_min, self.scan_max = float(xs[0]), float(xs[-1])
        self._pchip = PchipInterpolator(xs, ys, extrapolate=False)
        self._xs = np.asarray(xs, float)
        self._ys = np.asarray(ys, float)

    def bp(self, scan) -> np.ndarray:
        scan = np.atleast_1d(np.asarray(scan, dtype=float))
        out = np.empty_like(scan)
        inside = (scan >= self.scan_min) & (scan <= self.scan_max)
        if inside.any():
            out[inside] = self._pchip(scan[inside])
        if (~inside).any():
            out[~inside] = self._extrapolate(scan[~inside])
        return out

    def _extrapolate(self, scan: np.ndarray) -> np.ndarray:
        lo = scan < self.scan_min
        hi = ~lo
        out = np.empty_like(scan)
        if lo.any():
            slope = ((self._ys[1] - self._ys[0]) /
                     (self._xs[1] - self._xs[0]))
            out[lo] = self._ys[0] + slope * (scan[lo] - self._xs[0])
        if hi.any():
            slope = ((self._ys[-1] - self._ys[-2]) /
                     (self._xs[-1] - self._xs[-2]))
            out[hi] = self._ys[-1] + slope * (scan[hi] - self._xs[-1])
        return out

    @property
    def in_range(self):
        def _check(scan) -> np.ndarray:
            s = np.atleast_1d(np.asarray(scan, dtype=float))
            return (s >= self.scan_min) & (s <= self.scan_max)
        return _check


def ladder_fit_quality(scans: Sequence[int],
                       lengths: Sequence[float]) -> Dict[str, object]:
    """Leave-one-out accuracy of the sizing curve, in bp.

    A residual on the fit itself is identically zero (PCHIP passes through its
    anchors), so the honest measure is how well each anchor is predicted from
    the others.  Needs at least four anchors to say anything; below that the
    errors are reported as ``None`` rather than as a falsely precise zero.
    """
    scans = [int(s) for s in scans]
    lengths = [float(b) for b in lengths]
    n = len(scans)
    quality: Dict[str, object] = {
        "anchors": n,
        "rms_error_bp": None,
        "max_error_bp": None,
        "rms_error_rel": None,
    }
    if n < 4:
        return quality
    errors = []
    for hold in range(n):
        xs = [scans[i] for i in range(n) if i != hold]
        ys = [lengths[i] for i in range(n) if i != hold]
        try:
            curve = SizeCurve(xs, ys)
        except ValueError:
            continue
        errors.append(float(curve.bp([scans[hold]])[0] - lengths[hold]))
    if errors:
        err = np.abs(errors)
        denom = np.abs(lengths)
        quality["rms_error_bp"] = float(np.sqrt(np.mean(np.square(errors))))
        quality["max_error_bp"] = float(err.max())
        quality["rms_error_rel"] = float(np.sqrt(np.mean((err / denom) ** 2)))
    return quality


# --------------------------------------------------------------------------- #
# top-level sizing
# --------------------------------------------------------------------------- #
def size_trace(doc, ladder, *, ladder_channel: int = 4,
               sample_channel: int = 2, base_order: Optional[str] = None,
               ladder_prominence_frac: float = DEFAULT_LADDER_PROMINENCE_FRAC,
               ladder_height_frac: float = DEFAULT_LADDER_HEIGHT_FRAC,
               sample_prominence_frac: float = DEFAULT_SAMPLE_PROMINENCE_FRAC,
               sample_height_frac: float = DEFAULT_SAMPLE_HEIGHT_FRAC) -> SizingResult:
    """Size every sample peak of one loaded trace against *ladder*.

    ``ladder_channel`` / ``sample_channel`` are 1..4 *physical* channels in the
    plate's dye order (see :func:`analyzer_core.acgt_index_for_channel`); the
    default (ladder on Ch4, sample on Ch2) is the common MegaBACE ROX-ladder /
    FAM-sample layout but is a property of the kit, so callers should say.
    """
    if isinstance(ladder, str):
        ladder = load_ladder(ladder)
    order = (base_order or getattr(doc, "base_order", None) or DEFAULT_BASE_ORDER)
    warnings: List[str] = []

    acgt = np.asarray(doc.acgt, dtype=float)
    lad_idx = acgt_index_for_channel(order, ladder_channel)
    sam_idx = acgt_index_for_channel(order, sample_channel)
    if ladder_channel == sample_channel:
        warnings.append("ladder and sample channel are the same; sizing will "
                        "match the sample's own peaks to the ladder")

    ladder_trace = acgt[:, lad_idx]
    sample_trace = acgt[:, sam_idx]

    ladder_peaks = detect_peaks(ladder_trace,
                                min_prominence_frac=ladder_prominence_frac,
                                min_height_frac=ladder_height_frac)
    anchors = align_peaks([p.scan for p in ladder_peaks], ladder.lengths)
    quality = ladder_fit_quality([s for s, _ in anchors],
                                 [b for _, b in anchors])

    if len(anchors) < MIN_LADDER_ANCHORS:
        warnings.append(
            f"only {len(anchors)} ladder peak(s) matched the {ladder.name} "
            f"table ({len(ladder_peaks)} peak(s) found); results are not "
            "trustworthy")
    if ladder_peaks and anchors:
        matched = len(anchors) / float(len(ladder_peaks))
        if matched < 0.6:
            warnings.append(
                f"only {matched:.0%} of ladder peaks matched the table; check "
                "the ladder channel and the fragment lengths")

    curve = None
    if len(anchors) >= 2:
        curve = SizeCurve([s for s, _ in anchors], [b for _, b in anchors])
    elif anchors:
        warnings.append("not enough anchors to build a sizing curve")

    sample_peaks = detect_peaks(sample_trace,
                                min_prominence_frac=sample_prominence_frac,
                                min_height_frac=sample_height_frac)
    areas = _peak_areas(sample_trace, sample_peaks)
    base_letter = {"A": "A", "C": "C", "G": "G", "T": "T"}[order[sample_channel - 1]]
    well = getattr(doc, "well", "") or Path(str(getattr(doc, "path", ""))).stem
    file_name = Path(str(getattr(doc, "path", ""))).name

    rows: List[dict] = []
    if curve is not None and sample_peaks:
        scans = [p.scan for p in sample_peaks]
        lengths_bp = curve.bp(scans)
        in_range = curve.in_range(scans)
        for peak, area, bp, rng in zip(sample_peaks, areas, lengths_bp, in_range):
            rows.append({
                "file": file_name,
                "well": well,
                "channel": sample_channel,
                "base": base_letter,
                "scan": peak.scan,
                "length_bp": round(float(bp), 2),
                "height": round(float(peak.height), 3),
                "area": round(float(area), 3),
                "in_range": bool(rng),
            })

    return SizingResult(ladder=ladder, ladder_channel=ladder_channel,
                        sample_channel=sample_channel, base_order=order,
                        anchors=anchors, quality=quality, rows=rows,
                        well=str(well), file=file_name, warnings=warnings)
