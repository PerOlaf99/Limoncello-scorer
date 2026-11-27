"""SNP genotyping by CTCE with a four-peak internal standard -- no marks needed.

Each well carries the same internal standard (STD) in its own dye: the SNP
fragment as a heterozygote, so four bands in a fixed order --

    Homo1  ->  Homo2  ->  (long gap)  ->  Hetero1  ->  Hetero2

The heteroduplexes may run so close that they merge into one band.  The
sample is the same fragment in another dye, so its species run where the
standard's do, a few scans ahead (the dyes' mobility differs).  A sample is

    hom1   Homo1 only
    hom2   Homo2 only
    het    both homoduplexes (any ratio: allelic imbalance, AI) and, when the
           heteroduplexes resolve, the heteroduplex band(s)

The engine works per plate, because the plate is what is shared:

1. **Channels.**  The sample is the plate's strongest dye; the standard is the
   other channel whose band pattern most wells agree on.
2. **The standard's spacing** (Homo1->Homo2 = g1, Homo2->Hetero1 = g2,
   Hetero1->Hetero2 = g3) by a plate vote: every well proposes the ordered
   band sets among its strongest standard-channel peaks; the spacing most
   wells reproduce wins.  A proposal only counts when its bands are roughly
   equal in height (one fragment, one amount) and the sample has a big peak
   right by Homo1 or Homo2 (the sample *is* that fragment) -- which rejects
   primer and dye fronts.  g2 is the long gap; g3 is 0 when the pair merges.
3. **Placement per well**, with a soft gap penalty: spacing varies capillary
   to capillary (on the T9 plate Homo1->Homo2 runs 57-100 scans).
4. **Sample offset** per plate: where the sample's main peak sits relative to
   the nearest standard homoduplex.
5. **Measure** the sample at the four positions (peak height above the
   local floor, in noise units), and **call**.

Validated against Per's hand calls on 12 plates (rs1695, ABCC2, CYBA, GSTA1,
IL10; 1,138 wells): 94 % agreement with every threshold chosen with the plate
itself left out.  See docs/GENOTYPING_ENGINE.md.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np
from scipy.ndimage import minimum_filter1d, uniform_filter1d
from scipy.signal import find_peaks

from analyzer_core import acgt_index_for_channel

# ---------------------------------------------------------------- settings
CUT = 1500                 # scans; nothing before the injection front counts
BASELINE_WINDOW = 151      # scans; rolling-minimum baseline
VOTE_PEAKS = 14            # strongest standard-channel peaks a well proposes from
PLACE_PEAKS = 14
VOTE_PROMINENCE = 6.0      # x noise, peaks a well may propose
PLACE_PROMINENCE = 4.0
G1_RANGE = (45, 200)       # scans, Homo1 -> Homo2 (marks: 49-140 over six genes)
G2_RANGE = (120, 330)      # scans, Homo2 -> Hetero1, the long gap (marks: 135-291)
G2_MIN_RATIO = 1.3         # the long gap is at least this x g1
EQUIMOLAR_MIN = 0.3        # weakest band / strongest band in a proposal
ANCHOR_FRAC = 0.5          # sample peak by a homoduplex vs the biggest one in the standard's stretch
ANCHOR_WINDOW = (30, 8)    # scans before / after the standard band
VOTE_TOL = (0.12, 0.12, 0.25)   # relative gap tolerance between agreeing wells
VOTE_TOL_MIN = (8, 8, 4)        # scans
FOUR_BAND_SHARE = 0.6      # a 4-band pattern wins with this share of 3-band support
FOUR_BAND_STRENGTH = 0.75  # ...and this share of its vote strength
SHADOW_SCANS = 30          # a band with a peak over twice its height this close before it is a satellite
MIN_OFFSET = 3             # scans; a sample at the standard's own scan is a front or bleed, not a product
ANCHOR_BONUS = 2.0         # placement score for a sample peak beside a homoduplex
PARTNER_BONUS = 1.5        # ...and for the partner heteroduplex where the plate puts it
MAX_BLEED = 0.5            # never take off more than this x the standard
MIN_WELLS_AGREE = 0.20     # share of wells that must reproduce the standard's spacing
MERGED_BELOW = 8           # scans; a narrower heteroduplex pair is one band
GAP_SD = 0.18              # log-ratio spread of a well's gaps around the plate's
PARTNER_TOL = 0.3          # relative, where the second heteroduplex may sit
OFFSET_MAX_FRAC = 0.45     # sample offset must be under this x g1

# calling (chosen with each plate left out; stable across all 12 plates)
MIN_PRODUCT_SIGMA = 40.0   # stronger homoduplex below this -> no call
MINOR_WITH_HET = 0.05      # minor homoduplex share that counts with heteroduplexes
HET_SHARE = 0.20           # heteroduplex band vs the stronger homoduplex
TWO_HOMS = 0.40            # minor/major homoduplex: two alleles on its own
AI_BELOW = 0.33            # a het with minor/major below this is flagged AI


# ---------------------------------------------------------------- signal
def flatten(y: np.ndarray) -> np.ndarray:
    """The trace minus a smoothed rolling-minimum baseline."""
    y = np.asarray(y, float)
    return y - uniform_filter1d(minimum_filter1d(y, BASELINE_WINDOW), BASELINE_WINDOW)


def noise_sigma(y: np.ndarray) -> float:
    """Scan-to-scan noise (robust), so heights can be read in noise units."""
    d = np.diff(np.asarray(y, float))
    return float(1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2)) + 1e-9


def _strong_peaks(y, sig, prominence, k):
    p, pr = find_peaks(y, prominence=prominence * sig, distance=2)
    keep = p >= CUT
    p, pr = p[keep], pr["prominences"][keep]
    return np.sort(p[np.argsort(-pr)[:k]])


# ---------------------------------------------------------------- results
@dataclass
class PlateStandard:
    """The plate's standard: channels, spacing and how it was found."""
    sample_col: int
    is_col: int
    gaps: np.ndarray                 # g1, g2, g3 (g3 = 0: merged pair)
    wells_agree: int
    wells_total: int
    offset: float = 0.0              # sample scan - standard scan
    note: str = ""

    @property
    def merged(self) -> bool:
        return self.gaps[2] < MERGED_BELOW


@dataclass
class WellCall:
    well: str
    call: str = "no-call"            # hom1 / hom2 / het / no-call
    frac2: float = float("nan")      # Homo2 share of the homoduplexes, 0..1
    flags: List[str] = field(default_factory=list)
    reason: str = ""
    std_scans: List[int] = field(default_factory=list)
    homo1: float = 0.0               # sample peak heights, noise units
    homo2: float = 0.0
    hetero: List[float] = field(default_factory=list)
    sample_sigma: float = 0.0

    def as_row(self, run=""):
        return {"run": run, "well": self.well, "call": self.call,
                "frac2": None if np.isnan(self.frac2) else round(self.frac2, 4),
                "flags": ",".join(self.flags), "reason": self.reason,
                "std_scans": "/".join(str(s) for s in self.std_scans),
                "homo1": round(self.homo1, 1), "homo2": round(self.homo2, 1),
                "hetero": "/".join(f"{h:.1f}" for h in self.hetero)}


# ---------------------------------------------------------------- the vote
def _is_satellite(y, comb):
    """True when the 4th band is only the companion every band has (the
    shoulder ~20 scans behind each standard peak), not a second
    heteroduplex: Homo1 and Homo2 show the same companion at the same
    distance and relative height."""
    d = int(comb[3] - comb[2])
    r = y[comb[3]] / max(y[comb[2]], 1e-9)

    def companion(p):
        lo, hi = p + d - 2, p + d + 3
        return y[lo:hi].max() / max(y[p], 1e-9) if hi <= len(y) else 0.0
    return companion(comb[0]) >= 0.6 * r and companion(comb[1]) >= 0.6 * r


def _shadowed(y, p):
    """True when *p* is only a shoulder or satellite: a peak over twice its
    height sits within SHADOW_SCANS before it."""
    lo = max(0, p - SHADOW_SCANS)
    return lo < p - 2 and y[lo:p - 2].max() > 2.0 * y[p]


def _anchored(s, s_sig, bands):
    """The sample is this fragment, so within the stretch the standard spans
    the sample's biggest peak sits just before Homo1 or Homo2 (in a
    heterozygote the heteroduplexes may match it, rarely much more)."""
    a, b = ANCHOR_WINDOW
    near = max(s[max(0, bands[0] - a):bands[0] + b].max(),
               s[max(0, bands[1] - a):bands[1] + b].max())
    lo, hi = max(0, bands[0] - 40), min(len(s), bands[-1] + 40)
    return near >= max(20 * s_sig, ANCHOR_FRAC * s[lo:hi].max())


def _vote(ys, ss, n_bands, anchored=True):
    """Best spacing over all wells' proposals -> (value, gaps, n_wells) or None."""
    V, Wi, Hh = [], [], []
    for i, w in enumerate(ys):
        y, sig = ys[w]
        P = _strong_peaks(y, sig, VOTE_PROMINENCE, VOTE_PEAKS)
        s, s_sig, s_top = ss[w]
        for comb in itertools.combinations(P, n_bands):
            g = np.diff(comb).astype(float)
            if not (G1_RANGE[0] <= g[0] <= G1_RANGE[1]
                    and max(G2_RANGE[0], G2_MIN_RATIO * g[0]) <= g[1] <= G2_RANGE[1]):
                continue
            if n_bands == 4 and not (4 <= g[2] <= 0.8 * g[1]):
                continue
            if n_bands == 4 and _is_satellite(y, comb):
                continue
            if any(_shadowed(y, x) for x in comb[:3]):
                continue
            hs = [y[x] for x in comb]
            if min(hs) < EQUIMOLAR_MIN * max(hs):
                continue
            if anchored and not _anchored(s, s_sig, comb):
                continue
            V.append(g); Wi.append(i)
            Hh.append(np.log1p(max(min(hs), 0.0) / sig))
    if not V:
        return None
    V, Wi, Hh = np.array(V), np.array(Wi), np.array(Hh)
    rel = np.array(VOTE_TOL[:n_bands - 1]); floor = np.array(VOTE_TOL_MIN[:n_bands - 1])
    best = None
    for j in range(len(V)):
        t = np.maximum(floor, rel * V[j])
        m = np.all(np.abs(V - V[j]) <= t, axis=1)
        support = np.zeros(len(ys))
        np.maximum.at(support, Wi[m], Hh[m])
        val = support.sum()
        if best is None or val > best[0]:
            best = (val, np.median(V[m], axis=0), int((support > 0).sum()))
    return best


def _pattern(ys, ss, anchored=True):
    b4, b3 = _vote(ys, ss, 4, anchored), _vote(ys, ss, 3, anchored)
    # A four-band pattern must be nearly as strong as the three-band one: the
    # vote sums each well's weakest band, so a fourth band that is only a
    # satellite (a third of a real band) pulls it well below.
    # When the three-band winner is only the four-band pattern with Hetero1
    # skipped (Hetero2 is often the taller), the four-band pattern explains it.
    skips = (b4 is not None and b3 is not None
             and abs(b3[1][0] - b4[1][0]) <= max(8, 0.12 * b4[1][0])
             and abs(b3[1][1] - (b4[1][1] + b4[1][2])) <= max(8, 0.12 * b3[1][1]))
    if b4 is not None and (b3 is None or (b4[2] >= FOUR_BAND_SHARE * b3[2]
                                          and (skips or b4[0] >= FOUR_BAND_STRENGTH * b3[0]))):
        return b4[0], b4[1], b4[2]
    if b3 is None:
        return None
    return b3[0], np.array([b3[1][0], b3[1][1], 0.0]), b3[2]


def _partner(y, p3, he, g3, ref):
    """The second heteroduplex for a band placed as Hetero1: ``(scan, side)``
    with side +1 when *he* is Hetero1 (partner after it), -1 when *he* is
    really Hetero2 (partner before it), or ``None``."""
    tol = max(4, PARTNER_TOL * g3)
    after = [x for x in p3 if abs((x - he) - g3) <= tol and y[x] >= 0.2 * ref]
    before = [x for x in p3 if abs((he - x) - g3) <= tol and y[x] >= 0.2 * ref]
    a = max(after, key=lambda x: y[x]) if after else None
    b = max(before, key=lambda x: y[x]) if before else None
    if a is not None and (b is None or y[a] >= y[b]):
        return int(a), 1
    if b is not None:
        return int(b), -1
    return None


def place_standard(y, sig, gaps, s=None, s_sig=1.0) -> Optional[List[int]]:
    """The standard's bands in one well: [H1, H2, Het1(, Het2)] or None.

    Every ordered (H1, H2, Het) set among the well's strongest peaks is
    scored by its weakest and mean band height, minus a soft penalty for gaps
    off the plate's, plus a bonus when the sample sits beside a homoduplex
    and -- on a plate whose heteroduplexes separate -- when the partner
    heteroduplex is where the plate puts it.  A primer or dye front can
    mimic three bands; it rarely has the pair."""
    P = _strong_peaks(y, sig, PLACE_PROMINENCE, PLACE_PEAKS)
    if len(P) < 3:
        return None
    H = np.log1p(np.maximum(y[P], 0.0) / sig)
    G = np.asarray(gaps[:2], float)
    g3 = float(gaps[2])
    pair = g3 >= MERGED_BELOW
    p3 = find_peaks(y, prominence=3 * sig, distance=2)[0] if pair else None
    best = None
    shadowed = {int(x): _shadowed(y, x) for x in P}
    partners = {}
    for comb in itertools.combinations(range(len(P)), 3):
        pos = P[list(comb)]
        if any(shadowed[int(x)] for x in pos):
            continue
        pen = np.sum((np.log(np.diff(pos) / G) / GAP_SD) ** 2)
        hh = H[list(comb)]
        bonus = 0.0
        if s is not None and _anchored(s, s_sig, pos):
            bonus += ANCHOR_BONUS
        part = None
        if pair:
            # the partner's height floor is set by the weakest band; keyed on
            # Hetero1 and that floor, so it is computed once per pair
            ref = min(y[x] for x in pos)
            k = (int(pos[2]), round(float(ref), 3))
            if k not in partners:
                partners[k] = _partner(y, p3, int(pos[2]), g3, ref)
            part = partners[k]
            if part is not None:
                bonus += PARTNER_BONUS
        sc = hh.min() + 0.5 * hh.mean() - 0.5 * pen + bonus
        if best is None or sc > best[0]:
            best = (sc, [int(x) for x in pos], part)
    if best is None:
        return None
    _, hits, part = best
    if part is not None:
        x, side = part
        hits = hits + [x] if side > 0 else hits[:2] + [x, hits[2]]
    return hits


def unbleed(s, y, bands):
    """The sample channel with the standard's crosstalk taken off.

    The standard bleeds into the sample channel at its own scans, while a
    product sits several scans to one side.  So at the standard's band scans
    the sample channel holds bleed (plus at most a product's tail), and the
    smallest sample/standard ratio over the bands is the bleed fraction."""
    ratios = [max(s[p], 0.0) / y[p] for p in bands if y[p] > 0]
    if not ratios:
        return s
    k = min(min(ratios), MAX_BLEED)
    return s - k * np.clip(y, 0.0, None)


def _check_offset(ys, ss, gaps, n=32):
    """The sample offset a candidate standard implies, from up to *n* wells."""
    wells = list(ys)[::max(1, len(ys) // n)]
    placed = {w: place_standard(ys[w][0], ys[w][1], gaps, ss[w][0], ss[w][1]) for w in wells}
    return _sample_offset(placed, ss)


def _sample_offset(placed, ss):
    offs = []
    for w, hits in placed.items():
        if not hits:
            continue
        s = ss[w][0]
        lo, hi = max(0, hits[0] - 60), hits[1] + 60
        p = lo + int(np.argmax(s[lo:hi]))
        d = min((p - hits[0], p - hits[1]), key=abs)
        if abs(d) < OFFSET_MAX_FRAC * (hits[1] - hits[0]):
            offs.append(d)
    return float(np.median(offs)) if offs else 0.0


# ---------------------------------------------------------------- calling
def call_well(homo1, homo2, hetero) -> tuple:
    """(call, frac2, flags, reason) from sample heights in noise units."""
    he = max(hetero) if hetero else 0.0
    big = max(homo1, homo2)
    if big < MIN_PRODUCT_SIGMA:
        return ("no-call", float("nan"), [],
                f"no product: strongest homoduplex {big:.0f}x noise "
                f"(needs {MIN_PRODUCT_SIGMA:.0f})")
    frac2 = homo2 / (homo1 + homo2)
    minor = min(homo1, homo2) / big
    het_share = he / big
    if minor >= TWO_HOMS or (minor >= MINOR_WITH_HET and het_share >= HET_SHARE):
        flags = ["ai"] if minor < AI_BELOW else []
        return "het", frac2, flags, ""
    return ("hom1" if homo1 > homo2 else "hom2"), frac2, [], ""


def _channel_cols(arrs, base_order, is_channel, sample_channel):
    if sample_channel:
        sample_col = acgt_index_for_channel(base_order, sample_channel)
    else:
        tops = [int(np.argmax(np.percentile(a[CUT:], 99.9, axis=0))) for a in arrs.values()]
        sample_col = max(set(tops), key=tops.count)
    if is_channel:
        return sample_col, [acgt_index_for_channel(base_order, is_channel)]
    n = next(iter(arrs.values())).shape[1]
    return sample_col, [c for c in range(n) if c != sample_col]


def _channel_of_col(base_order, col):
    for ch in range(1, 5):
        if acgt_index_for_channel(base_order, ch) == col:
            return ch
    return col + 1


def find_plate_standard(arrs: Dict[str, np.ndarray], base_order="ACTG",
                        is_channel=None, sample_channel=None) -> Optional[PlateStandard]:
    """Channels and spacing of the plate's standard, or None."""
    arrs = {w: np.asarray(a, float) for w, a in arrs.items()}
    arrs = {w: a for w, a in arrs.items() if a.ndim == 2 and a.shape[0] > CUT + 500}
    if not arrs:
        return None
    sample_col, is_cols = _channel_cols(arrs, base_order, is_channel, sample_channel)
    ss = {}
    for w, a in arrs.items():
        s = flatten(a[:, sample_col])
        ss[w] = (s, noise_sigma(s), float(s[CUT:].max()))
    ys_by = {col: {w: (lambda y: (y, noise_sigma(y)))(flatten(a[:, col]))
                   for w, a in arrs.items()} for col in is_cols}
    # Every channel but the sample's is a candidate.  A candidate counts only
    # if enough wells agree on its spacing and, with a product to check
    # against, the sample sits *beside* it: the two dyes run at different
    # speeds, so a product leads or trails its standard by several scans,
    # while a bleed copy of the sample, a dye front or leftover crosstalk
    # sits at the same scan (offset ~0).  Of those, the strongest vote wins.
    # A plate whose PCR failed has no product to anchor on, so the standard
    # is then looked for on its own (every well will be a no-call with the
    # standard still placed).
    best, unanchored = None, False
    for anchored in (True, False):
        for col in is_cols:
            got = _pattern(ys_by[col], ss, anchored)
            if not got or got[2] < MIN_WELLS_AGREE * len(arrs):
                continue
            if anchored and abs(_check_offset(ys_by[col], ss, got[1])) < MIN_OFFSET:
                continue
            if best is None or got[0] > best[0]:
                best = (got[0], col, got[1], got[2])
        if best is not None:
            unanchored = not anchored
            break
    if best is None:
        return None
    _, col, gaps, n = best
    ps = PlateStandard(sample_col=sample_col, is_col=col, gaps=gaps,
                       wells_agree=n, wells_total=len(arrs))
    ps.note = (f"standard on Ch{_channel_of_col(base_order, col)}, sample on "
               f"Ch{_channel_of_col(base_order, sample_col)}; spacing "
               f"{gaps[0]:.0f}/{gaps[1]:.0f}/"
               + ("merged" if ps.merged else f"{gaps[2]:.0f}")
               + f" scans, found in {n} of {len(arrs)} wells")
    if unanchored:
        ps.note += " (no sample product to check it against)"
    return ps


def corrected_plate(docs, base_order="ACTG", matrix="peaks"):
    """``{well: acgt}`` with the dye overlap removed.

    *matrix* is a :class:`spectral.SpectralMatrix`, or "peaks" to measure the
    plate's own (:func:`spectral.matrix_from_peaks`).  Needs the documents'
    physical ``channels``.  Returns ``(arrays, matrix)``."""
    import spectral
    chans = {w: np.asarray(d.channels, float) for w, d in docs.items()}
    m = matrix
    if isinstance(matrix, str):
        m = spectral.matrix_from_peaks(chans, cut=spectral.PEAKS_CUT) or spectral.identity(
            next(iter(chans.values())).shape[1])
    arrs = {}
    for w, x in chans.items():
        y = m.correct(spectral._flat(x))
        a = np.zeros((len(y), 4))
        for ch in range(1, 5):
            a[:, acgt_index_for_channel(base_order, ch)] = y[:, ch - 1]
        arrs[w] = a
    return arrs, m


def genotype_plate(docs, base_order="ACTG", is_channel=None, sample_channel=None,
                   run_name="", matrix=None) -> Dict[str, dict]:
    """Genotype every well of one plate with no marks -> {well: row}.

    *docs* maps well -> a trace document with an ``acgt`` array (or the array
    itself).  *is_channel* / *sample_channel* (physical 1..4) fix the roles
    when known; otherwise they are found.  *matrix* (None: raw traces) is a
    dye-overlap matrix or "peaks" for the plate's own, applied first (see
    :func:`corrected_plate`; raw traces still score slightly better, see
    docs/GENOTYPING_ENGINE.md).  Rows carry ``plate_note``.
    """
    mnote = ""
    if matrix is not None:
        arrs, m = corrected_plate(docs, base_order, matrix)
        mnote = f"; dye overlap removed ({m.name})"
    else:
        arrs = {w: np.asarray(getattr(d, "acgt", d), float) for w, d in docs.items()}
    ps = find_plate_standard(arrs, base_order, is_channel, sample_channel)
    if ps is not None:
        ps.note += mnote
    out = {}
    if ps is None:
        for w in arrs:
            row = WellCall(w, reason="no internal standard found on this plate").as_row(run_name)
            row["plate_note"] = "no internal standard found on this plate"
            out[w] = row
        return out
    placed, ss = {}, {}
    for w, a in arrs.items():
        if a.ndim != 2 or a.shape[0] <= CUT + 100:
            placed[w] = None
            continue
        y = flatten(a[:, ps.is_col]); s = flatten(a[:, ps.sample_col])
        placed[w] = place_standard(y, noise_sigma(y), ps.gaps, s, noise_sigma(s))
        if placed[w]:
            s = unbleed(s, y, placed[w])
        ss[w] = (s, noise_sigma(s), float(s[CUT:].max()))
    ps.offset = _sample_offset({w: h for w, h in placed.items() if w in ss}, ss)
    ps.note += f"; sample runs {ps.offset:+.0f} scans from the standard"
    g1 = ps.gaps[0]
    half_h = max(3, int(0.3 * min(g1, 40)))
    for w in arrs:
        wc = WellCall(w)
        hits = placed.get(w)
        if not hits:
            wc.reason = "internal standard not found in this well"
        else:
            s, sig, _ = ss[w]
            wc.std_scans = hits
            wc.sample_sigma = sig

            def peak(c, half):
                c = int(round(c + ps.offset))
                lo, hi = max(0, c - half), min(len(s), c + half + 1)
                return float(max(s[lo:hi].max(), 0.0)) / sig if hi > lo else 0.0

            wc.homo1, wc.homo2 = peak(hits[0], half_h), peak(hits[1], half_h)
            if ps.merged:
                wc.hetero = [peak(hits[2] + ps.gaps[2] / 2.0, max(4, half_h) + int(ps.gaps[2]))]
                wc.flags.append("merged-het")
            else:
                wc.hetero = [peak(x, max(4, half_h)) for x in hits[2:]]
            call, frac2, flags, reason = call_well(wc.homo1, wc.homo2, wc.hetero)
            wc.call, wc.frac2, wc.reason = call, frac2, reason
            wc.flags = flags + [f for f in wc.flags if f not in flags]
        row = wc.as_row(run_name)
        row["plate_note"] = ps.note
        out[w] = row
    return out
