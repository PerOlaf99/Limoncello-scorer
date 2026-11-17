"""
Cimarron-3.12-compatible calling stage (clean room): processed 4-channel trace
-> bases, peak positions, qualities and good-sequence limits.

Follows the functional specification SPEC_calling.md (provenance of each
constant in CLEANROOM_rewrite.md):

* the trace is cut into 2048-scan windows (step 1900, last one end-aligned);
* per window: envelope peaks -> bands (trough, peak, trough); a fuzzy
  omission classifier drops artefact bands; a fuzzy gap splitter inserts bands
  into gaps that span several expected spacings; a second omission pass;
* letters from the dominant channel (N when two channels are within 80 %);
* fuzzy quality scores; windows are stitched; the read end is cut 50 scans
  before the trace end; good-sequence limits from smoothed qualities.

Input scale: the rules assume an *unflattened* trace whose envelope has a mean
of 0.5 over the read.  ``call_trace(..., flattened=False)`` takes such an
internal trace (it is rescaled to that mean); ``flattened=True`` takes a
display-flattened trace (Cimarron's .esd convention) and first undoes the
flattening approximately (fixed-iteration eigenvector estimate).

Positions: processed index i is raw scan ``start_tp - 1 + i``; the reported
position is ``i + start_tp + 1``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations

import numpy as np
from scipy.ndimage import uniform_filter1d

WIN_LEN = 2048
WIN_STEP = 1900
END_MARGIN = 50
EPS = np.finfo(float).eps

# ---------------------------------------------------------------------------
# Fuzzy helpers
# ---------------------------------------------------------------------------


def mf(x, pts):
    """Piecewise-linear membership through (x, y) breakpoints, constant outside."""
    xs = np.array([p[0] for p in pts], float)
    ys = np.array([p[1] for p in pts], float)
    return np.interp(np.asarray(x, float), xs, ys)


def very(m):
    return m * m


def somewhat(m):
    return np.sqrt(np.clip(m, 0, None))


def contrast(m):
    return np.where(m >= 0.5, np.sqrt(np.clip(m, 0, None)), m * m)


def _AND(*a):
    return np.minimum.reduce([np.asarray(x, float) for x in a])


def _OR(*a):
    return np.maximum.reduce([np.asarray(x, float) for x in a])


def _NOT(a):
    return 1.0 - np.asarray(a, float)


class Defuzzer:
    """Centroid of the pointwise max of strength-scaled output sets over a
    domain, evaluated on a fine grid (trapezoid rule)."""

    def __init__(self, sets, lo, hi, n=4001):
        self.x = np.linspace(lo, hi, n)
        self.sets = np.array([mf(self.x, s) for s in sets])      # (R, G)
        dx = np.diff(self.x)
        self.wt = np.zeros(n)
        self.wt[:-1] += dx / 2
        self.wt[1:] += dx / 2

    def __call__(self, strengths):
        s = np.atleast_2d(np.asarray(strengths, float))            # (B, R)
        agg = np.max(s[:, :, None] * self.sets[None, :, :], axis=1)  # (B, G)
        area = agg @ self.wt
        mom = agg @ (self.wt * self.x)
        return np.where(area > 0, mom / np.where(area > 0, area, 1.0), 0.0)


def round_half_up(x):
    return np.floor(np.asarray(x, float) + 0.5)


# ---------------------------------------------------------------------------
# Parameters (SPEC_calling §10, §11, §13)
# ---------------------------------------------------------------------------

OMIT_PARAMS = {
    1: dict(a=(-0.10, 0.33, 0.73, 1.07), t=(1.0, 1.1, 1.8), c=(1.70, 1.97), tau=(0.47, 1.00), kap=(0.70, 1.50)),
    2: dict(a=(-0.10, 0.00, 0.70, 1.10), t=(0.7, 1.0, 1.5), c=(1.46, 2.23), tau=(0.50, 1.07), kap=(0.87, 1.80)),
}
OMIT_DEFUZZ = Defuzzer([
    [(0.4596, 1), (1.3797, 1), (1.6864, 0)],
    [(1.3797, 0), (1.6864, 1), (2.2998, 1), (2.6065, 0)],
    [(2.2998, 0), (2.6065, 1), (3.5540, 1)],
], 0.4596, 3.5540)

SPLIT_DEFUZZ = Defuzzer([
    [(0.5079, 1), (1.5002, 1), (1.5069, 0)],
    [(1.5002, 0), (1.5069, 1), (2.5063, 1)],
], 0.5079, 2.5063)
SUBGAP_OFFSET = 0.6

_TRI = [(-0.100, 0.010, 0.125), (0.010, 0.125, 0.250), (0.125, 0.250, 0.375), (0.250, 0.375, 0.500),
        (0.375, 0.500, 0.625), (0.500, 0.625, 0.750), (0.625, 0.750, 0.875), (0.750, 0.875, 0.990),
        (0.875, 0.990, 1.105)]
QUAL_DEFUZZ = Defuzzer([[(a, 0), (b, 1), (c, 0)] for a, b, c in _TRI], 0.0, 1.0, n=10001)

N_FRACTION = 0.80


# ---------------------------------------------------------------------------
# Per-scan helper quantities (§6)
# ---------------------------------------------------------------------------

def helpers(trace):
    v = np.sort(trace, axis=1)
    v1, v3, v4 = v[:, 0], v[:, 2], v[:, 3]
    E = v4
    dom = np.argmax(trace, axis=1)
    v3e = np.maximum(v3, EPS)
    X = np.where(v4 < 0.1, 1.0 + v4 / np.sqrt(v3e), v4 / v3e)
    X = np.clip(X, 1.0, 2.0)
    den = v4 - v3
    with np.errstate(divide="ignore", invalid="ignore"):
        B = (np.maximum(v3, 0) - np.maximum(v1, 0)) / den
    B = np.where(den == 0, 0.5, np.minimum(B, 9.99))
    return E, dom, X, B


# ---------------------------------------------------------------------------
# Bands
# ---------------------------------------------------------------------------

@dataclass
class Bands:
    l: np.ndarray
    p: np.ndarray
    r: np.ndarray
    ins: np.ndarray            # inserted by the gap splitter
    rebuilt: bool = False

    def __len__(self):
        return len(self.p)

    @property
    def width(self):
        return self.r - self.l + (1 if self.rebuilt else 0)

    def take(self, keep):
        keep = np.asarray(keep)
        return Bands(self.l[keep].copy(), self.p[keep].copy(), self.r[keep].copy(), self.ins[keep].copy(),
                     self.rebuilt)


def detect_bands(E):
    """Three-state peak/trough walk over the envelope (§7)."""
    n = len(E)
    peaks, troughs = [], []
    state = 0          # 0 unknown, 1 rising, -1 falling
    for t in range(1, n):
        if E[t] > E[t - 1]:
            if state == -1:
                troughs.append(t - 1)
            state = 1
        elif E[t] < E[t - 1]:
            if state == 1:
                peaks.append(t - 1)
            state = -1
    peaks = np.array(peaks, int)
    troughs = np.array(troughs, int)
    if len(troughs) < 2 or len(peaks) == 0:
        return Bands(*(np.zeros(0, int),) * 3, np.zeros(0, bool))
    peaks = peaks[(peaks > troughs[0]) & (peaks < troughs[-1])]
    k = np.searchsorted(troughs, peaks)        # troughs[k-1] < peak < troughs[k]
    return Bands(troughs[k - 1], peaks, troughs[k], np.zeros(len(peaks), bool))


def median_int(g):
    g = np.sort(np.asarray(g, int))
    if len(g) == 0:
        return 0
    m = len(g)
    if m % 2:
        return int(g[m // 2])
    return int((g[m // 2 - 1] + g[m // 2]) // 2)


def gaps(b: Bands):
    """Left and right gaps (peak differences); the outer ones are the median gap."""
    d = np.diff(b.p)
    med = median_int(d) if len(d) else 0
    left = np.concatenate([[med], d]).astype(int)
    right = np.concatenate([d, [med]]).astype(int)
    return left, right


def rebuild(b: Bands) -> Bands:
    """Make adjacent band boundaries touch (§7): the wider band's boundary moves
    to the narrower band's boundary."""
    order = np.argsort(b.p, kind="stable")
    b = b.take(order)
    l, r = b.l.copy(), b.r.copy()
    for j in range(len(l) - 1):
        if r[j] != l[j + 1]:
            if (r[j] - l[j]) >= (r[j + 1] - l[j + 1]):
                r[j] = l[j + 1]
            else:
                l[j + 1] = r[j]
    return Bands(l, b.p.copy(), r, b.ins.copy(), True)


# ---------------------------------------------------------------------------
# Expected-spacing / width curves (§8)
# ---------------------------------------------------------------------------

def expected_curve(p, y):
    p = np.asarray(p, float)
    y = np.asarray(y, float)
    n = len(p)
    if n == 0:
        return np.zeros(0)
    if n < 4:
        return np.full(n, round_half_up(np.median(y)))
    V = np.column_stack([np.ones(n), p, p * p])
    coef = np.linalg.lstsq(V, y, rcond=None)[0]
    pred = round_half_up(V @ coef)
    s = np.sqrt(np.sum((y - V @ coef) ** 2) / (n - 3))
    inl = np.abs(pred - y) < s
    if inl.sum() >= 4:
        coef = np.linalg.lstsq(V[inl], pred[inl], rcond=None)[0]
        pred = round_half_up(V @ coef)
    a, bb, c = coef
    if c != 0:
        v = -bb / (2 * c)
        if c > 0:
            if v > p[-1]:
                v = p[-1]
            if p[0] <= v <= p[-1]:
                k = int(np.argmin(np.abs(p - v)))
                pred = np.where(p <= v, pred[k], pred)
        else:
            if p[0] <= v <= p[-1]:
                k = int(np.argmin(np.abs(p - v)))
                pred = np.where(p >= v, pred[k], pred)
    return pred


# ---------------------------------------------------------------------------
# Omission classifier (§10)
# ---------------------------------------------------------------------------

def omission(b: Bands, E, X, s_exp, which):
    """Return a keep mask."""
    n = len(b)
    if n == 0:
        return np.zeros(0, bool)
    P = OMIT_PARAMS[which]
    valley = np.minimum(E[b.l], E[b.r])
    mbar = float(np.mean(valley))
    h = E[b.p]
    gl, gr = gaps(b)
    s = np.maximum(s_exp.astype(float), 1.0)

    def phase(g):
        return np.where(g < s / 2, 0.5, np.mod(g, s) / s)
    phl, phr = phase(gl), phase(gr)
    on = np.maximum(mf(phl, [(0.2, 1), (0.5, 0), (0.8, 1)]), mf(phr, [(0.2, 1), (0.5, 0), (0.8, 1)]))
    a1, a2, a3, a4 = P["a"]
    offp = [(a1, 0), (a2, 1), (a3, 1), (a4, 0)]
    off = very(np.maximum(mf(phl, offp), mf(phr, offp)))
    x = X[b.p]
    t1, t2, t3 = P["t"]
    comp = mf(x, [(t1, 1), (t2, 0.5), (t3, 0)])
    c1, c2 = P["c"]
    clear = mf(x, [(c1, 0), (c2, 1)])
    tau1, tau2 = P["tau"]
    k1, k2 = P["kap"]
    if mbar > 0:
        tiny = mf(h, [(tau1 * mbar, 1), (tau2 * mbar, 0)])
        tall = somewhat(mf(h, [(k1 * mbar, 0), (k2 * mbar, 1)]))
    else:
        tiny = (h <= 0).astype(float)
        tall = (h > 0).astype(float)
    keep_r = _AND(clear, _OR(tall, on))
    amb_r = _AND(comp, _OR(tall, _AND(on, tiny)))
    drop_r = _AND(tiny, off)
    concl = OMIT_DEFUZZ(np.column_stack([keep_r, amb_r, drop_r]))
    cls = round_half_up(concl)
    return (cls == 1) | (cls == 2)


# ---------------------------------------------------------------------------
# Gap splitter (§11)
# ---------------------------------------------------------------------------

def _q(nG, nC):
    return 1.0 - np.sqrt(((5 - nG) ** 2 + (5 - nC) ** 2) / 50.0)


Q22 = _q(2, 2)


def split_and_insert(b: Bands, E, dom, s_exp, w_exp, g_idx, c_idx):
    n = len(b)
    if n < 2:
        return b
    gl, _ = gaps(b)
    wid = b.width.astype(float)
    letters = dom[b.p]
    s = np.maximum(s_exp.astype(float), 1.0)
    w = np.maximum(w_exp.astype(float), 1.0)
    j = np.arange(1, n)
    r_gap = gl[j] / s[j] - 1
    r_prevgap = np.where(j - 1 == 0, 0.0, gl[j - 1] / s[j] - 1)
    r_wid = wid[j] / w[j] - 1
    r_prevwid = wid[j - 1] / w[j] - 1
    rho = np.zeros(len(j))
    for t, jj in enumerate(j):
        a = max(0, jj - 5)
        pre = letters[a:jj]
        nG = int(np.sum(pre == g_idx))
        nC = int(np.sum(pre == c_idx))
        rp = min(1.0, _q(nG, nC) / Q22) ** 2
        Pg = float(np.sum(gl[a:jj]))
        Qw = w[jj] * (jj - a)
        fac = 1.0 if (Qw == 0 or Pg < Qw) else Qw / Pg
        rho[t] = fac * rp
    big = very(mf(r_gap, [(0.22, 0), (0.60, 1), (1.00, 1)]))
    small = contrast(mf(r_prevgap, [(-1.0, 1), (-0.40, 1), (0.16, 0)]))
    wj = contrast(mf(r_wid, [(-0.10, 0), (0.70, 1)]))
    wp = contrast(mf(r_prevwid, [(-0.10, 0), (0.70, 1)]))
    normal = _OR(_NOT(big), _AND(big, rho), _AND(big, small, _NOT(wj), _NOT(wp)))
    split = _OR(_AND(big, _OR(wj, wp)), _AND(big, _NOT(small), _NOT(rho)))
    concl = SPLIT_DEFUZZ(np.column_stack([normal, split]))
    do = round_half_up(concl) == 2
    L, P, R, I = [b.l], [b.p], [b.r], [b.ins]
    m = len(E)
    for t in np.where(do)[0]:
        jj = j[t]
        G = int(b.p[jj] - b.p[jj - 1])
        sp = s[jj - 1]
        k = int(np.floor(G / sp + SUBGAP_OFFSET))
        if k < 2:
            continue
        u = G // k
        if u < 1:
            continue
        for i in range(1, k):
            c = int(b.p[jj - 1] + i * u)
            lo = c - u // 2
            hi = lo + u - 1
            lo_c, hi_c = max(lo, 0), min(hi, m - 1)
            pk = c
            if hi_c > lo_c:
                f = E[lo_c:hi_c + 1] - E[lo_c:hi_c + 1].min()
                f0, f1 = f[:-1], f[1:]
                area = (f0 + f1) / 2
                x0 = np.arange(lo_c, hi_c, dtype=float)
                A = area.sum()
                if A > 0:
                    mom = np.sum(x0 * area + f0 / 6 + f1 / 3)
                    pk = int(round_half_up(mom / A))
            L.append(np.array([lo])); P.append(np.array([pk])); R.append(np.array([hi]))
            I.append(np.array([True]))
    nb = Bands(np.concatenate(L).astype(int), np.concatenate(P).astype(int), np.concatenate(R).astype(int),
               np.concatenate(I).astype(bool), b.rebuilt)
    order = np.argsort(nb.p, kind="stable")
    nb = nb.take(order)
    # drop exact duplicates of a peak position (an insert landing on a band)
    if len(nb) > 1:
        keep = np.concatenate([[True], np.diff(nb.p) > 0])
        nb = nb.take(np.where(keep)[0])
    return rebuild(nb)


# ---------------------------------------------------------------------------
# Quality (§13)
# ---------------------------------------------------------------------------

def _quad_feats(y):
    x = np.arange(1, len(y) + 1, dtype=float)
    V = np.column_stack([np.ones_like(x), x, x * x])
    coef, *_ = np.linalg.lstsq(V, y, rcond=None)
    res = y - V @ coef
    dof = max(len(y) - 3, 1)
    return np.array([coef[0], coef[1], coef[2], float(res @ res) / dof])


def _pearson(a, b):
    a = a - a.mean(); b = b - b.mean()
    d = np.sqrt((a @ a) * (b @ b))
    return float(a @ b / d) if d > 0 else np.nan


def shape_scores(b: Bands, E, sp):
    n = len(b)
    wid = b.width
    F = min(100, int(np.floor(0.5 + 3.714 * 2048 / (2 * np.pi * max(sp, 1e-9)))))
    F = max(F, 2)
    sig_ref = 2048.0 / (2 * np.pi * (F - 1))
    M = max(4, int((median_int(wid) + 1) // 2))
    xr = np.arange(M) - (M - 1) / 2.0
    ref = np.exp(-xr ** 2 / (2 * (sig_ref / 2) ** 2))
    rf = _quad_feats(ref)
    out = np.full(n, -1.0)
    for i in range(n):
        if wid[i] <= 4:
            continue
        lo, hi = max(b.l[i], 0), min(b.r[i], len(E) - 1)
        y = E[lo:hi + 1]
        if len(y) < 4:
            continue
        y = y - y.min()
        c1 = _pearson(rf, _quad_feats(y))
        c2 = _pearson(rf, _quad_feats(y[::-1]))
        v = np.nanmax([c1, c2]) if not (np.isnan(c1) and np.isnan(c2)) else -1.0
        out[i] = v
    return out


def qualities(b: Bands, E, X, B, s_exp, w_exp, sp):
    n = len(b)
    if n == 0:
        return np.zeros(0)
    h = E[b.p]
    x = X[b.p]
    bz = B[b.p]
    r = b.width / np.maximum(w_exp.astype(float), 1.0) - 1
    shape = shape_scores(b, E, sp)
    gl, gr = gaps(b)
    s = np.maximum(s_exp.astype(float), 1.0)
    Lg = np.maximum(gl, gr).astype(float)
    Sg = np.minimum(gl, gr).astype(float)

    def fix(v):
        return np.where(v < s / 4, s / 2 - v * s / 4, v)
    phL = np.mod(fix(Lg), s) / s
    phS = np.mod(fix(Sg), s) / s
    tiny = mf(h, [(0, 1), (0.02, 1), (0.03, 0)])
    small = mf(h, [(0.015, 0), (0.025, 1), (0.04, 1), (0.05, 0)])
    moder = mf(h, [(0.035, 0), (0.045, 1), (0.14, 1), (0.16, 0)])
    normal = mf(h, [(0.135, 0), (0.145, 1), (0.48, 1), (0.55, 0)])
    tall = mf(h, [(0.475, 0), (0.485, 1), (1.0, 1)])
    xpres = mf(x, [(1.0, 1), (1.35, 0)])
    xneg = mf(x, [(1.35, 0), (1.5, 1)])
    Wm = mf(r, [(-1, 0), (-0.4, 1), (0.4, 1), (1.0, 0)])
    Gm = mf(shape, [(0.4, 0), (0.5, 1), (1.0, 1)])
    Qm = mf(bz, [(0, 1), (0.2, 1), (0.4, 0)])
    og = [(0.3, 1), (0.5, 0), (0.7, 1)]
    Km = _AND(mf(phL, og), mf(phS, og))
    hok = _OR(moder, normal, tall)
    checks = [Qm, Wm, Gm, Km]

    def failing(i):
        terms = []
        for fail in combinations(range(4), i):
            terms.append(_AND(*[(_NOT(checks[k]) if k in fail else checks[k]) for k in range(4)]))
        return _OR(*terms)
    f1, f2, f3, f4 = failing(1), failing(2), failing(3), failing(4)
    all4 = _AND(*checks)
    filler = _AND(_NOT(Gm), hok, xneg, Wm, Qm, Km)
    sm_mod = _OR(small, moder)
    R = [
        _AND(tiny, _NOT(Km)),
        _AND(small, xpres, _NOT(Km)),
        _AND(tiny, Km),
        _AND(sm_mod, _OR(xpres, f3, f4)),
        _AND(sm_mod, _OR(_NOT(xneg), f2, f3)),
        _AND(hok, _OR(_NOT(Qm), _NOT(xneg), f2)),
        _OR(_AND(hok, Qm, xneg, f2), filler, _AND(small, _NOT(xpres), all4)),
        _AND(Qm, hok, xneg, f1),
        _AND(hok, xneg, all4),
    ]
    return 100.0 * QUAL_DEFUZZ(np.column_stack(R))


# ---------------------------------------------------------------------------
# Window pipeline (§5)
# ---------------------------------------------------------------------------

@dataclass
class WindowResult:
    p: np.ndarray       # global peak index
    letter: np.ndarray  # channel index, -1 = N
    q: np.ndarray
    sp40: float


def run_window(E, dom, X, B, trace, a, b, g_idx, c_idx, prev_sp):
    Ew = E[a:b]
    bands = detect_bands(Ew)
    if len(bands) < 2:
        return None
    gl, _ = gaps(bands)
    s1 = expected_curve(bands.p, gl)
    keep = omission(bands, Ew, X[a:b], s1, 1)
    bands = rebuild(bands.take(np.where(keep)[0]))
    if len(bands) < 2:
        return None
    gl, _ = gaps(bands)
    s2 = expected_curve(bands.p, gl)
    w2 = expected_curve(bands.p, bands.width)
    bands = split_and_insert(bands, Ew, dom[a:b], s2, w2, g_idx, c_idx)
    gl, _ = gaps(bands)
    s3 = expected_curve(bands.p, gl)
    keep = omission(bands, Ew, X[a:b], s3, 2)
    bands = rebuild(bands.take(np.where(keep)[0]))
    if len(bands) == 0:
        return None
    gl, _ = gaps(bands)
    s4 = expected_curve(bands.p, gl)
    w4 = expected_curve(bands.p, bands.width)
    sp = float(np.percentile(gl, 40)) if len(gl) else 10.0
    sp = max(sp, prev_sp) if prev_sp else sp
    q = qualities(bands, Ew, X[a:b], B[a:b], s4, w4, sp)
    gp = bands.p + a
    v = trace[gp]
    ev = E[gp]
    nN = np.sum(v >= N_FRACTION * ev[:, None], axis=1)
    letter = np.where(nN >= 2, -1, dom[gp])
    return WindowResult(gp, letter, q, sp)


def window_starts(m):
    if m <= WIN_LEN:
        return [0]
    out = []
    s = 0
    while True:
        if s + WIN_LEN >= m:
            out.append(m - WIN_LEN)
            break
        out.append(s)
        s += WIN_STEP
    return out


JOIN_ORDER = (4, 5, 3, 6, 2, 7, 1, 8, 0, 9)


def stitch(windows):
    acc_p = np.zeros(0, int); acc_l = np.zeros(0, int); acc_q = np.zeros(0)
    prev = None
    n_fallback = 0
    for w in windows:
        if w is None or len(w.p) == 0:
            continue
        if prev is None:
            acc_p, acc_l, acc_q = w.p.copy(), w.letter.copy(), w.q.copy()
            prev = w
            continue
        best = None
        for k in JOIN_ORDER:
            if k >= len(w.p):
                continue
            j = int(np.argmin(np.abs(prev.p - w.p[k])))
            d = abs(int(prev.p[j]) - int(w.p[k]))
            if d == 0:
                best = (k, j, 0)
                break
            if d < 4 and (best is None or d < best[2]):
                best = (k, j, d)
        if best is not None:
            k, j, _ = best
            cut = acc_p < prev.p[j]
            acc_p = np.concatenate([acc_p[cut], w.p[k:]])
            acc_l = np.concatenate([acc_l[cut], w.letter[k:]])
            acc_q = np.concatenate([acc_q[cut], w.q[k:]])
        else:
            n_fallback += 1
            last = acc_p[-1] if len(acc_p) else -10 ** 9
            sel = w.p > last + 3
            acc_p = np.concatenate([acc_p, w.p[sel]])
            acc_l = np.concatenate([acc_l, w.letter[sel]])
            acc_q = np.concatenate([acc_q, w.q[sel]])
        prev = w
    return acc_p, acc_l, acc_q, n_fallback


# ---------------------------------------------------------------------------
# Good-sequence limits (§15)
# ---------------------------------------------------------------------------

GOOD_W = (7, 9, 11, 13, 15)
GOOD_T = (0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80)
GOOD_WEIGHT = np.array([
    [0.8944, 0.9076, 0.9208, 0.9340, 0.9472, 0.9604, 0.9736],
    [0.8755, 0.8885, 0.9014, 0.9143, 0.9272, 0.9401, 0.9530],
    [0.8567, 0.8693, 0.8819, 0.8946, 0.9072, 0.9199, 0.9325],
    [0.8378, 0.8501, 0.8625, 0.8749, 0.8872, 0.8996, 0.9119],
    [0.8189, 0.8310, 0.8430, 0.8551, 0.8672, 0.8793, 0.8914],
])


def _box(v, w):
    n = len(v)
    before = (w + 1) // 2           # ceil(w/2) calls up to and including the current one
    after = w // 2
    c = np.concatenate([[0.0], np.cumsum(v)])
    i = np.arange(n)
    lo = np.clip(i - before + 1, 0, n)
    hi = np.clip(i + after + 1, 0, n)
    return c[hi] - c[lo]


def good_limits(q):
    q = np.asarray(q, float) / 100.0
    n = len(q)
    if n == 0:
        return 0, 0
    best = (-1, 0, 0)
    for ti, th in enumerate(GOOD_T):
        for wi, w in enumerate(GOOD_W):
            sq = _box(q, w)
            st = _box(np.full(n, th), w)
            good = sq >= st - 1e-12
            # longest run, first if tied
            ba, bb, run_a = -1, -1, None
            for i in range(n + 1):
                if i < n and good[i]:
                    if run_a is None:
                        run_a = i
                elif run_a is not None:
                    if ba < 0 or (i - 1 - run_a) > (bb - ba):
                        ba, bb = run_a, i - 1
                    run_a = None
            if ba < 0:
                continue
            a1, b1 = ba + 1, bb + 1
            score = int(np.floor(GOOD_WEIGHT[wi, ti] * (b1 - a1)))
            if score > best[0]:
                best = (score, a1, b1)
    if best[0] < 0:
        return 0, 0
    return max(best[1] - 2, 0), best[2] - 2


# ---------------------------------------------------------------------------
# Inversion of a flattened trace (§3.3, fallback only)
# ---------------------------------------------------------------------------

def invert_flattening(flat, iterations=100, taps=257, divisor=129.0):
    S = np.clip(flat, 0, None).max(axis=1)
    D = np.ones(len(S))
    for _ in range(iterations):
        D = uniform_filter1d(S * D, taps, mode="constant") * taps / divisor
        mu = D.mean()
        if mu <= 0:
            return flat.copy()
        D = D / mu
    return flat * D[:, None]


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

@dataclass
class CallResult:
    sequence: str
    positions: np.ndarray
    qualities: np.ndarray
    good_start: int
    good_end: int
    info: dict = field(default_factory=dict)


def call_trace(trace, dye_order, start_tp, flattened=False, clip=True, scale_to=0.5) -> CallResult:
    """Call bases on one processed trace (instrument channel order)."""
    tr = np.asarray(trace, float)
    if flattened:
        tr = invert_flattening(tr)
    if clip:
        tr = np.clip(tr, 0, None)
    env = tr.max(axis=1)
    mu = float(np.mean(env)) if len(env) else 0.0
    if scale_to and mu > 0:
        tr = tr * (scale_to / mu)
    m = len(tr)
    order = dye_order.upper()
    g_idx = order.find("G")
    c_idx = order.find("C")
    E, dom, X, B = helpers(tr)
    wins = []
    prev_sp = 0.0
    for s in window_starts(m):
        a, b = s, min(s + WIN_LEN, m)
        w = run_window(E, dom, X, B, tr, a, b, g_idx, c_idx, prev_sp)
        if w is not None:
            prev_sp = w.sp40
        wins.append(w)
    p, let, q, nfb = stitch(wins)
    keep = p <= (m - 1) - END_MARGIN
    p, let, q = p[keep], let[keep], q[keep]
    seq = "".join("N" if c < 0 else dye_order[c] for c in let)
    gs, ge = good_limits(q)
    return CallResult(seq, (p + start_tp + 1).astype(int), q, gs, ge,
                      info=dict(join_fallbacks=nfb, windows=len(wins)))


def call_processed(trace, dye_order, start_tp, meta=None, flattened=False):
    """(sequence, positions, qualities) for a processed trace."""
    r = call_trace(trace, dye_order, start_tp, flattened=flattened)
    return r.sequence, r.positions, r.qualities
