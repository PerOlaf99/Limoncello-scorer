"""
Cimarron-3.12-compatible signal stage for MegaBACE raw traces (clean room).

raw .rsd trace (instrument channel order) -> processed 4-channel trace in the
frame of Cimarron's .esd result file: row i of the output is raw scan
``start_tp - 1 + i`` (START TRACE PROCESSING), ``end_tp - start_tp`` rows.

Pipeline (each step follows the functional specification SPEC_signal.md; the
provenance of every constant is logged in CLEANROOM_rewrite.md):

1. read span (start_tp, end_tp) from the pre-run noise statistics and the
   end-of-run blob;
2. per-channel baseline: moving mean (101) of a moving minimum (101), on the
   span only, edges repeated;
3. colour separation with a 4x4 matrix (instrument table, or estimated from
   the run by direction clustering of band colour vectors);
4. homomorphic (cepstral) blind deconvolution in 2048-row segments with a
   Gaussian low-pass whose width follows the measured band spacing;
5. per-channel whole-row dye-mobility delays, constant in fixed blocks;
6. flattening: divide by a smoothed running-maximum envelope, per-channel
   gain and offset, clip at zero.

``process`` returns both the *internal* trace (after step 5, unflattened,
unclipped) and the *flattened* trace (.esd convention), plus the flattening
divisor and constants so that the flattening can be inverted exactly.
"""

from __future__ import annotations

import re
import struct
from dataclasses import dataclass, field

import numpy as np
from scipy.ndimage import (gaussian_filter1d, maximum_filter1d, minimum_filter1d,
                           uniform_filter1d)

SEG_LEN = 2048
SEG_STEP = 1900

# ---------------------------------------------------------------------------
# Instrument profiles
# ---------------------------------------------------------------------------

# PLAYER1 (MegaBACE 1000, ET terminators, channels T G C A): plate mean of the
# per-run separation matrices stored in the m13 .esd files. Rows = output
# channel, columns = raw channel (SPEC_signal §4).
PLAYER1_SEPARATION = np.array([
    [1.2530, -1.0595, -0.0040, -0.4181],
    [-0.0901, 1.0781, -0.0305, 0.0101],
    [0.2154, 0.0033, 1.1359, -2.0134],
    [-0.5275, 0.1718, -0.0932, 1.3399],
])

# Block edges (output-frame rows) of the dye-mobility delays (SPEC_signal §6).
DELAY_BLOCK_EDGES = (0, 160, 410, 900, 1400, 1900, 2400, 3150, 3900, 4650, 5400, 6150, 6900, 7650)
DELAY_BLOCK_STEP = 750

# PLAYER1 delay table: rows = blocks above, columns = channels in instrument
# order (T, G, C, A); blocks past the table use the last row.
PLAYER1_DELAYS = np.array([
    [-1, 7, 7, 7],
    [1, 8, 7, 7],
    [2, 7, 6, 7],
    [2, 7, 6, 7],
    [3, 7, 6, 7],
    [3, 8, 7, 7],
    [4, 9, 7, 7],
    [5, 9, 8, 8],
    [5, 9, 8, 9],
    [5, 9, 8, 9],
    [5, 9, 8, 9],
    [5, 9, 8, 10],
])

PLAYER1_GAIN = np.array([1.441, 1.319, 1.480, 1.404])
PLAYER1_OFFSET = np.array([-0.0155, 0.0235, -0.0139, 0.1127])
DEFAULT_GAIN = np.full(4, 1.41)
DEFAULT_OFFSET = np.zeros(4)
FLATTEN_SCALE = 0.60

# Estimate the separation matrix from the run even for a known instrument
# (the profile's table is then only the seed).  Measured on m13 even wells:
# end-to-end agreement 0.940 with the estimate vs 0.931 with the PLAYER1
# table and 0.934 with each well's own Cimarron matrix (CLEANROOM_rewrite.md).
ESTIMATE_KNOWN = True

PROFILES = {
    "PLAYER1": dict(dye_order="TGCA", separation=PLAYER1_SEPARATION, delays=PLAYER1_DELAYS,
                    gain=PLAYER1_GAIN, offset=PLAYER1_OFFSET),
}


# ---------------------------------------------------------------------------
# Footer metadata
# ---------------------------------------------------------------------------

_TAG = re.compile(rb"\x05(.)([A-Z][A-Z0-9 /_\-\(\)\.#]{1,40})\x00", re.S)


def read_meta(footer: bytes) -> dict:
    """Text fields of an .rsd footer (tag/value records: a 0x05 type byte, a
    one-byte length, a NUL-terminated key, then a 0x05/0x06 string value)."""
    out: dict = {}
    if not footer:
        return out
    for m in _TAG.finditer(footer):
        if m.group(1)[0] != len(m.group(2)) + 1:
            continue
        j = m.end()
        if j >= len(footer):
            continue
        t = footer[j]
        if t == 0x05 and j + 1 < len(footer):
            n, k = footer[j + 1], j + 2
        elif t == 0x06 and j + 3 <= len(footer):
            n, k = struct.unpack_from("<H", footer, j + 1)[0], j + 3
        else:
            continue
        val = footer[k:k + n].rstrip(b"\0").decode("latin-1")
        out.setdefault(m.group(2).decode(), val)
    return out


_CHANNEL_BASE = re.compile(rb"CHANNEL([1-4])\x00\x01.\x05\x05BASE\x00\x05\x02([ACGTacgt])\x00", re.S)


def dye_order_from_footer(footer: bytes) -> str | None:
    """Base letter of each instrument channel (e.g. "TGCA") from the
    CHANNEL1..4 / BASE records of an .rsd or .esd footer, or None."""
    found: dict = {}
    for m in _CHANNEL_BASE.finditer(footer or b""):
        found.setdefault(int(m.group(1)), m.group(2).decode().upper())
    if len(found) != 4:
        return None
    order = "".join(found[i] for i in range(1, 5))
    return order if sorted(order) == list("ACGT") else None


# ---------------------------------------------------------------------------
# 1. Read span
# ---------------------------------------------------------------------------

def read_span(raw: np.ndarray) -> tuple[int, int]:
    """(S, E) = START / END TRACE PROCESSING (SPEC_signal §2)."""
    raw = np.asarray(raw, dtype=float)
    n = len(raw)
    onset = None
    if n >= 1700:
        pre = raw[200:1200]
        m = np.median(pre, axis=0)
        sd = 1.4826 * np.median(np.abs(pre - m), axis=0)
        sd = np.where(sd > 0, sd, 1e-9)
        s = uniform_filter1d(((raw - m) / sd).max(axis=1), 11, mode="nearest")
        run = 0
        for i in range(1200, n):
            run = run + 1 if s[i] > 8 else 0
            if run >= 20:
                onset = i - 19
                break
    if onset is None:
        onset = min(2001, n // 5) + 9
    S = max(1, onset - 9)
    p = end_blob(raw)
    E = min(n - 248, p + 35) if p is not None else n - 70
    if n > S + 700:
        E = max(E, S + 600)
    E = min(E, n - 1)
    E = max(E, S + 1)
    return int(S), int(E)


def end_blob(raw: np.ndarray, search: int = 1500) -> int | None:
    """Peak row of the end-of-run blob (SPEC_signal §2.4), or None."""
    n = len(raw)
    if n < 3000:
        return None
    u = uniform_filter1d(np.asarray(raw, float).sum(axis=1), 51, mode="nearest")
    p = n - search + int(np.argmax(u[n - search:]))
    if u[p] < 3.0 * np.median(u[n - 2 * search:n - search]):
        return None
    if p >= n - 60:
        return None
    if u[p:].min() > 0.5 * u[p]:
        return None
    return p


# ---------------------------------------------------------------------------
# 2. Baseline
# ---------------------------------------------------------------------------

def remove_baseline(x: np.ndarray, window: int = 101) -> np.ndarray:
    """x - moving mean(window) of moving min(window), per column, edges repeated."""
    b = minimum_filter1d(x, window, axis=0, mode="nearest")
    b = uniform_filter1d(b, window, axis=0, mode="nearest")
    return x - b


# ---------------------------------------------------------------------------
# 3. Colour separation
# ---------------------------------------------------------------------------

def _mixing_from_separation(M: np.ndarray) -> np.ndarray:
    C = np.linalg.inv(M)
    return C / np.diag(C)[None, :]


def estimate_separation(bs: np.ndarray, seed: np.ndarray = PLAYER1_SEPARATION,
                        iterations: int = 10) -> np.ndarray:
    """Separation matrix from the run itself (SPEC_signal §4, run estimate):
    cluster the colour directions of the strongest band picks around the
    seed's dye spectra and re-estimate each spectrum as the median colour
    vector of the purer half of its cluster."""
    seed_mix = _mixing_from_separation(seed)
    tot = gaussian_filter1d(bs.sum(axis=1), 1.0)
    i = np.arange(1, len(tot) - 1)
    pk = i[(tot[i] > tot[i - 1]) & (tot[i] >= tot[i + 1])]
    if len(pk) < 40:
        return seed.copy()
    pk = pk[tot[pk] >= np.percentile(tot[pk], 60)]
    if len(pk) < 40:
        return seed.copy()
    v = np.clip(bs[pk], 1e-9, None)
    u = v / np.linalg.norm(v, axis=1, keepdims=True)
    mix = seed_mix.copy()
    prev = None
    for _ in range(iterations):
        cols = mix / np.linalg.norm(mix, axis=0, keepdims=True)
        cos = u @ cols
        lab = np.argmax(cos, axis=1)
        if prev is not None and np.array_equal(lab, prev):
            break
        prev = lab
        new = mix.copy()
        for k in range(4):
            sel = np.where(lab == k)[0]
            if len(sel) >= 10:
                ck = cos[sel, k]
                sel = sel[ck >= np.median(ck)]
            if len(sel) >= 5:
                col = np.median(v[sel] / v[sel, k:k + 1], axis=0)
                col[k] = 1.0
                new[:, k] = col
        if np.allclose(new, mix):
            mix = new
            break
        mix = new
    try:
        return np.linalg.inv(mix)
    except np.linalg.LinAlgError:
        return seed.copy()


# ---------------------------------------------------------------------------
# 4. Homomorphic deconvolution
# ---------------------------------------------------------------------------

def segment_starts(L: int, N: int = SEG_LEN, step: int = SEG_STEP) -> list[int]:
    """Segment starts 0, step, 2*step, ...; the last one is end-aligned."""
    if L <= N:
        return [0]
    starts = []
    s = 0
    while True:
        if s + N >= L:
            starts.append(L - N)
            break
        starts.append(s)
        s += step
    return starts


def segment_ownership(starts: list[int], L: int, N: int = SEG_LEN) -> list[tuple[int, int]]:
    """Rows [a, b) each segment supplies: from the middle of its overlap with
    the previous segment to the middle of its overlap with the next."""
    out = []
    for k, s in enumerate(starts):
        a = 0 if k == 0 else (s + starts[k - 1] + N) // 2
        b = L if k == len(starts) - 1 else (starts[k + 1] + s + N) // 2
        out.append((a, min(b, L)))
    return out


def _lifter(N: int) -> np.ndarray:
    m = np.arange(N // 2 + 1, dtype=float)
    w = np.ones_like(m)
    w[m <= 11] = 0.0
    r = (m >= 12) & (m <= 22)
    w[r] = 0.5 * (1.0 - np.cos(np.pi * (m[r] - 12) / 10.0))
    full = np.empty(N)
    full[:N // 2 + 1] = w
    full[N // 2 + 1:] = w[1:N // 2][::-1]
    return full


def _gaussian_lowpass(N: int, W: float) -> np.ndarray:
    k = np.arange(N)
    omega = 2 * np.pi * np.where(k <= N // 2, k, k - N) / N
    sigma = 2 * np.pi * (W - 1) / N
    return np.sqrt(2 * np.pi) / sigma * np.exp(-omega ** 2 / (2 * sigma ** 2))


def homomorphic_deconvolve(seg: np.ndarray, W: float) -> np.ndarray:
    """SPEC_signal §5.3 on each column of an (N, c) segment."""
    N = seg.shape[0]
    X = np.fft.fft(seg, axis=0)
    X = np.where(X == 0, np.finfo(float).eps, X)
    logmag = np.log(np.abs(X))
    phase = np.angle(X)
    q = np.real(np.fft.ifft(logmag, axis=0))
    q *= _lifter(N)[:, None]
    lifted = np.real(np.fft.fft(q, axis=0))
    Z = np.exp(lifted + 1j * phase) * _gaussian_lowpass(N, W)[:, None]
    return np.real(np.fft.ifft(Z, axis=0))


# Width of the spacing-estimation pass.  SPEC_signal §5.4 uses 135; with 135
# the band count misses merged bands late in the read and the spacing comes
# out 23 % long in the last segment (7 % in the first); 170 brings the bias to
# -5..+6 % in every segment (measured against Cimarron's call spacing).
FBW_FIRST = 170.0
FBW_CONSTANT = 1210.8
BAND_THRESH = 0.18
FBW_PASSES = 1        # spacing re-estimations before the final pass
FBW_SEG_W = None      # optional f(segment index, n segments, W, spacing) -> W (development hook)


def band_spacing(first_pass: np.ndarray, n_real: int, margin: int = 100) -> float:
    """Mean gap between band positions of a first-pass segment (§5.4.2-3)."""
    if n_real < 400:
        return 10.0
    e = np.clip(first_pass[:n_real], 0, None).max(axis=1)
    lam = gaussian_filter1d(maximum_filter1d(e, 15, mode="nearest"), 100.0, mode="nearest")
    r = e / np.maximum(lam, 1e-12)
    i = np.arange(max(1, margin), n_real - margin)
    if len(i) == 0:
        return 10.0
    pos = i[(r[i] > r[i - 1]) & (r[i] >= r[i + 1]) & (r[i] > BAND_THRESH)]
    if len(pos) < 10:
        return 10.0
    return float(np.mean(np.diff(pos)))


def deconvolve(sep: np.ndarray) -> tuple[np.ndarray, list[float]]:
    """Segment-wise deconvolution of the separated span (§5.2-5.4)."""
    L = sep.shape[0]
    N = SEG_LEN
    starts = segment_starts(L, N)
    owners = segment_ownership(starts, L, N)
    out = np.zeros_like(sep)
    widths = []
    for kseg, (s, (a, b)) in enumerate(zip(starts, owners)):
        seg = sep[s:s + N]
        n_real = len(seg)
        if n_real < N:
            seg = np.vstack([seg, np.repeat(seg[-1:], N - n_real, axis=0)])
        W = FBW_FIRST
        for _ in range(FBW_PASSES):
            sp = band_spacing(homomorphic_deconvolve(seg, W), n_real)
            W = float(np.clip(FBW_CONSTANT / max(sp, 3.0), 40.0, 400.0))
        if FBW_SEG_W is not None:
            W = float(FBW_SEG_W(kseg, len(starts), W, sp))
        widths.append(W)
        res = homomorphic_deconvolve(seg, W)
        out[a:b] = res[a - s:b - s]
    return out, widths


# ---------------------------------------------------------------------------
# 5. Block delays
# ---------------------------------------------------------------------------

def block_edges(L: int) -> list[int]:
    edges = [e for e in DELAY_BLOCK_EDGES if e < L]
    nxt = DELAY_BLOCK_EDGES[-1] + DELAY_BLOCK_STEP
    while nxt < L:
        edges.append(nxt)
        nxt += DELAY_BLOCK_STEP
    return edges


def apply_block_delays(x: np.ndarray, table: np.ndarray) -> np.ndarray:
    """delayed_c(t) = x_c(t - d[block(t), c]), indices clamped to the span."""
    L = x.shape[0]
    edges = block_edges(L) + [L]
    out = np.empty_like(x)
    for j in range(len(edges) - 1):
        a, b = edges[j], edges[j + 1]
        d = table[min(j, len(table) - 1)]
        t = np.arange(a, b)
        for c in range(x.shape[1]):
            out[a:b, c] = x[np.clip(t - d[c], 0, L - 1), c]
    return out


# ---------------------------------------------------------------------------
# 6. Flattening
# ---------------------------------------------------------------------------

def flatten_divisor(internal: np.ndarray) -> np.ndarray:
    """lambda(t): Gaussian(100)-smoothed 15-row running max of the envelope."""
    e = np.clip(internal, 0, None).max(axis=1)
    lam = gaussian_filter1d(maximum_filter1d(e, 15, mode="nearest"), 100.0, mode="nearest")
    return np.maximum(lam, 1e-6 * max(float(lam.max()), 1e-300))


def flatten(internal: np.ndarray, lam: np.ndarray, gain: np.ndarray, offset: np.ndarray) -> np.ndarray:
    F = np.clip(FLATTEN_SCALE * internal / lam[:, None], 0, None)
    return np.clip(gain[None, :] * F - offset[None, :], 0, None)


def unflatten(flattened: np.ndarray, lam: np.ndarray, gain: np.ndarray, offset: np.ndarray) -> np.ndarray:
    """Exact inverse of `flatten` wherever neither clip was active (internal
    value positive and gain*F - offset positive); clipped entries come back
    as zero."""
    F = np.where(flattened > 0, (flattened + offset[None, :]) / gain[None, :], 0.0)
    return np.clip(F, 0, None) * lam[:, None] / FLATTEN_SCALE


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

@dataclass
class SignalResult:
    trace: np.ndarray          # flattened, .esd convention, instrument channel order
    internal: np.ndarray       # before flattening/clipping, same frame
    start_tp: int
    end_tp: int
    divisor: np.ndarray        # lambda(t)
    gain: np.ndarray
    offset: np.ndarray
    separated: np.ndarray      # baseline-subtracted, separated span (raw amplitude)
    estimated: bool            # True when the instrument was not known
    info: dict = field(default_factory=dict)

    @property
    def flattened(self) -> np.ndarray:
        return self.trace

    @property
    def caller_trace(self) -> np.ndarray:
        """The flattened trace with the flattening divisor multiplied back
        (exactly: flattened * lambda / 0.60).  Keeps the channel gains, the
        offsets that suppress small cross-talk bumps, and the clipping of the
        .esd trace, but restores the internal amplitude course.  This is what
        the calling stage is fed."""
        return self.trace * (self.divisor / FLATTEN_SCALE)[:, None]


def process(raw: np.ndarray, dye_order: str | None = None, meta: dict | None = None,
            span: tuple[int, int] | None = None, separation: np.ndarray | None = None) -> SignalResult:
    """Raw (n, 4) trace in instrument channel order -> SignalResult."""
    raw = np.asarray(raw, dtype=float)
    meta = meta or {}
    machine = str(meta.get("MACHINE ID", "")).strip().upper()
    prof = PROFILES.get(machine)
    estimated = prof is None
    S, E = span if span is not None else read_span(raw)
    x = raw[S - 1:E - 1]
    bs = remove_baseline(x)
    if separation is not None:
        M, msrc = np.asarray(separation, float), "given"
    elif prof is not None and not ESTIMATE_KNOWN:
        M, msrc = prof["separation"], "table"
    else:
        seed = prof["separation"] if prof is not None else PLAYER1_SEPARATION
        M, msrc = estimate_separation(bs, seed=seed), "estimated"
    sep = bs @ M.T
    dec, widths = deconvolve(sep)
    delays = (prof or PROFILES["PLAYER1"])["delays"]
    internal = apply_block_delays(dec, delays)
    gain = (prof or PROFILES["PLAYER1"])["gain"]
    offset = (prof or PROFILES["PLAYER1"])["offset"]
    lam = flatten_divisor(internal)
    flat = flatten(internal, lam, gain, offset)
    info = dict(machine=machine, separation=M, separation_source=msrc, widths=widths,
                delays_source="table" if prof else "PLAYER1 default")
    return SignalResult(trace=flat, internal=internal, start_tp=S, end_tp=E, divisor=lam,
                        gain=np.asarray(gain, float), offset=np.asarray(offset, float),
                        separated=sep, estimated=estimated, info=info)


def process_rsd(path: str) -> tuple[SignalResult, dict, str | None]:
    """Read an .rsd file and process it. Returns (result, meta, dye_order)."""
    from .rsd_io import read_rsd
    r = read_rsd(path)
    meta = read_meta(r.footer_raw)
    order = dye_order_from_footer(r.footer_raw)
    return process(r.trace, order, meta), meta, order
