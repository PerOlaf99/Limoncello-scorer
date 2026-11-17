"""Cycling temperature capillary electrophoresis (CTCE): find the peaks of
interest by how they move between runs at different temperatures, and
quantify them.

Background (Ekstrøm, Bjørheim & Thilly, BMC Genetics 8:54, 2007; Ekstrøm,
Warren & Thilly, Electrophoresis 33:1162, 2012): a fragment with a low-melting
domain beside a GC clamp is run while the capillary temperature cycles.  The
partly melted molecules move more slowly, and how much depends on sequence
and duplex state, so the two homoduplexes and the two heteroduplexes of a
mixed sample separate.  Raising the temperature moves these duplex peaks
later; molecules without a melting domain in range (ssDNA, primers, the
labelled oligo) keep their relative position.  Below the window all duplexes
still co-migrate (one peak), in the window they resolve (four peaks for two
variants, up to n² for n), above it they are fully melted, stop moving apart
and merge again.

So, per well, from two or more runs of the same plate at different
temperatures (read from each .rsd, never assumed from the run number):

1. one channel only -- the sample channel (strongest signal, or the one
   given); an internal standard on another channel, when present, is used
   only to line the runs up;
2. peaks of that channel, with valley-bounded areas and a shape measure
   (duplex peaks have a sharp front and a tail; ssDNA is symmetric);
3. the runs lined up on a common anchor (internal standard bands, else the
   large early peak);
4. between each pair of neighbouring temperatures, every pair of peaks of
   similar height votes for a shift; the strongest non-zero shift that a
   group of peaks shares is the melting group, its members the peaks of
   interest; peaks that pair with no shift are static (ssDNA / primer);
5. per run, how many melting peaks are resolved (the regime), the run with
   the most of them as the measuring run, and areas there -- a satellite
   (stutter / +A) that moves with its parent is reported and, by default,
   added to its parent;
6. fractions: each melting peak's share, and for four peaks (two
   homoduplexes, then two heteroduplexes, in migration order) the minor
   variant fraction by mass action, p = (A_hom + A_het / 2) / A_total.

Headless and free of Tk, so it can be scripted and tested.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
from scipy.ndimage import minimum_filter1d, uniform_filter1d
from scipy.signal import find_peaks

# --------------------------------------------------------------------------- #
# settings
# --------------------------------------------------------------------------- #
MIN_REL_PROMINENCE = 0.03     # peaks: prominence >= 3 % of the channel maximum
SHAPE_LEVEL = 0.25            # front/tail widths measured at 25 % of the height
FOOT_LEVEL = 0.005            # an outer peak's area runs out to 0.5 % of its height
STATIC_WINDOW = (-60, 10)     # a peak that moved between -60 and +10 scans per step is static
                              # (ssDNA drifts slightly EARLIER as it warms; melting peaks move later)
PAIR_HEIGHT_RATIO = (0.4, 2.5)  # a peak and its partner at the next temperature
MIN_PAIR_HEIGHT = 0.08        # ... both at least 8 % of the well's maximum
ANCHOR_REL_HEIGHT = 0.25      # anchor = earliest peak above 25 % of the maximum
SATELLITE_GAP = (6, 28)       # a satellite sits this many scans after its parent
SATELLITE_MAX_RATIO = 0.9     # ... and is lower than it
NO_PRODUCT_RATIO = 0.15       # melting signal below 15 % of the plate median: no product
DUPLEX_ASYM = 1.3             # one temperature: tail/front >= this is a duplex ...
SSDNA_ASYM = 1.2              # ... and below this a symmetric (static) peak
TEMPLATE_SHIFT = (40, 0.15)   # a well's pattern may sit max(40 scans, 15 % of the offset) off ...
TEMPLATE_TOL = (5, 0.02)      # ... but its peaks keep their spacing within max(5, 2 % of offset)
SHAPE_MIN_HEIGHT = 0.15       # one temperature: duplex peaks at least 15 % of the maximum
SAME_TEMPERATURE = 0.5        # runs within 0.5 °C count as the same temperature


# --------------------------------------------------------------------------- #
# data
# --------------------------------------------------------------------------- #
@dataclass
class Peak:
    pos: int                  # apex, scans (aligned when inside a WellResult)
    height: float
    left: int                 # valley before
    right: int                # valley after
    area: float               # above the straight line between the two valleys
    front: int                # width before the apex at SHAPE_LEVEL of the height
    tail: int                 # width after the apex
    role: str = "unassigned"  # melting / static / satellite / anchor / other
    parent: Optional[int] = None   # index of the parent peak, for a satellite
    raw_pos: int = -1         # apex in the run's own scans (pos may be aligned)
    label: str = ""           # template peak name (P1, P2, ...) when user-guided

    @property
    def asym(self) -> float:
        """Tail / front width: ~1 for a symmetric (ssDNA) peak, > 1.3 for the
        sharp-front, long-tail shape of a melting duplex."""
        return self.tail / max(self.front, 1)


@dataclass
class RunTrace:
    path: Path
    run: str
    plate: str
    well: str
    temperature: Optional[float]
    start: str
    channel: int              # physical channel analysed, 1-based
    y: np.ndarray             # baseline-corrected sample channel
    offset: float = 0.0       # alignment: aligned = scale * scan + offset
    scale: float = 1.0
    anchor: Optional[float] = None   # the large early peak, raw scans
    peaks: List[Peak] = field(default_factory=list)


@dataclass
class WellResult:
    plate: str
    well: str
    runs: List[RunTrace]
    shifts: List[Optional[float]]       # melting-group shift between neighbouring temperatures
    regimes: List[str]                  # per run: no product / single / resolved n / merged
    best: Optional[int]                 # index of the measuring run
    areas: List[float] = field(default_factory=list)       # melting peaks at the best run
    fractions: List[float] = field(default_factory=list)
    fraction_sd: List[float] = field(default_factory=list)  # over replicate runs (same temperature)
    method: str = "change"              # change (temperatures) / shape (one temperature) / template
    minor_fraction: Optional[float] = None
    flags: List[str] = field(default_factory=list)

    def summary(self) -> Dict:
        best = self.runs[self.best] if self.best is not None else None
        return {
            "plate": self.plate, "well": self.well,
            "temperatures": "/".join("?" if r.temperature is None else f"{r.temperature:g}"
                                     for r in self.runs),
            "shifts": "/".join("" if s is None else f"{s:+.0f}" for s in self.shifts),
            "regimes": "/".join(self.regimes),
            "measured_at": "" if best is None else f"{best.temperature:g}",
            "n_peaks": len(self.areas),
            "areas": "/".join(f"{a:.0f}" for a in self.areas),
            "fractions": "/".join(f"{f:.3f}" for f in self.fractions),
            "fraction_sd": "/".join(f"{s:.3f}" for s in self.fraction_sd),
            "method": self.method,
            "minor_fraction": "" if self.minor_fraction is None else round(self.minor_fraction, 4),
            "flags": ",".join(self.flags),
        }


# --------------------------------------------------------------------------- #
# one trace
# --------------------------------------------------------------------------- #
_RUN_RE = re.compile(r"[_\- ]?run\s*\d+$", re.I)


def plate_of(run_name: str) -> str:
    """Runs of one plate differ only by a trailing RunNN."""
    return _RUN_RE.sub("", run_name)


def sample_channel(acgt: np.ndarray, order: str) -> int:
    """Physical channel (1-based) with the strongest signal."""
    x = np.asarray(acgt, float) - np.median(np.asarray(acgt, float)[:300], axis=0)
    phys = [x[:, "ACGT".index(b)] for b in order[:4]]
    return int(np.argmax([np.percentile(c, 99.9) for c in phys])) + 1


def channel_trace(acgt: np.ndarray, order: str, channel: int,
                  baseline_window: int = 151) -> np.ndarray:
    """One physical channel, baseline removed (running minimum, smoothed)."""
    y = np.asarray(acgt, float)[:, "ACGT".index(order[channel - 1])]
    base = uniform_filter1d(minimum_filter1d(y, baseline_window), baseline_window)
    return y - base


def _width(y: np.ndarray, i: int, level: float, step: int) -> int:
    lvl = y[i] * level
    j = i
    while 0 < j < len(y) - 1 and y[j] > lvl:
        j += step
    return abs(j - i)


def detect_peaks(y: np.ndarray, min_rel_prominence: float = MIN_REL_PROMINENCE) -> List[Peak]:
    """Peaks with valley-bounded areas and front/tail widths."""
    y = np.asarray(y, float)
    if y.size < 10 or y.max() <= 0:
        return []
    idx, _props = find_peaks(y, prominence=min_rel_prominence * y.max(), distance=4)
    if idx.size == 0:
        return []
    # valleys: the lowest point between neighbouring apexes (trace ends outside)
    bounds = [0]
    for a, b in zip(idx[:-1], idx[1:]):
        bounds.append(int(a + np.argmin(y[a:b + 1])))
    bounds.append(len(y) - 1)
    out = []
    for k, i in enumerate(idx):
        lo, hi = bounds[k], bounds[k + 1]
        if k == 0:                       # first peak: walk down to its own foot
            lo = i - _width(y, i, FOOT_LEVEL, -1)
        if k == len(idx) - 1:
            hi = i + _width(y, i, FOOT_LEVEL, +1)
        lo, hi = max(0, lo), min(len(y) - 1, hi)
        seg = y[lo:hi + 1]
        line = np.linspace(y[lo], y[hi], hi - lo + 1)
        area = float(np.clip(seg - line, 0, None).sum())
        out.append(Peak(pos=int(i), height=float(y[i]), left=lo, right=hi, area=area,
                        front=_width(y, i, SHAPE_LEVEL, -1), tail=_width(y, i, SHAPE_LEVEL, +1),
                        raw_pos=int(i)))
    return out


def load_run_trace(path, channel: Optional[int] = None) -> RunTrace:
    """A well's sample channel from one run, with the run's temperature."""
    import analyzer_core as ac
    path = Path(path)
    doc = ac.load_trace(path)
    order = ac.declared_dye_order(getattr(doc, "footer_fields", None)) or "TGCA"
    rp = getattr(doc, "run_params", None) or {}
    try:
        temp = float(str(rp.get("Temperature", "")).strip())
    except ValueError:
        temp = None
    ch = channel or sample_channel(doc.acgt, order)
    y = channel_trace(doc.acgt, order, ch)
    rt = RunTrace(path=path, run=path.parent.name, plate=plate_of(path.parent.name),
                  well=path.stem, temperature=temp, start=str(rp.get("Run start time", "")),
                  channel=ch, y=y)
    rt.peaks = detect_peaks(y)
    rt.anchor = anchor_scan(rt)          # before any alignment touches the peaks
    rt._doc = doc                          # kept for the internal-standard path
    rt._order = order
    return rt


# --------------------------------------------------------------------------- #
# lining runs up
# --------------------------------------------------------------------------- #
def anchor_scan(rt: RunTrace) -> Optional[float]:
    """The early peak (labelled primer / oligo): the earliest peak above a
    quarter of the channel maximum.  Not "the first big peak": in high-yield
    wells the primer peak is far smaller than the duplex peaks."""
    if not rt.peaks:
        return None
    top = max(p.height for p in rt.peaks)
    big = [p for p in rt.peaks if p.height >= ANCHOR_REL_HEIGHT * top]
    return float(min(big, key=lambda p: p.pos).pos) if big else None


def internal_standard_scans(rt: RunTrace, is_channel: int) -> Optional[List[float]]:
    """Bands of an internal standard on another channel (genotyping's
    equimolar quartet search), or None."""
    try:
        import genotyping
        y = channel_trace(rt._doc.acgt, rt._order, is_channel)
        hit = genotyping.find_is_quartet(y)
    except Exception:                                  # noqa: BLE001
        return None
    return [float(s) for s in hit[0]] if hit else None


def align_runs(runs: List[RunTrace], is_channel: Optional[int] = None) -> str:
    """Put every run on the first run's scale: offset and stretch from the
    internal standard when it is found in every run, else offset from the
    large early peak.  Returns which method was used."""
    if is_channel:
        bands = [internal_standard_scans(r, is_channel) for r in runs]
        if all(b and len(b) == len(bands[0]) and len(b) >= 2 for b in bands):
            ref = np.asarray(bands[0])
            for r, b in zip(runs, bands):
                b = np.asarray(b)
                if len(b) >= 3:
                    m, c = np.polyfit(b, ref, 1)
                else:
                    m, c = 1.0, float(np.median(ref - b))
                r.scale, r.offset = float(m), float(c)
            _apply_alignment(runs)
            return "internal standard"
    anchors = [r.anchor for r in runs]
    if anchors[0] is None or any(a is None for a in anchors):
        return "none"
    for r, a in zip(runs, anchors):
        r.scale, r.offset = 1.0, float(anchors[0] - a)
    _apply_alignment(runs)
    return "early peak"


def _apply_alignment(runs: List[RunTrace]) -> None:
    for r in runs:
        for p in r.peaks:
            p.pos = int(round(r.scale * p.pos + r.offset))
            p.left = int(round(r.scale * p.left + r.offset))
            p.right = int(round(r.scale * p.right + r.offset))


# --------------------------------------------------------------------------- #
# what moves
# --------------------------------------------------------------------------- #
def _similar(a: Peak, b: Peak) -> float:
    lo, hi = PAIR_HEIGHT_RATIO
    r = b.height / max(a.height, 1e-9)
    return 0.0 if not (lo <= r <= hi) else min(r, 1.0 / r)


def match_step(pa: Sequence[Peak], pb: Sequence[Peak], ymax: float, warmer: bool = True):
    """Peaks of a colder run (pa) against the next warmer run (pb).

    1. static: a peak pairs with one of similar height that moved by a
       STATIC_WINDOW amount (ssDNA, primer -- they hold their relative
       position, at most drifting a little earlier);
    2. melting: the rest are followed one by one in migration order (duplex
       species do not overtake each other), each pair moving later than the
       static window -- by its own amount, since every species has its own
       melting curve.
    Returns (static_pairs, melting_pairs); with warmer=False the direction is
    reversed (pb colder than pa)."""
    s_lo, s_hi = STATIC_WINDOW
    if not warmer:
        s_lo, s_hi = -s_hi, -s_lo
    big_a = [p for p in pa if p.height >= MIN_PAIR_HEIGHT * ymax]
    big_b = [p for p in pb if p.height >= MIN_PAIR_HEIGHT * ymax]
    cand = sorted(((_similar(x, y) - abs(y.pos - x.pos) / 200.0, x, y)
                   for x in big_a for y in big_b
                   if s_lo <= y.pos - x.pos <= s_hi and _similar(x, y) > 0),
                  key=lambda c: -c[0])
    used, static = set(), []
    for _s, x, y in cand:
        if id(x) in used or id(y) in used:
            continue
        used.update((id(x), id(y)))
        static.append((x, y))
    ra = sorted((p for p in big_a if id(p) not in used), key=lambda p: p.pos)
    rb = sorted((p for p in big_b if id(p) not in used), key=lambda p: p.pos)
    # order-preserving matching (longest-common-subsequence style), each pair
    # moving beyond the static window in the warming direction
    n, m = len(ra), len(rb)
    best = [[0.0] * (m + 1) for _ in range(n + 1)]
    take = [[False] * (m + 1) for _ in range(n + 1)]
    for i in range(n - 1, -1, -1):
        for j in range(m - 1, -1, -1):
            d = rb[j].pos - ra[i].pos
            ok = d > s_hi if warmer else d < s_lo
            sim = _similar(ra[i], rb[j]) if ok else 0.0
            choices = [best[i + 1][j], best[i][j + 1]]
            if sim > 0:
                choices.append(sim + best[i + 1][j + 1])
            best[i][j] = max(choices)
            take[i][j] = sim > 0 and best[i][j] == sim + best[i + 1][j + 1]
    melt, i, j = [], 0, 0
    while i < n and j < m:
        if take[i][j]:
            melt.append((ra[i], rb[j]))
            i, j = i + 1, j + 1
        elif best[i + 1][j] >= best[i][j + 1]:
            i += 1
        else:
            j += 1
    return static, melt


def _mark_satellites(peaks: List[Peak]) -> None:
    """A smaller peak a few scans after a melting peak is its satellite (Taq
    stutter / +A).  Decided within the run, so the last run's satellites are
    found too; a satellite that was itself tracked keeps role "satellite"."""
    melt = sorted((k for k, p in enumerate(peaks) if p.role == "melting"),
                  key=lambda k: peaks[k].pos)
    for k, p in enumerate(peaks):
        if p.role not in ("melting", "unassigned") or p.label:
            continue                     # a clicked peak is a peak of interest
        for j in melt:
            q = peaks[j]
            if j == k or q.role != "melting":
                continue
            gap = p.pos - q.pos
            if (SATELLITE_GAP[0] <= gap <= SATELLITE_GAP[1]
                    and p.height <= SATELLITE_MAX_RATIO * q.height):
                p.role, p.parent = "satellite", j
                break


# --------------------------------------------------------------------------- #
# user guidance: a template from picked peaks
# --------------------------------------------------------------------------- #
@dataclass
class Template:
    """Peaks of interest the user clicked, as offsets from the early peak,
    per run temperature: {temperature: [offset, ...]} in migration order."""
    offsets: Dict[float, List[float]] = field(default_factory=dict)
    n_wells: int = 0

    def at(self, temperature: Optional[float]) -> Optional[List[float]]:
        if temperature is None:
            return None
        for t, offs in self.offsets.items():
            if abs(t - temperature) <= SAME_TEMPERATURE:
                return offs
        return None


def template_from_picks(picks: Dict[str, Sequence[float]], channel: Optional[int] = None) -> Template:
    """{trace path: [picked scans]} -> Template.  Each pick is taken relative
    to its own run's early peak; wells picked at one temperature are pooled
    by peak order (median offset), so a few clicked wells teach the plate."""
    per_temp: Dict[float, List[List[float]]] = {}
    for path, scans in picks.items():
        scans = sorted(float(s) for s in scans)
        if not scans:
            continue
        rt = load_run_trace(path, channel)
        if rt.anchor is None or rt.temperature is None:
            continue
        key = next((t for t in per_temp if abs(t - rt.temperature) <= SAME_TEMPERATURE),
                   rt.temperature)
        per_temp.setdefault(key, []).append([s - rt.anchor for s in scans])
    tpl = Template()
    for temp, lists in per_temp.items():
        n = max(len(x) for x in lists)
        full = [x for x in lists if len(x) == n]          # wells with every peak clicked
        tpl.offsets[temp] = [float(np.median([x[i] for x in full])) for i in range(n)]
        tpl.n_wells += len(full)
    return tpl


def _apply_template(run: RunTrace, offsets: Sequence[float]) -> int:
    """Mark the run's peaks matching the template as the peaks of interest
    (P1, P2, ...); every other melting peak there is demoted.

    The template is matched as a pattern: the whole set may sit somewhat off
    in this well (one shift), but its peaks keep their spacing to each other.
    Returns how many template peaks were found."""
    if run.anchor is None or not offsets:
        return 0
    for p in run.peaks:
        if p.role in ("melting", "satellite"):
            p.role, p.parent = "other", None
    pool = [p for p in run.peaks if p.role != "anchor" and p.raw_pos > run.anchor + 20]
    first = run.anchor + offsets[0]
    reach = max(TEMPLATE_SHIFT[0], TEMPLATE_SHIFT[1] * abs(offsets[0]))
    best = (0, 0.0, 0.0, [])
    for cand in pool:
        s = cand.raw_pos - first
        if abs(s) > reach:
            continue
        chosen, used = [], set()
        for off in offsets:
            target = run.anchor + off + s
            tol = max(TEMPLATE_TOL[0], TEMPLATE_TOL[1] * abs(off))
            near = [p for p in pool if id(p) not in used and abs(p.raw_pos - target) <= tol]
            q = min(near, key=lambda q: (abs(q.raw_pos - target), -q.height)) if near else None
            chosen.append(q)
            if q is not None:
                used.add(id(q))
        n = sum(q is not None for q in chosen)
        h = sum(q.height for q in chosen if q is not None)
        score = (n, h, -abs(s))
        if score > best[:3]:
            best = (n, h, -abs(s), chosen)
    for i, q in enumerate(best[3]):
        if q is not None:
            q.role, q.label = "melting", f"P{i + 1}"
    return best[0]


def _classify_by_shape(run: RunTrace, ymax: float) -> None:
    """One temperature: a duplex has a sharp front and a tail, ssDNA is symmetric."""
    for p in run.peaks:
        if p.role != "unassigned" or p.height < SHAPE_MIN_HEIGHT * ymax:
            continue
        if run.anchor is not None and p.raw_pos <= run.anchor + 60:
            continue
        if p.asym >= DUPLEX_ASYM:
            p.role = "melting"
        elif p.asym < SSDNA_ASYM:
            p.role = "static"


# --------------------------------------------------------------------------- #
# one well
# --------------------------------------------------------------------------- #
def analyse_well(paths: Sequence, channel: Optional[int] = None,
                 is_channel: Optional[int] = None,
                 merge_satellites: bool = True,
                 template: Optional[Template] = None) -> WellResult:
    """CTCE analysis of one well over its runs at different temperatures."""
    runs = [load_run_trace(p, channel) for p in paths]
    # order by temperature (a hotter run done first must not flip the direction)
    runs.sort(key=lambda r: (r.temperature is None, r.temperature or 0.0, r.start))
    res = WellResult(plate=runs[0].plate, well=runs[0].well, runs=runs,
                     shifts=[], regimes=[], best=None)
    temps = [r.temperature for r in runs]
    distinct = sorted({t for t in temps if t is not None})
    one_temp = len(distinct) < 2 or (distinct[-1] - distinct[0]) <= SAME_TEMPERATURE
    if one_temp:
        res.method = "shape"
        res.flags.append("one temperature: peaks classified by shape"
                         + (f", {len(runs)} replicate runs" if len(runs) > 1 else ""))
    method = align_runs(runs, is_channel)
    if method == "none":
        res.flags.append("runs could not be lined up")
    ymax = max(float(r.y.max()) for r in runs)
    for r in runs:
        if r.anchor is None:
            continue
        aligned = r.scale * r.anchor + r.offset
        for p in r.peaks:
            if abs(p.pos - aligned) <= 2:
                p.role = "anchor"
    # 4. static peaks and melting peaks between neighbouring temperatures
    for k in range(0 if one_temp else len(runs) - 1):
        pa = [p for p in runs[k].peaks if p.role != "anchor"]
        pb = [p for p in runs[k + 1].peaks if p.role != "anchor"]
        warmer = not (temps[k] is not None and temps[k + 1] is not None and temps[k + 1] < temps[k])
        static, melt = match_step(pa, pb, ymax, warmer)
        for a, b in static:
            for p in (a, b):
                if p.role == "unassigned":
                    p.role = "static"
        if not melt:
            res.shifts.append(None)
            continue
        res.shifts.append(float(np.median([b.pos - a.pos for a, b in melt])))
        for a, b in melt:
            a.role = b.role = "melting"
    if one_temp:
        res.shifts = [None] * max(0, len(runs) - 1)
        for r in runs:
            _classify_by_shape(r, ymax)
    # user guidance wins where it applies: the clicked peaks, by position
    if template is not None:
        guided = 0
        for r in runs:
            offs = template.at(r.temperature)
            if offs:
                if _apply_template(r, offs) < len(offs):
                    res.flags.append(f"template: {len(offs) - sum(p.label != '' for p in r.peaks)} "
                                     f"peak(s) not found at {r.temperature:g} °C")
                guided += 1
        if guided:
            res.method = "template"
    for r in runs:
        _mark_satellites(r.peaks)
    # 5. regime per run, and the measuring run
    melt_heights = []
    for r in runs:
        core = [p for p in r.peaks if p.role == "melting"]
        melt_heights.append(max([p.height for p in core], default=0.0))
        if not core:
            res.regimes.append("none")
        elif len(core) == 1:
            res.regimes.append("single")
        else:
            res.regimes.append(f"resolved {len(core)}")
    counts = [sum(p.role == "melting" for p in r.peaks) for r in runs]
    if max(counts, default=0) == 0:
        res.flags.append("no melting peaks")
        return res
    # most resolved; ties to the run in the middle of the series
    order = sorted(range(len(runs)), key=lambda k: (-counts[k], abs(k - (len(runs) - 1) / 2)))
    res.best = order[0]
    if counts[-1] < max(counts) and len(runs) > 2 and counts[-1] <= 1 and res.best != len(runs) - 1:
        res.flags.append("merged again at the highest temperature")
    # 6. areas and fractions at the measuring run
    best = runs[res.best]
    core = [p for p in best.peaks if p.role == "melting"]
    sat = [p for p in best.peaks if p.role == "satellite"]
    areas = {id(p): p.area for p in core}
    if merge_satellites:
        for s in sat:
            parent = best.peaks[s.parent]
            areas[id(parent)] = areas.get(id(parent), 0.0) + s.area
    core.sort(key=lambda p: p.pos)
    res.areas = [areas[id(p)] for p in core]
    total = sum(res.areas)
    res.fractions = [a / total for a in res.areas] if total > 0 else []
    if len(core) == 4 and total > 0:
        hom1, hom2, het1, het2 = res.areas
        p = (hom1 + (het1 + het2) / 2.0) / total
        res.minor_fraction = min(p, 1.0 - p)
    # replicates: runs at the measuring temperature with the same peaks
    reps = [r for r in runs if r.temperature is not None and best.temperature is not None
            and abs(r.temperature - best.temperature) <= SAME_TEMPERATURE]
    if len(reps) > 1:
        fr = []
        for r in reps:
            a = _run_areas(r, merge_satellites)
            if len(a) == len(res.areas) and sum(a) > 0:
                fr.append([x / sum(a) for x in a])
        if len(fr) > 1:
            arr = np.asarray(fr)
            res.fractions = [float(x) for x in arr.mean(axis=0)]
            res.fraction_sd = [float(x) for x in arr.std(axis=0, ddof=1)]
    if sat:
        res.flags.append(f"{len(sat)} satellite peak(s) {'added to' if merge_satellites else 'kept apart from'} their parent")
    return res


def _run_areas(run: RunTrace, merge_satellites: bool) -> List[float]:
    """Areas of a run's melting peaks in migration order (satellites merged)."""
    core = sorted((p for p in run.peaks if p.role == "melting"), key=lambda p: p.pos)
    areas = {id(p): p.area for p in core}
    if merge_satellites:
        for s in (p for p in run.peaks if p.role == "satellite"):
            par = run.peaks[s.parent]
            if id(par) in areas:
                areas[id(par)] += s.area
    return [areas[id(p)] for p in core]


# --------------------------------------------------------------------------- #
# a plate
# --------------------------------------------------------------------------- #
def flag_no_product(results: Sequence[WellResult]) -> None:
    """Flag wells whose signal (anything but the early peak) is far below the
    plate's typical melting signal."""
    tops = [max((p.height for r in res.runs for p in r.peaks if p.role in ("melting", "satellite")),
                default=0.0) for res in results]
    good = [t for t in tops if t > 0]
    if not good:
        return
    med = float(np.median(good))
    for res in results:
        prod = max((p.height for r in res.runs for p in r.peaks if p.role != "anchor"), default=0.0)
        if prod < NO_PRODUCT_RATIO * med and "no product" not in res.flags:
            res.flags.insert(0, "no product")


def plate_wells(run_folders: Sequence) -> List[str]:
    """Wells present in every one of the run folders."""
    names = None
    for f in run_folders:
        stems = {p.stem for p in Path(f).glob("*.rsd")}
        names = stems if names is None else names & stems
    return sorted(names or [])


def group_runs(folders: Sequence) -> Dict[str, List[Path]]:
    """Run folders -> {plate: [run folders]}, runs of one plate together."""
    out: Dict[str, List[Path]] = {}
    for f in folders:
        f = Path(f)
        out.setdefault(plate_of(f.name), []).append(f)
    return out


def analyse_plate(run_folders: Sequence, channel: Optional[int] = None,
                  is_channel: Optional[int] = None, merge_satellites: bool = True,
                  wells: Optional[Sequence[str]] = None,
                  template: Optional[Template] = None) -> List[WellResult]:
    """Every well present in all the given runs of one plate."""
    run_folders = [Path(f) for f in run_folders]
    names = plate_wells(run_folders)
    if wells:
        names = [w for w in names if w in set(wells)]
    results = [analyse_well([f / f"{w}.rsd" for f in run_folders], channel, is_channel,
                            merge_satellites, template) for w in names]
    flag_no_product(results)
    return results
