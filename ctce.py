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
   variant fraction by mass action, p = (A_hom + A_het / 2) / A_total;
7. plate level (quantify_species): each species measured as a block (main
   peak + Taq stutter) on a near-flat baseline, Homo1 / Homo2 / heteroduplexes
   found as a pattern the plate learns in migration-time ratios, clipped
   peaks rebuilt from another channel's crosstalk.

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
    method: str = "change"              # change / shape / template / species
    minor_fraction: Optional[float] = None
    species: List[str] = field(default_factory=list)       # names of the areas (species mode)
    mf_by_temperature: Dict[float, float] = field(default_factory=dict)  # the check runs' MF
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
            "species": "/".join(self.species),
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
# species: main peak + stutter, quantified as a block (plate level)
# --------------------------------------------------------------------------- #
# Each duplex species is a cluster: the main peak, often a shoulder, and its
# Taq stutter about 15 scans later.  The species are found as a pattern the
# plate teaches itself, in migration-time ratios (stable from capillary to
# capillary, unlike scan offsets): Homo1 the largest species, Homo2 after
# Homo1's cluster at a ratio the strong wells share, then up to two
# heteroduplexes.  Peaks that sit at the same ratio to the early peak in most
# wells, without a stutter, are static (ssDNA) and never a species.
# Validated against hand-integrated areas of 010823 #19_20 (94 wells):
# Homo1 and Homo2 areas within 2-3 % (median), minor fraction within 0.03 in
# 85 % of the wells.
SPECIES_BASELINE = 1601       # quantify on a near-flat floor (wide running minimum)
SPECIES_WINDOW = (10, 30)     # a species block: 10 scans before its main peak, 30 after
SPECIES_HUMP = 101            # hump baseline: opening this wide (scans), under the species
SPECIES_STUTTER = 15          # Taq stutter: scans after its main peak (the plate learns its own)
SPECIES_PROMINENCE = 10       # peaks: prominence >= 10 x the noise
SPECIES_MIN_GAP = 40          # Homo2 lies beyond Homo1's own cluster (stutter, shoulders)
SPECIES_MAX_R3 = 1.075        # ... and before 1.075 x Homo1's migration time
SPECIES_R3_SIGMA = 0.015      # spread of Homo2 / Homo1 between capillaries
SPECIES_HET_FROM = 1.5        # hets lie 1.5 to 4 Homo1->Homo2 distances after Homo1
SPECIES_HET_TO = 4.0
SPECIES_HET_EQUAL = 0.7       # the two hets' areas should agree (log-ratio sigma ~x2)
SPECIES_HET_MIN = 25          # heteroduplex main peak >= 25 x the noise
HET_CLUSTER_GAP = 12          # a small bump this close before a bigger heteroduplex peak is neither's main
HET_DOUBLET_SPAN = 18         # a heteroduplex's +A stutter lies at most this many scans after it
MERGED_VALLEY = 0.2           # heteroduplexes whose valley stays above this have begun to merge:
                              # one area under one line (split at the valley only for display)
SPECIES_MAX_WIDTH = 4.0       # a species' main peak is at most 4 x Homo1's width ...
SPECIES_MAX_WIDTH_SCANS = 20  # ... or 20 scans at half height: broad bumps are not species
MERGED_HET_WIDTH = 1.0        # a lone heteroduplex block at least Homo1's width is both, merged
TEMPS_DISAGREE = 1.3          # MF at another temperature off by > x1.3: a peak may hide
RUNS_DISAGREE = 0.5           # runs of one well differing by > 50 % of the mean are flagged
HOMO1_OVER_STATIC = 3.0       # a "static" candidate is still Homo1 if >= 3x as tall as the best other one
SPECIES_BAND = 0.06           # Homo1 / early peak within +-0.06 of the plate's typical value
STATIC_SHARE = 0.6            # static: at one ratio to the early peak in >= 60 % of wells ...
STATIC_CHANGE_SHARE = 0.25    # ... or seen static between temperatures in >= 25 % of them
CLIP_LEVEL = 65000            # raw counts at or above this are clipped (16-bit ADC)
SPECIES_NAMES = ("Homo1", "Homo2", "Hetero1", "Hetero2")


@dataclass
class SpeciesPattern:
    """What a plate taught itself about where its species sit."""
    static: List[float] = field(default_factory=list)   # static peaks / early peak
    band: Optional[tuple] = None                        # Homo1 / early peak range
    r3: Optional[float] = None                          # Homo2 / Homo1
    stutter: float = SPECIES_STUTTER                    # stutter offset (scans); block = 2x
    temperature: Optional[float] = None


def noise_level(y: np.ndarray) -> float:
    """Robust noise of a trace (MAD of the first difference)."""
    d = np.diff(np.asarray(y, float))
    return max(float(1.4826 * np.median(np.abs(d)) / np.sqrt(2)), 1e-9)


def quant_trace(rt: RunTrace) -> np.ndarray:
    """The trace areas are measured on: the sample channel over a wide
    baseline, with clipped stretches rebuilt from the strongest unclipped
    channel (its crosstalk scaled by the gain fitted on the unclipped part of
    the same well) -- measuring an out-of-scale peak in another channel."""
    cached = getattr(rt, "_yq", None)
    if cached is not None:
        return cached
    doc = getattr(rt, "_doc", None)
    if doc is None:
        rt._yq = rt.y
        return rt.y
    order = rt._order
    raw = np.asarray(doc.acgt, float)
    y = channel_trace(raw, order, rt.channel, SPECIES_BASELINE)
    cs = "ACGT".index(order[rt.channel - 1])
    clip = raw[:, cs] >= CLIP_LEVEL
    rt._clipped = int(clip.sum())
    if clip.any():
        best = None
        for k in range(4):
            if k == cs or (raw[:, k] >= CLIP_LEVEL).any() or "ACGT"[k] not in order:
                continue
            yo = channel_trace(raw, order, order.index("ACGT"[k]) + 1, SPECIES_BASELINE)
            good = (~clip) & (y > 0.03 * CLIP_LEVEL) & (yo > 0)
            if good.sum() < 10:
                continue
            corr = float(np.corrcoef(y[good], yo[good])[0, 1])
            if corr <= 0.9:              # not the same dye's crosstalk
                continue
            gain = float((y[good] * yo[good]).sum() / (yo[good] ** 2).sum())
            if best is None or gain < best[0]:      # the strongest copy: least noise
                best = (gain, yo, "ACGT"[k])
        if best is not None:
            from scipy.ndimage import binary_dilation
            zone = binary_dilation(clip, iterations=4)
            y = y.copy()
            y[zone] = np.maximum(y[zone], best[0] * best[1][zone])
            rt._desat = (best[2], best[0])
    rt._yq = y
    return y


def species_trace(rt: RunTrace) -> np.ndarray:
    """The trace species are found and ranked on: the measuring trace with
    the broad hump under the species taken off, so a bump riding on Homo1's
    tail does not look like a block of its own."""
    cached = getattr(rt, "_ys", None)
    if cached is None:
        y = quant_trace(rt)
        cached = rt._ys = y - hump_baseline(y)
    return cached


def hump_baseline(y: np.ndarray, width: int = SPECIES_HUMP) -> np.ndarray:
    """The broad hump under the species (unresolved material, Homo1's long
    tail): a morphological opening as wide as several species blocks,
    smoothed, never above the trace.  Taken off before peaks are found and
    ranked; the blocks themselves are measured above a line joining their
    ends, which also follows a tail sloping under them."""
    from scipy.ndimage import maximum_filter1d, minimum_filter1d, uniform_filter1d
    y = np.asarray(y, float)
    base = uniform_filter1d(maximum_filter1d(minimum_filter1d(y, width), width), width // 2)
    return np.minimum(base, y)


def _species_peaks(y: np.ndarray, anchor: float, nz: float) -> List[int]:
    idx, _ = find_peaks(y, prominence=SPECIES_PROMINENCE * nz, distance=3)
    return [int(i) for i in idx if i > anchor * 1.05]


def _foot(y: np.ndarray, c: int) -> float:
    lo, hi = SPECIES_WINDOW
    return float(np.clip(y[max(0, c - lo):c + hi + 1], 0, None).sum())


def _lead(idx: Sequence[int], y: np.ndarray, c: int, span: int = 25) -> int:
    """The first big peak of the cluster that ends at *c* (a species' main
    peak can be lower than its shoulder or stutter)."""
    return min([j for j in idx if c - span <= j <= c and y[j] > 0.4 * y[c]], default=c)


def find_species(y: np.ndarray, anchor: float, pattern: Optional[SpeciesPattern] = None,
                 static_scans: Sequence[float] = ()) -> Optional[List[int]]:
    """Main-peak scans of [Homo1, Homo2, Hetero1?, Hetero2?] in one trace, or None.
    *static_scans*: peaks seen not to move between temperatures -- never a
    species."""
    pattern = pattern or SpeciesPattern()
    nz = noise_level(y)
    idx = [j for j in _species_peaks(y, anchor, nz) if j > anchor * 1.1]

    def allowed(c, use_static=True):
        if use_static and any(abs(c - s) <= 4 for s in static_scans):
            return False
        r = c / anchor
        if any(abs(r - s) <= 0.012 for s in pattern.static):
            return False
        return pattern.band is None or pattern.band[0] <= r <= pattern.band[1]
    def first_homo1(use_static):
        for c in sorted(idx, key=lambda j: -_foot(y, j)):
            if not allowed(c, use_static):
                continue
            lead = _lead(idx, y, c)
            if any(lead + SPECIES_MIN_GAP <= j <= lead * SPECIES_MAX_R3 and y[j] < y[c] for j in idx):
                return lead
        return None
    h2 = first_homo1(True)
    # Peaks move over one another between temperatures: Homo1 can sit where
    # the other run had a peak and be paired with it as "static".  A
    # candidate far bigger than anything left is still Homo1 (a static peak
    # that big would have been learnt plate-wide).
    loose = first_homo1(False)
    if loose is not None and loose != h2 and (
            h2 is None or y[loose] >= HOMO1_OVER_STATIC * y[h2]):
        h2 = loose
    if h2 is None:
        return None
    c3 = [j for j in idx if h2 + SPECIES_MIN_GAP <= j <= h2 * SPECIES_MAX_R3]
    # a bump on the rise of a much taller peak is not Homo2: its footprint
    # (10 before to 30 after) borrows that peak's area (a 234-high bump beat
    # a 7433-high Homo2 11 scans after it); a stutter is never twice its main
    c3 = [j for j in c3 if not y[j + 1:j + SPECIES_WINDOW[1] + 1].max(initial=0) > 2 * y[j]] or c3
    if pattern.r3 is None:
        score = lambda j: _foot(y, j)                                   # noqa: E731
    else:
        score = lambda j: _foot(y, j) * np.exp(                         # noqa: E731
            -0.5 * ((j / h2 - pattern.r3) / SPECIES_R3_SIGMA) ** 2)
    h3 = max(c3, key=score)
    lead3 = _lead(idx, y, h3)
    if lead3 >= h2 + SPECIES_MIN_GAP:
        h3 = lead3
    # heteroduplexes: the two crossed strand pairings (Watson of one variant
    # with Crick of the other, and the reverse) -- equal amounts in theory,
    # so they are chosen as the pair with the most similar, largest blocks
    st = stutter_of(y, h2, pattern.stutter)
    # the heteroduplexes' stutter sits closer than Homo1's (040723 A01: Homo1
    # +24, Homo2 +17, heteroduplexes +16 / +13): take Homo2's, the nearest kin
    st = stutter_of(y, h3, st) if y[h3] > 0 else st
    d = h3 - h2
    lo4 = max(h3 + SPECIES_MIN_GAP, h2 + SPECIES_HET_FROM * d)
    het_end = h2 + SPECIES_HET_TO * d
    c4 = sorted(j for j in idx if lo4 <= j <= het_end and y[j] > SPECIES_HET_MIN * nz)
    # A heteroduplex shows as a doublet: its main peak and its Taq +A
    # stutter 8-16 scans on.  The spacing inside a doublet overlaps the
    # spacing between the two heteroduplexes (021123: 8-16 vs 11-16), so a
    # gap rule cannot tell them apart; pairing the peaks in order can -- the
    # first with the next one within HET_DOUBLET_SPAN, then the next free
    # one -- which matches every operator mark (Run01 A05: 2720 + 2733 and
    # 2770 + 2782, not 2720 and its own stutter).  A small bump just before
    # a bigger peak is neither (040723 A01 Hetero2 at 2865, not 2854).
    c4 = [j for j in c4 if not any(0 < k - j <= HET_CLUSTER_GAP and y[j] < 0.4 * y[k] for k in c4)]
    clusters: List[List[int]] = []
    for j in c4:
        if clusters and len(clusters[-1]) < 2 and j - clusters[-1][0] <= HET_DOUBLET_SPAN:
            clusters[-1].append(j)
        else:
            clusters.append([j])
    # A cluster is a stutter only when it starts about one stutter distance
    # after another cluster *and* is clearly smaller -- the two
    # heteroduplexes themselves can sit a stutter apart and are about as
    # large as each other.
    top = {id(c): max(y[j] for j in c) for c in clusters}
    keep = []
    for c in clusters:
        stutter_of_prev = any(st - 4 <= c[0] - p[-1] <= st + 6 and top[id(c)] < 0.6 * top[id(p)]
                              for p in clusters if p is not c and p[-1] < c[0])
        if not stutter_of_prev:
            keep.append(c)
    # a cluster's main peak: its first big peak (a small bump just before
    # the main one is not it -- 040723 A01 Hetero2 at 2865, not 2854)
    lead = {id(c): next(j for j in c if y[j] > 0.4 * top[id(c)]) for c in keep}
    area = {id(c): float(np.clip(y[max(0, c[0] - 6):c[-1] + 7], 0, None).sum()) for c in keep}
    best, best_s = [], 0.0
    for i, a in enumerate(keep):
        for b in keep[i + 1:]:
            if min(area[id(a)], area[id(b)]) <= 0:
                continue
            q = np.log(area[id(b)] / area[id(a)])
            sc = (area[id(a)] + area[id(b)]) * np.exp(-0.5 * (q / SPECIES_HET_EQUAL) ** 2)
            if sc > best_s:
                best, best_s = [lead[id(a)], lead[id(b)]], sc
    if not best and keep:                                   # one het (the other merged / lost)
        best = [lead[id(max(keep, key=lambda c: area[id(c)]))]]
    hets = best
    return [h2, h3] + sorted(hets)


def _front(ys: np.ndarray, c: int) -> int:
    """Foot of a main peak's sharp front (hump-free trace), at most 10 scans
    back: walking down the (smoothed) front until it is flat (below 1 % of
    the peak) or turns up again (the valley after a peak before it)."""
    sm = np.convolve(ys, np.ones(5) / 5, mode="same")
    top = ys[max(0, c - 2):c + 3].max()
    k = c
    while k > max(1, c - SPECIES_WINDOW[0] - 4) and sm[k] >= 0.01 * top and sm[k - 1] <= sm[k]:
        k -= 1
    return k


def _back(ys: np.ndarray, c: int, stutter: float) -> int:
    """End of a species block: past its stutter -- the first peak 8 to
    stutter+6 scans on (a bump beyond it is not the species') -- and down
    the stutter's tail until flat (below 1 % of the main peak) or turning up
    again; at most twice the stutter offset."""
    top = ys[max(0, c - 2):c + 3].max()
    cap = min(len(ys) - 2, c + int(round(2 * stutter)))
    sm = np.convolve(ys, np.ones(3) / 3, mode="same")
    lo, hi = c + 8, min(len(sm) - 2, c + int(round(stutter)) + 6)
    # the stutter is the tallest peak near the stutter offset, not the first
    # bump after the main peak: a shoulder 7-8 scans on ended the block in
    # the valley before an equally tall Taq +A cluster (Homo1 at half size)
    near = max(lo, c + int(round(stutter)) - 6)
    tops = [i for i in range(near, hi + 1)
            if sm[i] >= sm[i - 1] and sm[i] >= sm[i + 1] and sm[i] > 0.1 * top]
    if tops:
        k = max(tops, key=lambda i: sm[i])
    else:
        tops = [i for i in range(lo, hi + 1)
                if sm[i] >= sm[i - 1] and sm[i] >= sm[i + 1] and sm[i] > 0.1 * top]
        k = tops[0] if tops else (lo + int(np.argmax(sm[lo:hi + 1])) if hi >= lo else c)
    # down the stutter's tail -- through its own sub-peaks (a doublet), up
    # to a real valley (below 20 % of the stutter peak) or flat ground
    ktop = sm[k]
    while k < cap and sm[k] >= 0.01 * top:
        if sm[k + 1] > sm[k] and sm[k] < 0.2 * ktop:
            break
        k += 1
    return k


def _single_area(ys: np.ndarray, c: int, stutter: float) -> float:
    """Block area of one species on its own (hump-free trace, line baseline)."""
    lo, hi = _front(ys, c), _back(ys, c, stutter)
    return block_area(ys, lo, hi)


def species_blocks(y: np.ndarray, mains: Sequence[int],
                   stutter: float = SPECIES_STUTTER) -> List[tuple]:
    """(start, stop, base at start, base at stop) of each species block.

    A block runs from the foot of its main peak's front to the end of its
    stutter's tail (at most twice the stutter offset), and is measured above
    a straight line -- as a drag between two points would be.  Neighbours
    that would overlap split at the valley after the first one's stutter:
    if that valley is back near the floor each keeps its own line (a bump
    between them is neither's); if they have begun to merge they share one line
    from the first's start to the last's end, split by a drop at the valley."""
    y = np.asarray(y, float)
    ys = y - hump_baseline(y)
    hi_w = int(round(2 * stutter))
    bounds = [[_front(ys, c), _back(ys, c, stutter)] for c in mains]
    joined = [False] * len(mains)
    lo_w = SPECIES_WINDOW[0]
    for k in range(len(mains) - 1):
        c, n = mains[k], mains[k + 1]
        if n - c >= lo_w + hi_w or n - c < 6:
            continue
        a = c + min(int(round(stutter)) + 5, int(0.6 * (n - c)))
        v = int(a + np.argmin(y[a:n - 1])) if n - 1 > a else (c + n) // 2
        low = min(ys[max(0, c - 2):c + 3].max(), ys[max(0, n - 2):n + 3].max())
        bounds[k][1] = v
        # Neighbours that have begun to merge (the valley stays above
        # MERGED_VALLEY of the smaller peak) are measured as one area under
        # one line, split by a drop at the valley only to show each -- the
        # operator's rule: Hetero1 and Hetero2 each count half towards the
        # homoduplexes, so their split does not matter, and two lines cut
        # at a raised valley lose what lies below it (more even plates:
        # robust CV 021123 0.23 -> 0.21, 010823 0.25 -> 0.23, 040723
        # 0.22 -> 0.20).  Resolved ones (valley near the floor) keep their
        # own lines, so ground between them is neither's.
        if ys[v] > MERGED_VALLEY * low:
            joined[k] = True
            bounds[k + 1][0] = v + 1
        else:
            bounds[k + 1][0] = max(v + 1, bounds[k + 1][0]) if ys[v] <= 0.05 * low else v + 1
    out: List[tuple] = []
    k = 0
    while k < len(mains):
        g = k
        while g < len(mains) - 1 and joined[g]:
            g += 1
        lo, hi = bounds[k][0], bounds[g][1]
        a, b = _block_ends(y, lo, hi)
        line = lambda x: a + (b - a) * (x - lo) / max(hi - lo, 1)   # noqa: E731
        for j in range(k, g + 1):
            out.append((bounds[j][0], bounds[j][1], line(bounds[j][0]), line(bounds[j][1])))
        k = g + 1
    return out


def species_bounds(y: np.ndarray, mains: Sequence[int],
                   stutter: float = SPECIES_STUTTER) -> List[List[int]]:
    """[start, stop] of each species block (see species_blocks)."""
    return [[lo, hi] for lo, hi, _a, _b in species_blocks(y, mains, stutter)]


def stutter_of(ys: np.ndarray, h2: int, default: float = SPECIES_STUTTER) -> float:
    """This well's stutter offset: the tallest peak 8-35 scans after Homo1
    (hump-free trace), else the plate's value."""
    seg = ys[h2 + 8:h2 + 36]
    if len(seg) and seg.max() > 0.1 * ys[h2]:
        k = int(np.argmax(seg))
        if 0 < k < len(seg) - 1:
            return float(np.clip(k + 8, 10, 30))
    return default


def _block_ends(y: np.ndarray, lo: int, hi: int) -> tuple:
    """Baseline levels at a block's two ends (each the mean of 3 scans)."""
    return (float(np.mean(y[max(0, lo - 1):lo + 2])), float(np.mean(y[max(0, hi - 1):hi + 2])))


def block_area(y: np.ndarray, lo: int, hi: int, a: Optional[float] = None,
               b: Optional[float] = None) -> float:
    """Area of y[lo..hi] above the straight line from level *a* at its start
    to *b* at its end (default: the trace's own ends, each the mean of 3
    scans) -- what a drag between two points measures, so a broad hump
    under the species is left out."""
    seg = np.asarray(y[lo:hi + 1], float)
    if len(seg) < 3:
        return 0.0
    if a is None or b is None:
        a, b = _block_ends(y, lo, hi)
    return float(np.clip(seg - np.linspace(a, b, len(seg)), 0, None).sum())


def species_areas(y: np.ndarray, mains: Sequence[int],
                  stutter: float = SPECIES_STUTTER) -> List[float]:
    """Block area of each species (main peak + shoulder + stutter) above its
    baseline line (see species_blocks)."""
    return [block_area(y, lo, hi, a, b) for lo, hi, a, b in species_blocks(y, mains, stutter)]


def species_minor_fraction(areas: Sequence[float]) -> Optional[float]:
    """Mass action over [Homo1, Homo2, Het...]: p = (A_Homo2 + A_het / 2) / A_total."""
    a = list(areas) + [0.0] * (4 - len(areas))
    total = sum(a)
    if total <= 0:
        return None
    p = (a[1] + (a[2] + a[3]) / 2.0) / total
    return min(p, 1.0 - p)


def _fwhm(y: np.ndarray, c: int) -> int:
    """Width at half height, each side stopping at a valley: a peak merged
    with its neighbour (two heteroduplexes, a doublet) never falls to half
    height between them and read as a broad bump (021123 Run02 A01 lost its
    Hetero2), while a true broad bump has no valley inside."""
    sm = np.convolve(y, np.ones(3) / 3, mode="same")
    lvl = y[c] * 0.5
    w = 0
    for step in (-1, 1):
        j = c
        while 0 < j < len(y) - 1 and sm[j] > lvl and not (
                abs(j - c) > 2 and sm[j + step] > sm[j] and sm[j] < 0.9 * y[c]):
            j += step
        w += abs(j - c)
    return w


def present_species(ys: np.ndarray, mains: Sequence[int]):
    """(names, scans) of the species that are really there.

    Homo1 always; Homo2 and each heteroduplex only if its main peak stands
    SPECIES_HET_MIN x the noise above the trace and is a peak, not a broad
    bump (at most SPECIES_MAX_WIDTH x Homo1's width).  At a low minor
    fraction Homo2 (~p^2) drops out first, leaving Homo1 and the two
    heteroduplexes."""
    nz = noise_level(ys)
    w1 = max(_fwhm(ys, mains[0]), 2)
    names, kept = [SPECIES_NAMES[0]], [mains[0]]
    for name, c in zip(SPECIES_NAMES[1:], mains[1:]):
        if ys[c] >= SPECIES_HET_MIN * nz and \
                _fwhm(ys, c) <= max(SPECIES_MAX_WIDTH * w1, SPECIES_MAX_WIDTH_SCANS):
            names.append(name)
            kept.append(c)
    return names, kept


def merged_heteroduplex(names: Sequence[str], blocks: Sequence[tuple]) -> List[str]:
    """Rename a lone heteroduplex block "Hetero1+2" when it is the two
    heteroduplexes merged: wider than Homo1's own block (main + stutter).
    A merged block holds both and counts once; a lone, narrow one lost its
    partner, whose area it then stands in for (minor_fraction_named)."""
    names = list(names)
    hets = [i for i, n in enumerate(names) if n.startswith("Hetero")]
    if len(hets) == 1 and blocks:
        w1 = blocks[0][1] - blocks[0][0]
        i = hets[0]
        if blocks[i][1] - blocks[i][0] >= MERGED_HET_WIDTH * w1:
            names[i] = "Hetero1+2"
    return names


def minor_fraction_named(areas: Dict[str, float]) -> Optional[float]:
    """The minor fraction from the species present:
    p = (A_Homo2 + (A_Hetero1 + A_Hetero2) / 2) / A_total.

    It holds whatever the share of heteroduplexes: a well with little
    template is still amplifying when the program ends and forms few
    heteroduplexes, so p rests on the homoduplexes, Homo2 / (Homo1 + Homo2);
    at a low fraction Homo2 (~p^2) drops out and p rests on the
    heteroduplexes.  Neither Homo2 nor a heteroduplex seen: None -- the
    minor species is below detection and no number is forced out of noise.
    One heteroduplex lost: its partner's area stands in for it (the two
    crossed pairings are equal in theory)."""
    hets = [areas[n] for n in ("Hetero1", "Hetero2") if n in areas]
    h1, h2 = areas.get("Homo1", 0.0), areas.get("Homo2", 0.0)
    if "Hetero1+2" in areas:                    # merged: both in one block
        het_sum = areas["Hetero1+2"]
        total = h1 + h2 + het_sum
        if total <= 0:
            return None
        p = (h2 + het_sum / 2.0) / total
        return min(p, 1.0 - p)
    if not hets and "Homo2" not in areas:
        return None
    het_each = float(np.mean(hets)) if hets else 0.0
    total = h1 + h2 + 2 * het_each
    if total <= 0:
        return None
    p = (h2 + het_each) / total
    return min(p, 1.0 - p)


def learn_species_pattern(runs: Sequence[RunTrace]) -> SpeciesPattern:
    """Teach the plate its pattern from one run per well (same temperature):
    the static peaks, Homo1's band, and Homo2's ratio to Homo1."""
    traces = [(r, species_trace(r)) for r in runs if r.anchor]
    pat = SpeciesPattern(temperature=next((r.temperature for r in runs), None))
    if not traces:
        return pat
    # static: same ratio to the early peak in most wells, no stutter, not the largest species
    bins = np.arange(1.05, 1.6, 0.004)
    hits = np.zeros(len(bins))
    for r, y in traces:
        nz = noise_level(y)
        ip = _species_peaks(y, r.anchor, nz)
        if not ip:
            continue
        top = max(_foot(y, i) for i in ip)
        ratios = np.array([i / r.anchor for i in ip
                           if y[i] > 300 * nz and _foot(y, i) < 0.5 * top
                           and not any(8 <= j - i <= 24 and 0.15 * y[i] < y[j] < y[i] for j in ip)
                           and not any(8 <= i - j <= 24 and y[j] > y[i] for j in ip)])
        if len(ratios):
            hits += np.array([(np.abs(ratios - c) <= 0.008).any() for c in bins])
    for k in np.argsort(-hits):
        if hits[k] < STATIC_SHARE * len(traces):
            break
        if all(abs(bins[k] - s) > 0.02 for s in pat.static):
            pat.static.append(round(float(bins[k]), 4))
    # ... and wherever the temperature change showed static peaks in a good
    # share of the wells (they are only paired where their heights agree)
    moved = [r for r, _y in traces if getattr(r, "_static_scans", None) is not None]
    if moved:
        hs = np.zeros(len(bins))
        for r in moved:
            ratios = np.array([s / r.anchor for s in _static_scans(r) if s > r.anchor * 1.05])
            if len(ratios):
                hs += np.array([(np.abs(ratios - c) <= 0.006).any() for c in bins])
        for k in np.argsort(-hs):
            if hs[k] < STATIC_CHANGE_SHARE * len(moved):
                break
            if all(abs(bins[k] - s) > 0.02 for s in pat.static):
                pat.static.append(round(float(bins[k]), 4))
    # Homo1 band and Homo2 ratio from the strongest third of the wells
    first = []
    for r, y in traces:
        sp = find_species(y, r.anchor, SpeciesPattern(static=pat.static), _static_scans(r))
        if sp:
            first.append((_foot(y, sp[0]), sp[0], r, y))
    if not first:
        return pat
    first.sort(key=lambda t: -t[0])
    strong = first[:max(5, len(first) // 3)]
    mid = float(np.median([h2 / r.anchor for _f, h2, r, _y in strong]))
    pat.band = (mid - SPECIES_BAND, mid + SPECIES_BAND)
    # Homo2: the dominant block after Homo1's own cluster, in each strong well
    r3s, stut = [], []
    for _f, h2, r, y in strong:
        ip = _species_peaks(y, r.anchor, noise_level(y))
        cand = [j for j in ip if h2 + SPECIES_MIN_GAP <= j <= h2 * SPECIES_MAX_R3
                and y[j] > 0.03 * y[h2]
                and not any(abs(j / r.anchor - s) <= 0.012 for s in pat.static)]
        if cand:
            j = max(cand, key=lambda j: _foot(y, j))
            r3s.append(_lead(ip, y, j) / h2)
        # the stutter: the tallest peak 8-35 scans after Homo1
        after = [j for j in ip if h2 + 8 <= j <= h2 + 35 and y[j] > 0.1 * y[h2]]
        if after:
            stut.append(max(after, key=lambda j: y[j]) - h2)
    if r3s:
        pat.r3 = round(float(np.median(r3s)), 4)
    if len(stut) >= 3:
        pat.stutter = float(np.clip(np.median(stut), 10, 30))
    return pat


def _static_scans(run: RunTrace) -> List[float]:
    """Scans of the peaks the temperature change showed to be static."""
    return list(getattr(run, "_static_scans", ()) or ())


def _label_species(run: RunTrace, mains: Sequence[int], stutter: float = SPECIES_STUTTER,
                   names: Optional[Sequence[str]] = None) -> None:
    """Mark the species on the run's peaks for display: the main peak red and
    labelled, the rest of its block as its satellites; the early peak is the
    primer peak."""
    for p in run.peaks:
        if p.role in ("melting", "satellite"):
            p.role, p.parent = "other", None
        if p.role == "anchor":
            p.label = "primer peak"
    lo_w, hi_w = SPECIES_WINDOW[0], int(round(2 * stutter))
    for name, c in zip(names or SPECIES_NAMES, mains):
        near = [k for k, p in enumerate(run.peaks) if abs(p.raw_pos - c) <= 3]
        if near:
            k = min(near, key=lambda k: abs(run.peaks[k].raw_pos - c))
        else:                                   # found on the quantification trace only
            run.peaks.append(Peak(pos=int(round(run.scale * c + run.offset)), height=float(run.y[c]),
                                  left=c - lo_w, right=c + hi_w, area=0.0, front=0, tail=0,
                                  raw_pos=c))
            k = len(run.peaks) - 1
        run.peaks[k].role, run.peaks[k].label = "melting", name
        for j, p in enumerate(run.peaks):
            if j != k and c + 3 < p.raw_pos <= c + hi_w and p.role not in ("anchor", "melting"):
                p.role, p.parent = "satellite", k


def _measure_run(run: RunTrace, pat: SpeciesPattern):
    """(names, areas) of one run's species, measured with *pat*, or None."""
    ys = species_trace(run)
    mains = find_species(ys, run.anchor, pat, _static_scans(run))
    if not mains:
        return None
    names, kept = present_species(ys, mains)
    st = stutter_of(ys, kept[0], pat.stutter)
    yq = quant_trace(run)
    blocks = species_blocks(yq, kept, st)
    return merged_heteroduplex(names, blocks), [block_area(yq, *bl) for bl in blocks]


def _measure_other_temperatures(res: WellResult, t_meas, runs_at, temps, pats, names, areas) -> None:
    """Measure the well at every other temperature too (first run there),
    kept only where that temperature resolves the same species clearly
    (each >= 1 % of the total) -- too cold they have not separated, too hot
    they merged again.  Stored for the plate-level check."""
    res._temp_check = {}
    if res.minor_fraction is None or not temps or len(temps) < 2 or sum(areas) <= 0:
        return
    for t in temps:
        if t_meas is not None and abs(t - t_meas) <= SAME_TEMPERATURE:
            continue
        for run in [x for x in runs_at(res, t) if x.anchor][:1]:
            got = _measure_run(run, pats[t])
            if not got:
                continue
            n2, a2 = got
            mf2 = minor_fraction_named(dict(zip(n2, a2)))
            if (mf2 is None or sorted(n2) != sorted(names)
                    or min(a2) < 0.01 * sum(a2) or min(areas) < 0.01 * sum(areas)):
                continue
            res.mf_by_temperature[t] = mf2
            res._temp_check[t] = dict(zip(n2, (a / sum(a2) for a in a2)))


def _flag_temperature_disagreement(todo: Sequence[WellResult], t_meas) -> None:
    """Peaks move over one another as the temperature changes: a peak
    resolved at one temperature can hide under a species block at another
    and inflate it -- a big one changes the fraction.  The same DNA must give
    the same fraction at every temperature *relative to the plate* (a
    temperature can shift the whole plate, e.g. a broad peak sitting under
    Homo2 on every well at 53 deg C), so a well is flagged when its fraction
    relative to the plate median differs between temperatures by more than
    TEMPS_DISAGREE; the species whose share moved most is named."""
    meas = [r.minor_fraction for r in todo if r.minor_fraction is not None]
    if not meas:
        return
    med = {None: float(np.median(meas))}
    temps = sorted({t for r in todo for t in r.mf_by_temperature})
    for t in temps:
        v = [r.mf_by_temperature[t] for r in todo if t in r.mf_by_temperature]
        if len(v) >= 5:
            med[t] = float(np.median(v))
    for res in todo:
        if res.minor_fraction is None or not getattr(res, "_temp_check", None):
            continue
        rel = res.minor_fraction / med[None]
        share = dict(zip(res.species, (a / sum(res.areas) for a in res.areas)))
        notes = []
        for t, mf2 in res.mf_by_temperature.items():
            if t not in med or med[t] <= 0:
                continue
            rel2 = mf2 / med[t]
            if max(rel, rel2) / max(min(rel, rel2), 1e-9) <= TEMPS_DISAGREE:
                continue
            s2 = res._temp_check[t]
            # the species that moved most, after the plate-wide shift
            common = [n for n in share if n in s2 and min(share[n], s2[n]) > 0]
            moved = max(common, key=lambda n: abs(np.log(share[n] / s2[n])), default=None)
            where = f", {moved} {share[moved]:.3f} vs {s2[moved]:.3f} of the total" if moved else ""
            notes.append(f"{t:g} °C gives {mf2:.3f} (plate {med[t]:.3f}){where}")
        if notes:
            tm = f"{t_meas:g} °C" if t_meas is not None else "the measured run"
            res.flags.append(f"temperatures disagree (MF {res.minor_fraction:.3f} at {tm}, plate "
                             f"{med[None]:.3f}; " + "; ".join(notes)
                             + "): a peak may hide under a species at one of them")


def quantify_species(results: Sequence[WellResult]) -> Optional[SpeciesPattern]:
    """Plate level: find Homo1 / Homo2 / heteroduplexes in every well as a
    pattern the plate teaches itself, and measure them as blocks (main peak
    + stutter).  Replaces the areas, fractions and minor fraction of each
    well not guided by a template.  Measured at the temperature where the
    most wells resolve all four species (ties: the warmer), replicate runs
    there averaged; peaks the temperature change showed to be static are
    never a species."""
    todo = [r for r in results if r.method != "template" and r.runs]
    if not todo:
        return None
    for res in todo:                     # what the temperatures showed to be static
        if res.method == "change":
            for run in res.runs:
                run._static_scans = [p.raw_pos for p in run.peaks if p.role == "static"]
    temps = sorted({round(x.temperature, 1) for r in todo for x in r.runs if x.temperature is not None})

    def runs_at(res, t):
        if t is None:
            return list(res.runs)
        return [x for x in res.runs if x.temperature is not None
                and abs(x.temperature - t) <= SAME_TEMPERATURE]
    # the measuring temperature: where the most wells resolve all four species
    best_t, best_n, pattern = None, -1, SpeciesPattern()
    pats: Dict[Optional[float], SpeciesPattern] = {}
    for t in (temps or [None]):
        firsts = [rr[0] for rr in (runs_at(r, t) for r in todo) if rr]
        pat = learn_species_pattern(firsts)
        pat.temperature = t
        pats[t] = pat
        n = sum(1 for x in firsts if x.anchor and len(
            find_species(species_trace(x), x.anchor, pat, _static_scans(x)) or []) == 4)
        if n >= best_n:                  # ties: the warmer run
            best_t, best_n, pattern = t, n, pat
    t_meas = best_t
    pattern.temperature = t_meas

    def temp_of(run):
        if run.temperature is None or not pats:
            return next(iter(pats), None)
        return min(pats, key=lambda t: abs((t if t is not None else 0) - run.temperature))

    for res in todo:
        # Measured at the temperature that resolves the species best
        # (replicate runs there averaged): at the other temperature the
        # heteroduplexes merge or shift onto other peaks, and averaging it in
        # widened the spread over 96 replicate wells (robust CV 0.19 -> 0.34).
        found = []
        for run in [x for x in runs_at(res, t_meas) if x.anchor]:
            pat = pats.get(temp_of(run), pattern)
            ys = species_trace(run)
            mains = find_species(ys, run.anchor, pat, _static_scans(run))
            if not mains:
                continue
            names, kept = present_species(ys, mains)
            st = stutter_of(ys, kept[0], pat.stutter)
            yq = quant_trace(run)
            blocks = species_blocks(yq, kept, st)
            areas = [block_area(yq, *bl) for bl in blocks]
            names = merged_heteroduplex(names, blocks)
            run._species_overlay = {
                "channel": run.channel, "yq": yq, "hump": yq - ys,
                "blocks": [(n,) + tuple(bl) for n, bl in zip(names, blocks)]}
            found.append((run, kept, names, areas,
                          minor_fraction_named(dict(zip(names, areas)))))
        res.flags = [f for f in res.flags if not f.startswith(("no melting", "satellite"))
                     and "satellite peak(s)" not in f]
        if not found:
            res.flags.append("species: Homo1/Homo2 not found")
            continue
        # the run shown: the one with the strongest Homo1
        run, mains, names, areas, _p = max(found, key=lambda f: f[3][0])
        _label_species(run, mains, pattern.stutter, names)
        for x, m, n, _a, p in found:
            k = res.runs.index(x)
            if k < len(res.regimes):
                res.regimes[k] = ", ".join(n) if p is not None else "below detection"
        res.flags = [f.replace("one temperature: peaks classified by shape", "one temperature")
                     for f in res.flags]
        res.best = res.runs.index(run)
        res.method = "species"
        res.species = list(names)
        res.areas = areas
        total = sum(areas)
        res.fractions = [a / total for a in areas] if total > 0 else []
        res.fraction_sd = []
        # A number only where every run sees the minor species: a fraction
        # forced out of noise in one run is not reported (below detection).
        # The runs are at different temperatures and do not resolve the
        # heteroduplexes equally well: a well is measured from the runs that
        # see its minor species clearly (see present_species).  Where none
        # does, it is below detection -- no number forced out of noise.
        n_runs = len([x for x in runs_at(res, t_meas) if x.anchor])
        mfs = [f[4] for f in found if f[4] is not None]
        if not mfs:
            res.minor_fraction = None
            res.flags.append("below detection: neither Homo2 nor a heteroduplex is seen clearly")
        else:
            if len(mfs) < n_runs:
                res.flags.append(f"minor species seen in {len(mfs)} of {n_runs} runs")
            res.minor_fraction = float(np.mean(mfs))
            if len(mfs) > 1:
                res.fraction_sd = [float(np.std(mfs, ddof=1))]
                spread = (max(mfs) - min(mfs)) / max(res.minor_fraction, 1e-9)
                if spread > RUNS_DISAGREE:
                    res.flags.append(f"runs disagree: minor fraction {', '.join(f'{p:.3f}' for p in mfs)}")
        _measure_other_temperatures(res, t_meas, runs_at, temps, pats, names, areas)
        if "Homo2" not in names and res.minor_fraction is not None:
            res.flags.append("Homo2 not seen (low minor fraction): fraction from the heteroduplexes")
        clipped = [x for x, *_ in found if getattr(x, "_desat", None)]
        if clipped:
            ch, gain = clipped[0]._desat
            res.flags.append(f"clipped peak measured in channel {ch} (x{gain:.1f})")
        elif any(getattr(x, "_clipped", 0) for x, *_ in found):
            res.flags.append("clipped peak (no clean channel to measure it in)")
    _flag_temperature_disagreement(todo, t_meas)
    return pattern


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
                  template: Optional[Template] = None,
                  species: bool = True) -> List[WellResult]:
    """Every well present in all the given runs of one plate.  With
    *species* (the default), wells not guided by a template are quantified
    as Homo1 / Homo2 / heteroduplex blocks learnt over the plate."""
    run_folders = [Path(f) for f in run_folders]
    names = plate_wells(run_folders)
    if wells:
        names = [w for w in names if w in set(wells)]
    results = [analyse_well([f / f"{w}.rsd" for f in run_folders], channel, is_channel,
                            merge_satellites, template) for w in names]
    if species:
        quantify_species(results)
    flag_no_product(results)
    return results
