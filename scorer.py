#!/usr/bin/env python3
"""Limoncello scorer -- genotype calls, run QC and a small HTML report.

This is the piece the repository name ("Limoncello-scorer") promises but did
not ship.  It is completely headless (no tkinter, no matplotlib) and works on
two kinds of input:

* **Picked-peak rows** -- the table exported by ``genotyping.py`` (CSV / JSON /
  XLSX), which is also the ML training library.  :func:`score_rows` adds a
  genotype ``call``, the minor-allele ``fraction`` and a ``confidence`` to every
  peak; :func:`well_summary` turns a sample's rows into one genotype record.
* **Trace files** -- ``.rsd`` / ``.scf`` / ``.ab1`` / text, base-called with the
  bundled caller.  :func:`well_qc` returns per-well read statistics (length,
  Qmean/Qmin, N, peak spacing, per-channel signal) and a pass/fail flag.

Confidence is a **documented heuristic**, not a calibrated probability: it
combines how far the allele balance sits from the decision boundary, how well
the two allele peaks stand out of the noise, and whether a CTC-CE internal
standard (HOM/HET duplexes) was picked.  See ``CONFIDENCE_NOTE``.

CLI::

    python scorer.py qc    example_data/M13/*.rsd --report run.html
    python scorer.py peaks picks.csv --out scored.csv
    python scorer.py ladders
    python scorer.py size  well.rsd --ladder geneflo1000_rox --out sizes.csv
    python scorer.py check
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import statistics
import sys
from pathlib import Path
from typing import Iterable, List, Optional

CHANNEL_ORDER = "ACGT"

# A two-allele position is called heterozygous when the minor peak is between
# these fractions of the summed area; outside the band it is a homozygote.
GENOTYPE_HET_LO = 0.15
GENOTYPE_HET_HI = 0.85
# Peaks closer than this many scans count as one allelic position (a single
# repeat, ~one base) -- matches genotyping.HET_WINDOW.
ALLELE_WINDOW_MIN = 8.0
ALLELE_WINDOW_SPACING = 1.6
# Peak area/height below this fraction of the cluster's largest is treated as
# baseline noise rather than an allele.
MIN_ALLELE_FRAC = 0.03

# Read QC gates.  Quality values from the bundled caller are signal heights on
# a per-preset scale, so only scale-free measures gate pass/fail.
MIN_READ_LEN = 50
MAX_N_FRACTION = 0.05
MAX_LOW_SIGNAL_FRAC = 0.35

CONFIDENCE_NOTE = (
    "Confidence 0-100 is a heuristic: 70% weighting on the distance of the "
    "minor-allele fraction from the call boundary (so a clean 50/50 het or a "
    "clean 100/0 hom scores high and a ~20/80 borderline scores low), plus 30% "
    "for allele separation (both peaks standing clear of the cluster's "
    "smallest).  It is not a calibrated probability."
)

CALLS = ("hom-major", "hom-1", "hom-2", "het", "uncertain", "no-call")

# genotyping.export_rows() names the peak size columns height_V / area_Vscan;
# hand-made training tables often just say height / area.  Accept both, plus
# the column the picker writes in memory.
_AREA_KEYS = ("area", "area_Vscan", "area_V")
_HEIGHT_KEYS = ("height", "height_V", "height_Vscan")


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def _num(value, default=0.0) -> float:
    try:
        if value is None or value == "":
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _pick(row, keys, default=0.0) -> float:
    """First present, parseable value among *keys* (alias-tolerant)."""
    for k in keys:
        if k in row and row[k] not in (None, ""):
            return _num(row[k], default)
    return default


def peak_area(row) -> float:
    return _pick(row, _AREA_KEYS)


def peak_height(row) -> float:
    return _pick(row, _HEIGHT_KEYS)


def _norm_key(row) -> tuple:
    return (str(row.get("file", "")), str(row.get("well", "")),
            int(_num(row.get("channel"), 0)))


def is_main(row) -> bool:
    return str(row.get("kind", "main")).lower() in ("main", "")


def has_ambiguous_satellite(cluster: Iterable) -> bool:
    """True when a ``+A`` in this cluster may really be a minor allele.

    ``genotyping._shoulders`` must tag a trailing same-channel peak as either
    ``+A`` or a second main, and a heterozygote's minor allele one base
    downstream is the case it gets wrong: tagged ``+A`` it is excluded from
    ``group_alleles``, so the position scores as a confident homozygote with
    fraction 0.0.  When the peak is a substantial fraction of its parent the
    picker records that, and the only honest answer is "unresolved".
    """
    for rec in cluster:
        if not is_main(rec) and rec.get("ambiguous"):
            return True
    return False


def median_peak_spacing(rows: Iterable) -> float:
    """Median scan gap between consecutive main peaks on one channel."""
    scans = sorted(_num(r.get("scan")) for r in rows if is_main(r))
    gaps = [b - a for a, b in zip(scans, scans[1:]) if b > a]
    return float(statistics.median(gaps)) if gaps else float(ALLELE_WINDOW_MIN)


def allele_window(spacing: Optional[float]) -> float:
    sp = spacing if spacing and spacing > 0 else ALLELE_WINDOW_MIN
    return max(ALLELE_WINDOW_MIN, ALLELE_WINDOW_SPACING * sp)


# --------------------------------------------------------------------------- #
# genotype calling
# --------------------------------------------------------------------------- #
def group_alleles(rows: List[dict], spacing: Optional[float] = None) -> List[List[dict]]:
    """Group main peaks into allelic positions.

    Two main peaks belong together when they sit on the same channel within one
    repeat of each other (the spacing-scaled heterozygote window).  Rows for
    different samples/wells/channels never mix.
    """
    window = allele_window(spacing)
    out: List[List[dict]] = []
    for key in sorted({_norm_key(r) for r in rows if is_main(r)}):
        ch = sorted((r for r in rows if is_main(r) and _norm_key(r) == key),
                    key=lambda r: _num(r.get("scan")))
        for rec in ch:
            placed = False
            for cl in out:
                if (cl and _norm_key(cl[-1]) == key
                        and abs(_num(cl[-1].get("scan"))
                                - _num(rec.get("scan"))) <= window):
                    cl.append(rec)
                    placed = True
                    break
            if not placed:
                out.append([rec])
    return out


def _duplex_labels(cluster: Iterable) -> List[str]:
    return [str(r.get("duplex") or "").upper() for r in cluster]


def genotype_call(cluster: List[dict]) -> str:
    """Return one of ``CALLS`` for an allelic cluster.

    Uses the CTC-CE duplex labels when present (any HET peak proves a
    heterozygote), otherwise the minor-allele fraction of the two largest areas.
    """
    labels = _duplex_labels(cluster)
    has_hom = any(lab.startswith("HOM") for lab in labels)
    if any(lab.startswith("HET") for lab in labels):
        return "het"
    if len(cluster) < 2:
        return "hom-major" if has_hom else "no-call"
    areas = sorted((peak_area(r) for r in cluster), reverse=True)
    top = areas[:2]
    if sum(top) <= 0:
        return "no-call"
    if top[0] <= 0 or top[1] < MIN_ALLELE_FRAC * top[0]:
        # second peak is just noise: one allele only, so a HOM standard label
        # is the only thing that lets us call it without a reference.
        return "hom-major" if has_hom else "no-call"
    frac = top[1] / sum(top)
    if GENOTYPE_HET_LO <= frac <= GENOTYPE_HET_HI:
        return "het"
    return "hom-major"


def minor_fraction(cluster: List[dict]) -> float:
    """small / (small + large) of the two largest main-peak areas (0..0.5)."""
    areas = sorted((peak_area(r) for r in cluster), reverse=True)
    if len(areas) < 2 or sum(areas[:2]) <= 0:
        return 0.0
    return min(areas[0], areas[1]) / (areas[0] + areas[1])


def cluster_confidence(cluster: List[dict], call: str) -> float:
    """Heuristic 0-100 confidence for *call* (see :data:`CONFIDENCE_NOTE`)."""
    if call == "no-call":
        return 0.0
    labels = _duplex_labels(cluster)
    if any(lab.startswith("HET") for lab in labels):
        return 80.0                          # an internal-standard HET is direct proof
    if any(lab.startswith("HOM") for lab in labels):
        return 80.0                          # matching a HOM standard is direct too
    areas = sorted((peak_area(r) for r in cluster), reverse=True)
    if len(areas) < 2 or areas[0] <= 0:
        return 0.0
    frac = minor_fraction(cluster)
    if call == "het":
        # 1.0 at a clean 50/50, 0.0 at the 15/85 boundary.
        boundary = min(abs(frac - GENOTYPE_HET_LO), abs(GENOTYPE_HET_HI - frac))
        half = (GENOTYPE_HET_HI - GENOTYPE_HET_LO) / 2.0
        distance_term = min(1.0, (boundary + half - abs(frac - 0.5)) / half) \
            if half > 0 else 0.0
        distance_term = max(0.0, distance_term)
    else:                                   # hom-major
        distance_term = min(1.0, max(0.0, (GENOTYPE_HET_LO - frac) / GENOTYPE_HET_LO))
    # both allele peaks must stand clear of the cluster floor
    floor = min(a for a in areas[:2])
    sep_term = min(1.0, floor / max(areas[0], 1e-9) / 0.5)
    conf = 70.0 * distance_term + 30.0 * sep_term
    return round(max(0.0, min(100.0, conf)), 1)


def score_rows(rows: List[dict]) -> List[dict]:
    """Return copies of *rows* with ``call``, ``minor_fraction`` and
    ``confidence`` filled in.

    Main peaks take their cluster's call.  +A satellite rows have no allelic
    area of their own, so they inherit the call/confidence of the nearest main
    peak in the same sample (they are a size marker, not a second allele)."""
    spacing = median_peak_spacing(rows)
    clusters = group_alleles(rows, spacing)
    by_id = {}
    for cl in clusters:
        call = genotype_call(cl)
        conf = cluster_confidence(cl, call)
        frac = minor_fraction(cl)
        for rec in cl:
            by_id[id(rec)] = (call, frac, conf)
    mains = [(r, _norm_key(r), _num(r.get("scan")))
             for r in rows if is_main(r)]

    def nearest_main(rec):
        key, scan = _norm_key(rec), _num(rec.get("scan"))
        best, best_d = None, None
        for m, mkey, mscan in mains:
            if mkey != key:
                continue
            d = abs(mscan - scan)
            if best_d is None or d < best_d:
                best, best_d = m, d
        return best if (best is not None and best_d is not None
                        and best_d <= allele_window(spacing)) else None

    def satellite_parent(rec):
        """The main peak a satellite shadows, from the picker's own record.

        A satellite sits up to ~1.6 spacings behind its parent, which is
        outside the allelic window ``nearest_main`` uses, so the search finds
        nothing and the ambiguity flag would be dropped.  The picker records
        the parent exactly, so prefer that and only fall back for rows
        exported by an older ``genotyping``.
        """
        parent = rec.get("parent_scan")
        if parent not in (None, ""):
            for m, mkey, mscan in mains:
                if (mkey == _norm_key(rec)
                        and mscan == _num(parent)):
                    return m
        return nearest_main(rec)

    out = []
    for r in rows:
        rec = dict(r)
        if is_main(r):
            call, frac, conf = by_id.get(id(r), ("no-call", _num(r.get("fraction")), 0.0))
            rec["call"], rec["confidence"] = call, conf
            rec["minor_fraction"] = round(frac, 4)
        else:
            donor = nearest_main(rec)
            call, frac, conf = by_id.get(id(donor), ("no-call", 0.0, 0.0))
            rec["call"] = "satellite-" + call if call != "no-call" else "satellite"
            rec["confidence"] = conf
            rec["minor_fraction"] = ""
        out.append(rec)

    # An ambiguous satellite cannot be scored on its own -- it was excluded from
    # the cluster -- so it downgrades the main peak it shadows.  Doing it here
    # rather than in genotype_call keeps the cluster functions pure.
    for rec, donor in ((r, satellite_parent(r)) for r in rows
                       if not is_main(r)):
        if rec.get("ambiguous") and donor is not None:
            for scored in out:
                if scored is rec:
                    continue
                if (is_main(scored) and _norm_key(scored) == _norm_key(donor)
                        and _num(scored.get("scan")) == _num(donor.get("scan"))):
                    if scored["call"] in ("het", "hom-major", "hom-1", "hom-2"):
                        scored["call"] = "uncertain"
                        scored["confidence"] = 0.0
                        scored["ambiguous_reason"] = \
                            "satellite may be a minor allele"
    return out


def well_summary(rows: List[dict]) -> dict:
    """One genotype record per (file, well) from its scored rows."""
    groups: dict = {}
    for r in rows:
        groups.setdefault((str(r.get("file", "")), str(r.get("well", ""))), []).append(r)
    out = []
    for (file, well), rs in sorted(groups.items()):
        scored = score_rows(rs)
        positions = [r for r in scored if is_main(r)]
        calls = [r["call"] for r in positions if r["call"] in ("het", "hom-major")]
        n_het = sum(1 for c in calls if c == "het")
        n_hom = sum(1 for c in calls if c == "hom-major")
        n_uncertain = sum(1 for r in positions if r["call"] == "uncertain")
        if n_uncertain and not (n_het or n_hom):
            overall = "uncertain"
        elif calls and n_het and n_hom:
            overall = "mixed"
        elif n_het:
            overall = "het"
        elif n_hom:
            overall = "hom-major"
        else:
            overall = "no-call"
        confs = [r["confidence"] for r in positions if r["call"] != "no-call"]
        clusters = group_alleles(rs, median_peak_spacing(rs))
        fracs = [minor_fraction(cl) for cl in clusters
                 if genotype_call(cl) in ("het", "hom-major")]
        out.append({
            "file": file, "well": well, "n_peaks": len(positions),
            "n_het": n_het, "n_hom": n_hom, "n_uncertain": n_uncertain,
            "genotype": overall,
            "mean_confidence": round(statistics.mean(confs), 1) if confs else 0.0,
            "mean_minor_fraction": round(statistics.mean(fracs), 4) if fracs else 0.0,
        })
    return out


# --------------------------------------------------------------------------- #
# rs1695 / CTC-CE positions
# --------------------------------------------------------------------------- #
# Heterozygosity rests on the *heteroduplex* peaks, not on the ratio of the two
# homoduplexes.  Heteroduplex only forms when two mismatched strands re-anneal,
# so its presence is direct proof of heterozygosity.  A homoduplex ratio cannot
# make that distinction: a severe allelic imbalance looks exactly like a
# homozygote.  Measured on the 96-well T9 rs1695 plate, gating on the ratio
# miscalls in both directions -- A12/G09/G10 are hets with big H2 that read as
# homo2, while E03/G03 are homo1 whose minor peak drags the ratio to ~0.82.
#
# The call therefore needs BOTH heteroduplexes clearly present and comparable.
# A one-sided heteroduplex is an artifact (stutter, or a partially re-annealed
# duplex) and is deliberately not read as evidence of heterozygosity -- B12 is
# a clean homo2 with a 119/24 sigma one-sided tail.
#
# T9_MIN_HET_SIGMA and T9_HET_RATIO_MIN are load-bearing.  They were fitted on
# the plate and G09 sits in a genuine dead zone: C07 outranks it on both
# measures (27.4/29.1 sigma, ratio 0.94) yet is a real homo2, so no cutoff
# separates the two.  Relaxing T9_MIN_HET_SIGMA to catch G09 costs 7 correctly
# called homozygotes.
T9_MIN_HET_SIGMA = 30.0        # each heteroduplex must clear this, in sigmas
T9_HET_RATIO_MIN = 0.5         # min/max of the two heteroduplexes
T9_MIN_DOMINANT_SIGMA = 40.0   # a well below this everywhere is no-call
# A second allele has to be a band you can see, not merely an area above zero.
# Deliberately the same 30 sigma the heteroduplexes must clear: one floor for
# "this band is real", applied to every duplex a het call rests on.
#
# Fitted on the rs1695 plate. All 44 reference hets clear 30 sigma on *both*
# homoduplexes -- the tightest is G10 at 46.0 / 1044.9, then A12 at 74.2 /
# 3300.5 -- so this costs no true call there. It fixes RS1695_N2 A04, which
# called het off heteroduplex peaks at 40 and 56 sigma while its first
# homoduplex read 7.5: 0.7% of the real peak, sitting in the very window a het
# needs a second allele in. It cannot create a het either, since a well with no
# heteroduplex evidence is rejected before this is reached.
T9_MIN_SECOND_ALLELE_SIGMA = 30.0
T9_MIN_TOTAL_AREA = 3000.0
T9_AI_DEVIATION = 0.25         # a het outside f1 in (0.25, 0.75) is AI


def t9_allele_fraction(hom1, hom2, het1, het2) -> float:
    """A (variant) fraction of one rs1695 position.

    Each heteroduplex is one A strand and one G strand, so half of its area
    belongs to A.  Returns 0.0 when the position carries no sample.
    """
    total = hom1 + hom2 + het1 + het2
    if total <= 0:
        return 0.0
    return (hom1 + 0.5 * (het1 + het2)) / total


def t9_call(hom1, hom2, het1, het2=None, sigmas=None) -> tuple:
    """Genotype for one rs1695 position -> ``(call, fraction, flags)``.

    *sigmas* is the optional ``(hom1, hom2, het1, het2)`` peak significance
    above the trace noise.  When given it drives the detection thresholds; when
    omitted the call rests on areas alone, which is all a hand-built peak table
    carries.

    ``call`` is one of ``CALLS``.  ``flags`` is a set and may contain ``"ai"``
    for a heterozygote far off 50/50.

    *het2* may be ``None``/omitted when the two heteroduplexes co-migrate and so
    measure as a single band.  That band carries the area of both heteroduplex
    strands together, so it is passed as *het1* and the het test falls back to
    requiring that one merged band to be present above the noise floor.  This is
    not a relaxation of the split-sample test, which is untouched: the separated
    case still requires two comparable heteroduplex bands, because a single band
    sitting where only one heteroduplex belongs is a shoulder, not a
    heterozygote.
    """
    # A merged heteroduplex is measured as one band holding both strands'
    # area.  Keep that as a distinct state (it is what lets the het test relax
    # below) while the arithmetic downstream still sees a plain number.
    merged_het = het2 is None
    het2 = 0.0 if merged_het else het2

    sigmas = tuple(x for x in (sigmas or ()) if x is not None)
    # For global dominance, only consider homoduplex sigmas - het bands shouldn't
    # make a weak homoduplex sample appear dominant (prevents het bands standing
    # in for missing homs).
    hom_sigmas = tuple(s for i, s in enumerate(sigmas) if i < 2 and s is not None)
    dom = max(hom_sigmas) if hom_sigmas else (max(sigmas) if sigmas else None)
    total = hom1 + hom2 + het1 + het2
    s1 = sigmas[0] if len(sigmas) > 0 else None
    s2 = sigmas[1] if len(sigmas) > 1 else None
    s3 = sigmas[2] if len(sigmas) > 2 else None
    s4 = sigmas[3] if len(sigmas) > 3 else None
    lo, hi = ((min(s3, s4), max(s3, s4))
              if (s3 is not None and s4 is not None) else (None, None))

    # A well with no usable sample is a no-call, not a homozygote.
    if dom is not None and dom < T9_MIN_DOMINANT_SIGMA:
        return "no-call", 0.0, set()
    if dom is None and total < T9_MIN_TOTAL_AREA:
        return "no-call", 0.0, set()

    frac = t9_allele_fraction(hom1, hom2, het1, het2)
    flags = set()

    # Het needs both heteroduplexes, present and comparable.  When they
    # co-migrate there is one merged band instead, which is equally conclusive
    # that both alleles are present -- that band cannot form from a single
    # allele, so its presence alone carries the evidence.
    if merged_het:
        both_present = (het1 > 0 and
                        (s3 is None or s3 >= T9_MIN_HET_SIGMA))
    elif s3 is not None:
        both_present = (lo >= T9_MIN_HET_SIGMA and hi > 0
                        and lo / hi >= T9_HET_RATIO_MIN)
    else:
        # No significance available: fall back to the two heteroduplex areas
        # being comparable to each other, which is the same idea on the scale a
        # peak table actually carries.
        both_present = het1 > 0 and het2 > 0 and min(het1, het2) / max(het1, het2) \
            >= T9_HET_RATIO_MIN

    # A heterozygote builds *both* homoduplexes -- one per allele -- so the
    # second homoduplex band has to carry area.  Without that requirement a het
    # can rest entirely on a band in the heteroduplex neighbourhood while no
    # second allele exists, which is exactly what the Taq A conformers do on
    # CYBA: they are a second sub-peak beside the one real allele, not a second
    # allele.  Measured across every het call carrying a reference genotype,
    # 0 of 31 genuine hets have an empty second homoduplex while 5 of 11 false
    # ones do, so this costs no true call and rejects half the false ones.
    # Also require at least one homoduplex to be significant when sigma info
    # is available; a strong het band should not create a het call when both
    # homs = below threshold.
    hom1_ok_sig = hom1 > 0 and (s1 is None or s1 >= T9_MIN_DOMINANT_SIGMA)
    hom2_ok_sig = hom2 > 0 and (s2 is None or s2 >= T9_MIN_DOMINANT_SIGMA)
    # Both homoduplexes must clear T9_MIN_SECOND_ALLELE_SIGMA.  A het's evidence
    # is the *minor* allele, so it arrives in the weaker of the two homoduplex
    # windows -- the one most likely to be a ripple rather than a product, which
    # is why the floor is applied to both and not to the dominant band.  Note
    # this is measured in significance, not area: the 1+2 shape is ordinary on
    # this assay (RS1695 A12 and G10 are reference hets with hom1 area 0), so
    # area cannot carry the test, while the bands themselves are plainly there.
    # It is deliberately ``hom*_ok_sig``'s sibling rather than those predicates,
    # which also require non-zero area.
    if both_present and hom2 > 0 and (
            s1 is None or s2 is None
            or (s1 >= T9_MIN_SECOND_ALLELE_SIGMA
                and s2 >= T9_MIN_SECOND_ALLELE_SIGMA)):
        if frac < T9_AI_DEVIATION or frac > 1.0 - T9_AI_DEVIATION:
            flags.add("ai")
        return "het", frac, flags

    # Homozygote: the labelled homoduplex says which allele, so unlike the
    # generic path this can distinguish the two homozygotes.
    #
    # A hom allele has to be *there*.  ``dom`` above only proves that *some*
    # band is strong, and the two are not the same claim: a strong
    # heteroduplex will happily carry a homozygous call that its own homoduplex
    # does not support.  ABCC2 H12 reads hom1 at 8.5 sigma against het1 at 44,
    # and was called hom-1 on the strength of a band belonging to a different
    # molecule.  And H01 read hom-2 from area2 25282 carrying only 16.9 sigma --
    # a segment-width artefact, since D07 has the same shape (area2 29354 at
    # 13.0 sigma) and is a hom-1 purely because its real allele happened to be
    # the larger area.  Ordering by area alone is not a substitute for asking
    # whether the allele is present.
    hom1_ok = hom1 > 0 and (s1 is None or s1 >= T9_MIN_DOMINANT_SIGMA)
    hom2_ok = hom2 > 0 and (s2 is None or s2 >= T9_MIN_DOMINANT_SIGMA)
    if hom1_ok and (not hom2_ok or hom1 >= hom2):
        return "hom-1", frac, flags
    if hom2_ok:
        return "hom-2", frac, flags
    if hom1_ok:
        return "hom-1", frac, flags
    return "no-call", 0.0, flags


# --------------------------------------------------------------------------- #
# read QC (needs the basecaller)
# --------------------------------------------------------------------------- #
def well_qc(doc, settings=None) -> dict:
    """Per-well read QC statistics from a (base-called) trace document."""
    from analyzer_core import AnalysisSettings, run_basecall
    settings = settings or AnalysisSettings()
    if not getattr(doc, "sequence", ""):
        run_basecall(doc, settings)
    quals = [float(q) for q in (doc.qualities or [])]
    seq = doc.sequence or ""
    peaks = [int(p) for p in (doc.peak_positions or [])]
    gaps = [b - a for a, b in zip(peaks, peaks[1:]) if b > a]
    acgt = getattr(doc, "acgt", None)
    chan_max = []
    if acgt is not None and len(acgt):
        chan_max = [round(float(acgt[:, i].max()), 4) for i in range(min(4, acgt.shape[1]))]
    qmean = statistics.mean(quals) if quals else 0.0
    qmin = min(quals) if quals else 0.0
    qmed = statistics.median(quals) if quals else 0.0
    q10 = sorted(quals)[max(0, len(quals) // 10 - 1)] if quals else 0.0
    # qualities are the caller's normalized peak heights, *not* Phred scores,
    # and their absolute scale differs per preset (raw_peaks ~ hundreds, the
    # spacing callers ~ single digits).  Only scale-independent measures are
    # safe for pass/fail, so we report the absolute values but gate on the
    # fraction of calls sitting below a quarter of the read's median signal.
    low_signal = (sum(1 for q in quals if q < 0.25 * qmed) / len(quals)) if quals else 1.0
    spacing_med = statistics.median(gaps) if gaps else 0.0
    spacing_cv = ((statistics.pstdev(gaps) / spacing_med) if len(gaps) > 1
                  and spacing_med > 0 else 0.0)
    n_frac = seq.count("N") / max(len(seq), 1)
    passed = bool(len(seq) >= MIN_READ_LEN and n_frac < MAX_N_FRACTION
                  and low_signal <= MAX_LOW_SIGNAL_FRAC)
    return {
        "well": getattr(doc, "well", Path(str(getattr(doc, "path", "?"))).stem),
        "source": getattr(doc, "source", "?"),
        "length": len(seq),
        "n_calls": len(peaks),
        "q_mean": round(qmean, 2),
        "q_min": round(qmin, 2),
        "q_10pct": round(q10, 2),
        "low_signal_frac": round(low_signal, 4),
        "n_bases": seq.count("N"),
        "mean_spacing": round(spacing_med, 2),
        "spacing_cv": round(spacing_cv, 4),
        "channel_max": chan_max,
        "preset": settings.basecaller,
        "pass": passed,
    }


# --------------------------------------------------------------------------- #
# HTML report
# --------------------------------------------------------------------------- #
def _table(headers: List[str], rows: List[List], bad_col: Optional[int] = None) -> str:
    head = "".join(f"<th>{html.escape(str(h))}</th>" for h in headers)
    body = []
    for r in rows:
        tds = []
        for i, c in enumerate(r):
            cls = ""
            if bad_col is not None and i == bad_col:
                cls = " class='pass'" if c is True else " class='fail'"
                c = "PASS" if c is True else "FAIL"
            tds.append(f"<td{cls}>{html.escape(str(c))}</td>")
        body.append("<tr>" + "".join(tds) + "</tr>")
    return (f"<table><thead><tr>{head}</tr></thead>"
            f"<tbody>{''.join(body)}</tbody></table>")


def html_report(qc_rows: List[dict] = (), summary_rows: List[dict] = (),
                title: str = "Limoncello scorer report") -> str:
    """Self-contained HTML report (no external assets) for QC and genotype rows."""
    parts = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        f"<title>{html.escape(title)}</title>",
        "<style>body{font-family:system-ui,Arial,sans-serif;margin:24px;"
        "color:#1a1c22}h1,h2{color:#22350e}table{border-collapse:collapse;"
        "margin:8px 0 24px;font-size:13px}th,td{border:1px solid #d7dbe0;"
        "padding:4px 8px;text-align:left}th{background:#f2f4f7}"
        ".pass{color:#0a7d33;font-weight:bold}.fail{color:#b3261e;"
        "font-weight:bold}.note{color:#555;font-size:12px}</style></head><body>",
        f"<h1>{html.escape(title)}</h1>",
        f"<p class='note'>{html.escape(CONFIDENCE_NOTE)}</p>",
    ]
    if qc_rows:
        parts.append("<h2>Read QC</h2>")
        headers = ["well", "source", "length", "q_mean", "q_min", "q_10pct",
                   "low_signal_frac", "n_bases", "mean_spacing", "spacing_cv",
                   "preset", "pass"]
        parts.append(_table(headers, [[r.get(h) for h in headers] for r in qc_rows],
                            bad_col=headers.index("pass")))
    if summary_rows:
        parts.append("<h2>Genotypes</h2>")
        headers = ["well", "n_peaks", "n_het", "n_hom", "genotype",
                   "mean_confidence", "mean_minor_fraction"]
        parts.append(_table(headers, [[r.get(h) for h in headers] for r in summary_rows]))
    if not qc_rows and not summary_rows:
        parts.append("<p>No rows.</p>")
    parts.append("</body></html>")
    return "\n".join(parts)


def sizing_html(results, title: str = "Limoncello fragment sizing report") -> str:
    """Self-contained HTML report for fragment-length sizing results.

    *results* is an iterable of ``fragment_sizing.SizingResult``.  Each well
    gets a fit summary (ladder, channels, anchors, leave-one-out error, any
    warnings) followed by its sized peaks.
    """
    parts = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        f"<title>{html.escape(title)}</title>",
        "<style>body{font-family:system-ui,Arial,sans-serif;margin:24px;"
        "color:#1a1c22}h1,h2{color:#22350e}table{border-collapse:collapse;"
        "margin:8px 0 24px;font-size:13px}th,td{border:1px solid #d7dbe0;"
        "padding:4px 8px;text-align:left}th{background:#f2f4f7}"
        ".note{color:#555;font-size:12px}.warn{color:#b3261e;font-size:12px}"
        "</style></head><body>",
        f"<h1>{html.escape(title)}</h1>",
    ]
    results = list(results)
    if not results:
        parts.append("<p>No results.</p>")
    for res in results:
        header = html.escape(res.well or res.file or "?")
        parts.append(f"<h2>Well {header}</h2>")
        q = res.quality
        err = q.get("rms_error_bp")
        err_txt = "n/a (need >=4 anchors)" if err is None else f"{err:.2f} bp"
        fit = [
            ["Ladder", f"{res.ladder.name} ({res.ladder.dye})"],
            ["Channels", f"ladder Ch{res.ladder_channel}, sample Ch{res.sample_channel}"
                         f" (order {res.base_order})"],
            ["Ladder anchors matched", len(res.anchors)],
            ["Leave-one-out RMS error", err_txt],
            ["Sample peaks", len(res.rows)],
        ]
        parts.append(_table(["field", "value"], fit))
        for warning in res.warnings:
            parts.append(f"<p class='warn'>! {html.escape(warning)}</p>")
        headers = ["scan", "length_bp", "height", "area", "in_range"]
        parts.append(_table(
            headers,
            [[r.get(h) for h in headers] for r in res.rows]))
    parts.append("<p class='note'>length_bp outside the ladder span is "
                 "extrapolated (in_range = False) and less trustworthy.</p>")
    parts.append("</body></html>")
    return "\n".join(parts)


# --------------------------------------------------------------------------- #
# persistence helpers
# --------------------------------------------------------------------------- #
def read_rows(path: str) -> List[dict]:
    """Read scored/picked rows from CSV, TSV or JSON."""
    p = Path(path)
    if p.suffix.lower() == ".json":
        data = json.loads(p.read_text(encoding="utf-8"))
        return [dict(r) for r in data]
    delim = "\t" if p.suffix.lower() in (".tsv", ".txt") else ","
    with open(p, newline="", encoding="utf-8-sig") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter=delim)]


def write_rows(path: str, rows: List[dict]) -> None:
    """Write rows as CSV or JSON (suffix-driven)."""
    p = Path(path)
    if not rows:
        p.write_text("[]\n" if p.suffix.lower() == ".json" else "",
                     encoding="utf-8")
        return
    if p.suffix.lower() == ".json":
        p.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
        return
    with open(p, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _cmd_check(_args) -> int:
    from analyzer_core import environment_report
    return environment_report()


def _cmd_qc(args) -> int:
    from analyzer_core import AnalysisSettings, load_trace
    settings = AnalysisSettings(basecaller=args.preset)
    rows = []
    for path in args.inputs:
        try:
            doc = load_trace(Path(path))
            rows.append(well_qc(doc, settings))
            print(f"{rows[-1]['well']:<10} {rows[-1]['length']:>5} bp  "
                  f"Qmean {rows[-1]['q_mean']:>5}  "
                  f"{'PASS' if rows[-1]['pass'] else 'FAIL'}")
        except Exception as exc:                            # noqa: BLE001
            print(f"{path}: {exc}", file=sys.stderr)
    if args.out:
        write_rows(args.out, rows)
        print(f"wrote {args.out} ({len(rows)} well(s))")
    if args.report:
        Path(args.report).write_text(html_report(qc_rows=rows), encoding="utf-8")
        print(f"wrote {args.report}")
    return 0 if rows else 1


def _cmd_peaks(args) -> int:
    rows = read_rows(args.input)
    scored = score_rows(rows)
    summary = well_summary(rows)
    if args.out:
        write_rows(args.out, scored)
        print(f"wrote {args.out} ({len(scored)} row(s))")
    if args.report:
        Path(args.report).write_text(
            html_report(summary_rows=summary), encoding="utf-8")
        print(f"wrote {args.report}")
    if not args.out and not args.report:
        for r in summary:
            print(f"{r['well']:<10} {r['genotype']:<10} "
                  f"conf {r['mean_confidence']:>5}  peaks {r['n_peaks']}")
    return 0


def _cmd_ladders(_args) -> int:
    from fragment_sizing import (BUILTIN_LADDERS, KNOWN_MEGABACE_STANDARDS,
                                 list_ladders)
    print("Bundled ladders (use by name with `size --ladder`):")
    for key in list_ladders():
        data = BUILTIN_LADDERS[key]
        lengths = data["lengths"]
        span = f"{min(lengths)}-{max(lengths)} bp, {len(lengths)} fragments"
        print(f"  {key:<18} {data['name']:<22} {data['dye']:<6} {span}")
    print("\nMegaBACE-compatible standards (fill in lengths from your kit insert):")
    for std in KNOWN_MEGABACE_STANDARDS:
        res = "" if std["resolution_bp"] is None else f", {std['resolution_bp']} bp resolution"
        print(f"  {std['name']}: {std['dye']}, {std['range_bp']}{res}")
    print("\nA custom ladder is a comma-separated length list or a JSON file:")
    print('  python scorer.py size well.rsd --lengths 50,100,150,200')
    print('  python scorer.py size well.rsd --ladder my_ladder.json')
    return 0


def _cmd_size(args) -> int:
    from analyzer_core import load_trace
    from fragment_sizing import load_ladder, size_trace
    try:
        ladder = load_ladder(args.lengths or args.ladder)
    except (ValueError, OSError, json.JSONDecodeError) as exc:  # noqa: BLE001
        print(f"ladder: {exc}", file=sys.stderr)
        return 2

    results = []
    for path in args.inputs:
        try:
            doc = load_trace(Path(path))
            res = size_trace(doc, ladder, ladder_channel=args.ladder_channel,
                             sample_channel=args.sample_channel,
                             base_order=args.base_order)
        except Exception as exc:                            # noqa: BLE001
            print(f"{path}: {exc}", file=sys.stderr)
            continue
        results.append(res)
        err = res.quality.get("rms_error_bp")
        err_txt = "n/a" if err is None else f"{err:.2f} bp"
        print(f"{res.well or res.file:<10} {len(res.rows):>3} peak(s)  "
              f"anchors {len(res.anchors):>2}  LOO RMS {err_txt}")
        for warning in res.warnings:
            print(f"    ! {warning}", file=sys.stderr)

    rows = [r for res in results for r in res.rows]
    if args.out:
        write_rows(args.out, rows)
        print(f"wrote {args.out} ({len(rows)} row(s))")
    if args.report:
        Path(args.report).write_text(sizing_html(results), encoding="utf-8")
        print(f"wrote {args.report}")
    return 0 if results else 1


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="limoncello-scorer",
                                 description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="command", required=True)

    q = sub.add_parser("qc", help="read QC for one or more trace files")
    q.add_argument("inputs", nargs="+", help=".rsd / .scf / .ab1 / text traces")
    q.add_argument("--preset", default="mb1000_accuracy",
                   help="basecaller preset (default mb1000_accuracy)")
    q.add_argument("--out", help="write the QC rows here (.csv or .json)")
    q.add_argument("--report", help="write a self-contained HTML report here")
    q.set_defaults(func=_cmd_qc)

    p = sub.add_parser("peaks", help="score exported picked-peak rows")
    p.add_argument("input", help="picked-peak table (.csv / .tsv / .json)")
    p.add_argument("--out", help="write the scored rows here (.csv or .json)")
    p.add_argument("--report", help="write an HTML genotype report here")
    p.set_defaults(func=_cmd_peaks)

    s = sub.add_parser("size", help="size fragment peaks against a ladder")
    s.add_argument("inputs", nargs="+", help=".rsd / .scf / .ab1 / text traces")
    s.add_argument("--ladder", default="genescan500_rox",
                   help="ladder name, JSON path or lengths (default genescan500_rox)")
    s.add_argument("--lengths",
                   help="comma-separated ladder lengths (overrides --ladder)")
    s.add_argument("--ladder-channel", type=int, default=4,
                   help="ladder channel 1..4 (default 4)")
    s.add_argument("--sample-channel", type=int, default=2,
                   help="sample channel 1..4 (default 2)")
    s.add_argument("--base-order", default="ACTG",
                   help="plate dye order, e.g. ACTG (default)")
    s.add_argument("--out", help="write sized peak rows here (.csv or .json)")
    s.add_argument("--report", help="write an HTML sizing report here")
    s.set_defaults(func=_cmd_size)

    ld = sub.add_parser("ladders", help="list ladders available to `size`")
    ld.set_defaults(func=_cmd_ladders)

    c = sub.add_parser("check", help="report which parts are available")
    c.set_defaults(func=_cmd_check)
    return ap


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
