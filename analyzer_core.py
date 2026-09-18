"""
Core data loading and processing for the Limoncello CE Analyzer.
Wraps best_basecaller configs and optional Cimarron-style DSP stages.
"""
from __future__ import annotations

import struct
import sys
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import Optional, List, Dict, Any, Tuple

import numpy as np

# Locate best_basecaller package
_ROOT = Path(__file__).resolve().parent
for _p in (
    _ROOT.parent / "BEST_BASECALLER_RELEASE",
    Path("/home/per/Nedlastinger/BEST_BASECALLER_RELEASE"),
    _ROOT.parent / "best_basecaller" / "best_basecaller",
    Path("/home/workdir/artifacts/BEST_BASECALLER_RELEASE"),
    Path("/home/workdir/artifacts/best_basecaller/best_basecaller"),
):
    if (_p / "cimarron_basecaller").is_dir() or (_p / "cimarron_basecaller").exists():
        sys.path.insert(0, str(_p))
        break

try:
    from cimarron_basecaller.rsd_io import read_rsd, to_acgt_trace
except ImportError:
    read_rsd = to_acgt_trace = None  # type: ignore

try:
    from cimarron_basecaller.scf_io import read_scf
except ImportError:
    read_scf = None  # type: ignore

try:
    from cimarron_basecaller.spacing_caller import track_bases
except ImportError:
    track_bases = None  # type: ignore

try:
    from configs import CONFIGS as BB_CONFIGS
except ImportError:
    BB_CONFIGS = {}


# ---------------------------------------------------------------------------
# Analysis parameter set (Sequence Analyzer "knobs")
# ---------------------------------------------------------------------------
@dataclass
class AnalysisSettings:
    """Mirrors typical CE Sequence Analyzer / basecall settings."""

    # Basecaller version
    basecaller: str = "pos_bonus07"  # pos_bonus07 | pos_profile | hz_soften | raw_peaks

    # Channel / dye
    base_order: str = "TGCA"  # instrument dye order → ACGT columns

    # Baseline
    baseline_method: str = "percentile"  # percentile | none
    baseline_window: int = 151

    # Smoothing
    smooth_enable: bool = True
    smooth_window: int = 5  # odd Savitzky-Golay style window if used later

    # Spectral separation
    spectral_enable: bool = True
    position_adaptive_spectral: bool = True

    # Mobility
    mobility_enable: bool = True

    # Band filter / deconv
    use_gaussian_reconstruction: bool = True
    gaussian_recon_segment_size: int = 384
    gaussian_recon_noise_reg: float = 0.05
    use_multipass_wiener: bool = False

    # Spacing tracker
    use_combined_channel_score: bool = True
    channel_peak_bonus: float = 0.7
    pullback_weight: float = 0.008
    ema_alpha: float = 0.08
    window_frac_lo: float = 0.70
    window_frac_hi: float = 1.30
    local_norm_window: int = 1800
    pullback_profile_enable: bool = True
    pullback_frac: float = 0.33
    pullback_start: float = 0.008
    pullback_end: float = 0.001

    # Display / signal region
    signal_start: int = 0
    signal_end: int = 0  # 0 = full

    # View
    view_mode: str = "processed"  # raw | baseline | processed | called

    def to_track_kwargs(self) -> Dict[str, Any]:
        """Map settings → track_bases kwargs (best_basecaller)."""
        kw: Dict[str, Any] = dict(
            use_gaussian_reconstruction=self.use_gaussian_reconstruction,
            gaussian_recon_segment_size=self.gaussian_recon_segment_size,
            gaussian_recon_noise_reg=self.gaussian_recon_noise_reg,
            use_multipass_wiener=self.use_multipass_wiener,
            use_combined_channel_score=self.use_combined_channel_score,
            channel_peak_bonus=self.channel_peak_bonus,
            pullback_weight=self.pullback_weight,
            ema_alpha=self.ema_alpha,
            window_frac=(self.window_frac_lo, self.window_frac_hi),
            local_norm_window=self.local_norm_window,
            baseline_window=self.baseline_window,
            position_adaptive_spectral=self.position_adaptive_spectral and self.spectral_enable,
            local_hardzone_deconv=False,
        )
        if self.pullback_profile_enable:
            kw["pullback_profile"] = (
                self.pullback_frac,
                self.pullback_start,
                self.pullback_end,
            )
        # Merge named config defaults if present
        if self.basecaller in BB_CONFIGS:
            base = dict(BB_CONFIGS[self.basecaller])
            base.update({k: v for k, v in kw.items() if v is not None})
            # ensure flags
            base.setdefault("local_hardzone_deconv", False)
            base.setdefault("use_multipass_wiener", self.use_multipass_wiener)
            return base
        return kw


BASECALLER_VERSIONS = {
    "pos_bonus07": "Best dual-aware (recommended)",
    "pos_profile": "Max matched_bp (longer tail)",
    "hz_soften": "Mild mid-zone less deconv",
    "raw_peaks": "Minimal processing + envelope peaks",
}


@dataclass
class TraceDocument:
    """One loaded well / file."""

    path: Path
    well: str
    raw: np.ndarray  # (n, 4) instrument order before map
    base_order: str
    acgt: np.ndarray  # (n, 4) A,C,G,T
    current: Optional[np.ndarray] = None  # raw current (µA x 10 for RSD)
    sequence: str = ""
    peak_positions: List[int] = field(default_factory=list)
    qualities: List[float] = field(default_factory=list)
    settings_used: Optional[AnalysisSettings] = None
    source: str = "rsd"  # rsd | abi | scf | text
    meta: str = ""       # free-form instrument/source info for exports

    @property
    def n_scans(self) -> int:
        return int(self.acgt.shape[0])

    @property
    def current_ua(self) -> Optional[np.ndarray]:
        """RSD 'current' column is stored in tenths of µA -> divide by 10.
        Spurious tail samples are clamped to robust percentiles."""
        if self.current is None:
            return None
        c = np.asarray(self.current, dtype=float) / 10.0
        if c.size:
            lo = float(np.nanpercentile(c, 0.5))
            hi = float(np.nanpercentile(c, 99.9))
            c = np.clip(c, lo, hi)
        return c


RFD_EXTS = (".rsd", ".RSD")
SCF_EXTS = (".scf", ".SCF")
ABI_EXTS = (".abi", ".ab1")
TEXT_EXTS = (".txt", ".tsv", ".csv", ".dat")
SUPPORTED_EXTS = RFD_EXTS + SCF_EXTS + ABI_EXTS + TEXT_EXTS
# Only instrument trace formats are auto-discovered in folders; text/CSV
# exports (Text/, *.txt, reports) are intentionally NOT listed.
DISCOVER_EXTS = RFD_EXTS + SCF_EXTS + ABI_EXTS

# Instrument-folder convention for which wells to show per data type (informational)
SOURCE_LABELS = {
    "rsd": "RSD (raw + current)",
    "scf": "SCF (std.)",
    "abi": "ABI / ThermoFisher (.ab1)",
    "text": "Text / CSV trace",
}


def _isnum(token: str) -> bool:
    try:
        float(token)
        return True
    except ValueError:
        return False


# Instrument logs that are numeric tables but not traces
TEXT_IGNORE_STEMS = ("cderr", "totalerr", "rferr", "runerr", "matrixerr")


def is_likely_trace_file(path: Path) -> bool:
    """Light filter so logs/reports (CDErr.txt, TotalErr.txt, *_report.txt)
    never pollute the file list. Known instrument log names are excluded
    outright; any other text file must contain at least one all-numeric row
    (the same rule the loader enforces on the full file)."""
    ext = path.suffix.lower()
    if ext not in {e.lower() for e in TEXT_EXTS}:
        return True
    if path.stem.lower().lstrip("_ ").startswith(TEXT_IGNORE_STEMS):
        return False
    try:
        with open(path, "r", errors="ignore") as f:
            for i, ln in enumerate(f):
                if i > 255:
                    break
                toks = ln.replace(",", " ").split()
                if toks and all(_isnum(t) for t in toks):
                    return True
            return False
    except OSError:
        return False


def discover_files(folders: List[Path]) -> List[Path]:
    files: List[Path] = []
    for folder in folders:
        if not folder.is_dir():
            continue
        for ext in DISCOVER_EXTS:
            files.extend(sorted(folder.glob(f"*{ext}")))
    # unique, preserve sorted order per folder
    seen = set()
    out = []
    for f in files:
        k = str(f.resolve())
        if k not in seen:
            seen.add(k)
            out.append(f)
    return out


def find_run_folders(root: Path, max_depth: int = 6) -> List[Path]:
    """Return the data folders under ``root`` (each directly containing trace
    files). A folder that directly holds trace files is itself a run folder;
    otherwise its subfolders are scanned (bounded by ``max_depth``). Adding a
    big project folder like an OY run collection therefore loads every run."""
    ext_lower = {e.lower() for e in DISCOVER_EXTS}
    out: List[Path] = []

    def walk(d: Path, depth: int):
        if depth > max_depth:
            return
        try:
            entries = sorted(d.iterdir(), key=lambda p: p.name)
        except OSError:
            return
        has_data = any(e.is_file() and e.suffix.lower() in ext_lower
                       for e in entries)
        if has_data:
            out.append(d)
            return  # data folder — don't descend further
        for e in entries:
            if e.is_dir() and not e.name.startswith("."):
                walk(e, depth + 1)

    walk(root, 0)
    return out


def discover_rsd(folders: List[Path]) -> List[Path]:
    """Legacy alias — RSD only."""
    files: List[Path] = []
    for folder in folders:
        if not folder.is_dir():
            continue
        for ext in RFD_EXTS:
            files.extend(sorted(folder.glob(f"*{ext}")))
    seen = set()
    out = []
    for f in files:
        k = str(f.resolve())
        if k not in seen:
            seen.add(k)
            out.append(f)
    return out


def _reorder(ch: np.ndarray, base_order: str) -> np.ndarray:
    """Reorder columns from instrument order string -> A,C,G,T."""
    order = base_order.upper()
    if sorted(order) != ["A", "C", "G", "T"]:
        raise ValueError(f"Unexpected base_order {base_order!r}")
    base_to_col = {b: i for i, b in enumerate(order)}
    out = np.empty_like(ch)
    for target, base in enumerate("ACGT"):
        out[:, target] = ch[:, base_to_col[base]]
    return out


# ---------------------------------------------------------------------------
# ABIF (.ab1 / .abi from Applied Biosystems / ThermoFisher) reader
# ---------------------------------------------------------------------------
def _abif_entries(raw: bytes) -> Dict[str, Tuple[int, int, int, int, int]]:
    """Parse an ABIF directory.  Standard ABI layout first, then the
    alternate .abd variant (numElements at 16:20), whichever yields records."""
    candidates = [
        (6, 26),   # standard ABIF
        (16, 26),  # alternate .abd variant
    ]
    for num_off, dir_off in candidates:
        try:
            doff = struct.unpack(">I", raw[dir_off:dir_off + 4])[0]
            nelem = struct.unpack(">I", raw[num_off:num_off + 4])[0]
        except Exception:
            continue
        entries: Dict[str, Tuple[int, int, int, int, int]] = {}
        for i in range(min(nelem, 500)):
            off = doff + i * 28
            if off + 28 > len(raw):
                break
            name = raw[off:off + 4].decode("ascii", "replace")
            try:
                num, etype, esize, ne, dsize, d_off, _ = struct.unpack(
                    ">IHHIIII", raw[off + 4:off + 28]
                )
            except Exception:
                break
            if not name.isidentifier():
                break
            if 0 <= d_off and d_off + dsize <= len(raw) and dsize > 0:
                entries[f"{name}{num}"] = (etype, esize, ne, dsize, d_off)
        if any(k.startswith("DATA") for k in entries):
            return entries
    return {}


def _abif_payload(raw: bytes, rec: Tuple[int, int, int, int, int]) -> np.ndarray:
    etype, esize, ne, dsize, doff = rec
    dtype = {2: ">i1", 4: ">i2", 5: ">i4", 6: ">i8", 7: ">f4", 8: ">f8"}.get(etype)
    buf = raw[doff:doff + dsize]
    if dtype is None:
        return np.frombuffer(buf, dtype=">i2").astype(np.float64)
    return np.frombuffer(buf, dtype=dtype).astype(np.float64)


def load_ab1(path: Path, base_order: str = "ACGT") -> TraceDocument:
    raw = path.read_bytes()
    if raw[:4] != b"ABIF":
        raise ValueError(f"{path.name}: not an ABIF file")
    entries = _abif_entries(raw)
    if not any(k.startswith("DATA") for k in entries):
        raise ValueError(f"{path.name}: no ABIF DATA records found")

    ch = np.column_stack([_abif_payload(raw, entries[f"DATA{i}"]) for i in range(1, 5)])

    # FWO_1 holds the run's dye/base-order string (channel c -> base).
    order = base_order.upper()
    for key in ("FWO_1", "FWO_2"):
        rec = entries.get(key)
        if rec:
            s = bytes(raw[rec[4]:rec[4] + rec[3]]).decode("ascii", "replace").strip()
            if len(s) == 4 and sorted(s.upper()) == ["A", "C", "G", "T"]:
                order = s.upper()
                break

    # Peak locations + called bases (PBAS/PLOC).
    seq, peaks = "", []
    for pk_key, pb_key in (("PLOC_1", "PBAS_1"), ("PLOC_2", "PBAS_2")):
        if pk_key in entries and pb_key in entries:
            p = _abif_payload(raw, entries[pk_key]).astype(np.int64)
            pr = entries[pb_key]
            seq = bytes(raw[pr[4]:pr[4] + pr[3]]).decode("ascii", "replace")
            peaks = list(p)
            break

    acgt = _reorder(ch, order)
    return TraceDocument(
        path=path,
        well=path.stem,
        raw=acgt.copy(),
        base_order=order,
        acgt=acgt,
        sequence=seq,
        peak_positions=peaks,
        source="abi",
        meta=f"ABIF .ab1 (ABI/ThermoFisher), dye-order {order}",
    )


def load_scf_trace(path: Path, base_order: str = "ACGT") -> TraceDocument:
    if read_scf is None:
        raise RuntimeError("scf_io (best_basecaller) not available for SCF files")
    scf = read_scf(str(path))
    tr = np.asarray(scf.trace, dtype=float)
    peaks = list(np.asarray(scf.peak_indices, dtype=np.int64)) if scf.peak_indices is not None and len(scf.peak_indices) else []
    seq = ""
    if getattr(scf, "bases", None):
        seq = "".join(ch for ch in scf.bases if ch in "ACGTN")
    return TraceDocument(
        path=path,
        well=path.stem,
        raw=tr.copy(),
        base_order="ACGT",
        acgt=tr,
        sequence=seq,
        peak_positions=peaks,
        source="scf",
        meta="SCF (standard / Agilent-style)",
    )


def load_text_trace(path: Path, base_order: str = "ACGT") -> TraceDocument:
    """Generic numeric trace: 4 cols = A,C,G,T; >=5 cols = scan + channels;
    >=6 cols = scan + channels + current (µA)."""
    txt = path.read_text(encoding="utf-8", errors="replace")
    rows: List[List[float]] = []
    header_lines: List[str] = []
    for ln in txt.splitlines():
        if not ln.strip():
            continue
        stripped = ln.lstrip()
        if stripped.startswith(("#", ";", "//")):
            header_lines.append(ln.strip())
            continue
        nums: List[float] = []
        for tok in ln.replace(",", " ").split():
            try:
                nums.append(float(tok))
            except ValueError:
                continue
        if nums:
            rows.append(nums)
    if not rows:
        raise ValueError(f"{path.name}: no numeric rows found")
    w = max(len(r) for r in rows)
    arr = np.array([r + [0.0] * (w - len(r)) for r in rows], dtype=float)
    if arr.shape[1] >= 6:
        ch = arr[:, 1:5]
        cur = arr[:, 5]
    elif arr.shape[1] == 5:
        ch = arr[:, 1:5]
        cur = None
    elif arr.shape[1] == 4:
        ch = arr[:, :4]
        cur = None
    else:
        raise ValueError(f"{path.name}: need >=4 numeric columns, got {arr.shape[1]}")
    order = base_order.upper()
    acgt = _reorder(ch, order)
    meta = "; ".join(header_lines[:3])
    return TraceDocument(
        path=path,
        well=path.stem,
        raw=ch.copy(),
        base_order=order,
        acgt=acgt,
        current=None if cur is None else np.asarray(cur, float),
        source="text",
        meta=meta or "Text/CSV trace",
    )


def load_trace(path: Path, base_order: str = "TGCA") -> TraceDocument:
    """Load any supported trace format into a TraceDocument.

    base_order default is the RSD dye order ('TGCA'); non-RSD
    formats carry their own order where applicable."""
    suffix = path.suffix.lower()
    try:
        if suffix in {e.lower() for e in RFD_EXTS}:
            return load_rsd(path, base_order)
        if suffix in {e.lower() for e in SCF_EXTS}:
            return load_scf_trace(path)
        if suffix in {e.lower() for e in ABI_EXTS}:
            return load_ab1(path, base_order)
        if suffix in {e.lower() for e in TEXT_EXTS}:
            return load_text_trace(path, base_order)
    except ValueError:
        raise
    raise ValueError(f"{path.name}: unsupported format {suffix}")


def load_rsd(path: Path, base_order: str = "TGCA") -> TraceDocument:
    if read_rsd is None or to_acgt_trace is None:
        raise RuntimeError(
            "cimarron_basecaller not found. Place BEST_BASECALLER_RELEASE next to this app."
        )
    rsd = read_rsd(str(path))
    # rsd may be dict or object
    if isinstance(rsd, dict):
        ch = np.column_stack([rsd[f"ch{i}"] for i in range(1, 5)]) if "ch1" in rsd else rsd.get("channels")
        cur = rsd.get("current")
        acgt, order = to_acgt_trace(ch, base_order=base_order)
        return TraceDocument(
            path=path,
            well=path.stem,
            raw=np.asarray(ch, float),
            base_order=order,
            acgt=np.asarray(acgt, float),
            current=None if cur is None else np.asarray(cur, float),
            source="rsd",
            meta="RSD",
        )

    try:
        acgt, order = to_acgt_trace(rsd, base_order=base_order)
    except Exception:
        acgt, order = to_acgt_trace(rsd.trace, base_order=base_order)
    cur = getattr(rsd, "current_raw", None)
    if cur is None:
        cur = getattr(rsd, "current", None)
    meta = "RSD"
    for block, label in (("header_raw", "header"), ("footer_raw", "footer")):
        b = getattr(rsd, block, None)
        if b:
            try:
                meta += f" | {label}: {b[:28]!r}"
            except Exception:
                pass
    return TraceDocument(
        path=path,
        well=path.stem,
        raw=acgt.copy(),
        base_order=order,
        acgt=acgt,
        current=None if cur is None else np.asarray(cur, float),
        source="rsd",
        meta=meta,
    )


def run_basecall(doc: TraceDocument, settings: AnalysisSettings) -> TraceDocument:
    """Run selected basecaller; updates sequence + peak_positions."""
    if track_bases is None:
        raise RuntimeError("track_bases not available")

    order = "ACGT"
    kw = settings.to_track_kwargs()
    # mobility / spectral off if disabled
    if not settings.mobility_enable:
        kw["mobility_shifts"] = (0, 0, 0, 0)
    if not settings.spectral_enable:
        kw["position_adaptive_spectral"] = False
        # identity matrix if supported
        try:
            kw["spectral_separation_matrix"] = np.eye(4)
        except Exception:
            pass

    if settings.basecaller == "raw_peaks":
        # minimal: local maxima on envelope
        env = doc.acgt.max(axis=1)
        peaks = []
        for i in range(2, len(env) - 2):
            if env[i] >= env[i - 1] and env[i] > env[i + 1] and env[i] > 0.05 * env.max():
                if not peaks or i - peaks[-1] >= 5:
                    peaks.append(i)
        seq = "".join(order[int(np.argmax(doc.acgt[p]))] for p in peaks)
        doc.sequence = seq
        doc.peak_positions = peaks
        doc.qualities = [float(env[p]) for p in peaks]
        doc.settings_used = settings
        return doc

    seq, quals, bands = track_bases(doc.acgt, base_order=order, **kw)
    doc.sequence = seq
    doc.peak_positions = [int(b.position) for b in bands]
    doc.qualities = list(quals) if quals is not None else []
    doc.settings_used = settings
    return doc


def display_trace(doc: TraceDocument, settings: AnalysisSettings) -> np.ndarray:
    """Return (n,4) array for plotting according to view_mode."""
    if settings.view_mode == "raw":
        return doc.raw
    # For baseline/processed without full intermediate API, show acgt
    # (full pipeline intermediates would need deeper hooks)
    return doc.acgt
