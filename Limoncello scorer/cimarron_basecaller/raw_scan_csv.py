"""
Reader for the tab-separated "raw scan" export of a MegaBACE .rsd file
(the format produced by the instrument software's own CSV/TXT export --
not a reverse-engineered binary format, so this is on very solid ground).

File structure observed:
  - ~35 lines of "Key : Value" metadata (run info, chemistry, base order,
    base caller used, etc.), separated by blank lines into sections
  - A column header line: "Scan  Channel1  Channel2  Channel3  Channel4  Current"
  - Tab-separated data rows: scan_index, ch1, ch2, ch3, ch4, current

Notably the metadata includes "Base order : TGCA" -- i.e. which base each
of Channel1..Channel4 corresponds to for THIS run (the dye/channel
assignment is a per-run calibration choice, not a fixed constant).
"""

from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np


@dataclass
class RawScanFile:
    metadata: dict
    base_order: str          # e.g. "TGCA" -> channel1=T, channel2=G, channel3=C, channel4=A
    trace: np.ndarray        # (n_scans, 4) float, columns in Channel1..4 order (NOT yet reordered to ACGT)
    current_ua: np.ndarray   # (n_scans,) float, current in microamps (corrected)
    scan_rate_hz: float = 1.75  # per your note; not always present/reliable in the header


def read_raw_scan_csv(path: str, scan_rate_hz: float = 1.75, current_scale: float = 0.1) -> RawScanFile:
    """Parse a MegaBACE raw-scan text/CSV export.

    current_scale: correction factor for the Current column (per your note
    that raw values like 138 should be read as 13.8 uA -- i.e. divide by 10).
    """
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    metadata = {}
    base_order = None
    header_line_idx = None

    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.lower().startswith("scan") and "channel1" in stripped.lower():
            header_line_idx = i
            break
        if ":" in stripped and not stripped.startswith("\t"):
            key, _, val = stripped.partition(":")
            key = key.strip()
            val = val.strip()
            if key:
                metadata[key] = val
                if key.lower() == "base order":
                    base_order = val

    if header_line_idx is None:
        raise ValueError(f"{path}: could not find 'Scan  Channel1 ...' header row")
    if base_order is None:
        raise ValueError(
            f"{path}: no 'Base order' field found in metadata -- can't determine "
            "which channel is which base without it."
        )

    data_rows = []
    for line in lines[header_line_idx + 1:]:
        parts = [p for p in line.strip().split("\t") if p != ""]
        if len(parts) < 6:
            continue
        try:
            values = [float(x) for x in parts[:6]]
        except ValueError:
            continue
        data_rows.append(values)

    arr = np.array(data_rows)
    if arr.shape[0] == 0:
        raise ValueError(f"{path}: no data rows parsed")

    scan_idx = arr[:, 0]
    trace = arr[:, 1:5]
    current = arr[:, 5] * current_scale

    # Sanity check: scan index should be a simple 0..n-1 running count.
    expected = np.arange(len(scan_idx))
    if not np.array_equal(scan_idx, expected):
        raise ValueError(
            f"{path}: scan index column isn't a clean 0..n-1 sequence -- "
            "file may be malformed or use a different layout than expected."
        )

    return RawScanFile(
        metadata=metadata, base_order=base_order, trace=trace,
        current_ua=current, scan_rate_hz=scan_rate_hz,
    )


def to_acgt_trace(raw: RawScanFile) -> tuple[np.ndarray, str]:
    """Reorder raw.trace's 4 columns from Channel1..4 order into standard
    A, C, G, T column order, using raw.base_order (e.g. "TGCA" means
    channel index 0=T, 1=G, 2=C, 3=A).

    Returns (trace_acgt, "ACGT") ready for cimarron_basecaller.basecall().
    """
    channel_to_base = list(raw.base_order.upper())
    if sorted(channel_to_base) != ["A", "C", "G", "T"]:
        raise ValueError(f"Unexpected base_order {raw.base_order!r}; expected a permutation of ACGT")

    base_to_channel_idx = {b: i for i, b in enumerate(channel_to_base)}
    out = np.empty_like(raw.trace)
    for target_pos, base in enumerate("ACGT"):
        out[:, target_pos] = raw.trace[:, base_to_channel_idx[base]]
    return out, "ACGT"
