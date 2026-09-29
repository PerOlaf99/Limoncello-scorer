"""
Reader for the binary MegaBACE .rsd (raw scan data) format.

This layout was determined empirically by comparing A01.rsd against a
matching known-good text export of the SAME sample's raw scan data
(A01.txt, tab-separated Scan/Channel1-4/Current columns) -- i.e. by using
the instrument software's own text export as ground truth to determine
the binary encoding. This is standard file-format interoperability work
(comparable to e.g. examining your own exported data to write a
compatible reader), not reverse-engineering any algorithm.

Validated: for A01.rsd (9647 scans), Channels 1-4 decode to an EXACT,
byte-for-byte match against the known text export for all 9647 rows.
The Current column differs by +/-1 on 75/9647 rows (0.8%), almost
certainly just rounding/smoothing applied somewhere in the text-export
path -- irrelevant to base calling, which only uses Channels 1-4.

Layout:
    bytes[0:4]      4-byte header (meaning not yet determined; doesn't
                    affect trace data -- skip it)
    bytes[4:-16]    N records of 5x little-endian int32:
                    [Channel1, Channel2, Channel3, Channel4, Current]
    bytes[-16:]     16-byte footer (meaning not yet determined)

Channel-to-base mapping is a per-run dye/instrument calibration setting
(the "Base order" field seen in text exports, e.g. "TGCA"), not encoded
in the .rsd binary itself as far as has been determined -- it's constant
across all wells of a given run/plate, so it only needs to be supplied
once per batch, not per file.
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np

RECORD_SIZE = 20   # 5 x int32
HEADER_SIZE = 4
FOOTER_SIZE = 16


@dataclass
class RSDFile:
    header_raw: bytes
    footer_raw: bytes
    trace: np.ndarray       # (n_scans, 4), Channel1..4 order
    current_raw: np.ndarray  # (n_scans,) raw int32 "current" column, uncorrected


def read_rsd(path: str, max_plausible_value: float = 2.0e5) -> RSDFile:
    """Parse a .rsd file.

    The trailing footer is NOT a fixed size (empirically it varies -- for
    A01.rsd the real scan table is exactly 9647 records, matching that
    file's own text-export metadata "Number of lines : 9647", while the
    file's total size implies 9708 records if a fixed 16-byte footer is
    assumed). Real channel values are always non-negative and modestly
    bounded (observed raw range ~0-20000); the footer/trailing region
    decodes to implausible values (negative, or absurdly large) as soon as
    you cross into it. So: decode every record the file's size allows,
    then truncate at the first record containing a negative value or a
    value above `max_plausible_value`.

    This is a heuristic boundary detector, not a fully reverse-engineered
    trailer structure -- if you have a matching text export with a
    "Number of lines" field for a given run, that's a more authoritative
    scan count to cross-check against (see raw_scan_csv.py).
    """
    with open(path, "rb") as f:
        data = f.read()

    n = len(data)
    body_len = n - HEADER_SIZE
    n_records_max = body_len // RECORD_SIZE
    if n_records_max <= 0:
        raise ValueError(f"{path}: file too short ({n} bytes) for the expected layout")

    header_raw = data[:HEADER_SIZE]
    ints = np.frombuffer(
        data[HEADER_SIZE:HEADER_SIZE + n_records_max * RECORD_SIZE], dtype="<i4"
    ).reshape(n_records_max, 5)

    channels = ints[:, 0:4]
    implausible = (channels < 0) | (channels > max_plausible_value)
    bad_rows = np.where(implausible.any(axis=1))[0]
    n_real = int(bad_rows[0]) if len(bad_rows) else n_records_max

    footer_raw = data[HEADER_SIZE + n_real * RECORD_SIZE:]
    trace = ints[:n_real, 0:4].astype(np.float64)
    current_raw = ints[:n_real, 4].astype(np.float64)

    return RSDFile(header_raw=header_raw, footer_raw=footer_raw, trace=trace, current_raw=current_raw)


def to_acgt_trace(rsd: RSDFile, base_order: str) -> tuple[np.ndarray, str]:
    """Reorder rsd.trace's 4 columns (Channel1..4) into standard A,C,G,T
    column order, given the run's base_order string (e.g. "TGCA" means
    Channel1=T, Channel2=G, Channel3=C, Channel4=A).

    base_order must be supplied by the caller -- it's a per-run dye/
    instrument calibration setting, not stored in the .rsd binary itself
    (as determined so far); check the plate's other exported files (CSV/
    TXT export, or Basecall.ini-adjacent run metadata) for this run's value.
    """
    channel_to_base = list(base_order.upper())
    if sorted(channel_to_base) != ["A", "C", "G", "T"]:
        raise ValueError(f"Unexpected base_order {base_order!r}; expected a permutation of ACGT")

    base_to_channel_idx = {b: i for i, b in enumerate(channel_to_base)}
    out = np.empty_like(rsd.trace)
    for target_pos, base in enumerate("ACGT"):
        out[:, target_pos] = rsd.trace[:, base_to_channel_idx[base]]
    return out, "ACGT"
