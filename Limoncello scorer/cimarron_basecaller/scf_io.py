"""
Reader/writer for SCF (Standard Chromatogram Format) files, version 3.

SCF is an open, publicly documented format for DNA sequencing chromatogram
data, specified at https://staden.sourceforge.net/scf-rfc.html (Dear &
Staden, "A standard file format for data from DNA sequencing instruments",
DNA Sequence 3, 107-110, 1992). It is NOT a proprietary GE/Amersham
format -- it's the standard interchange format the Cimarron software
itself can export to (see AutoBaseCall.exe's -SCF flag).

This module implements SCF version 3.00/3.10 as specified in the public
RFC: 128-byte header, per-channel sample blocks (delta-delta encoded),
and a Base structure per called base with peak index + per-channel
probability bytes.
"""

from __future__ import annotations
from dataclasses import dataclass, field
import struct
import numpy as np

SCF_MAGIC = 0x2E736366  # ".scf" packed as described in the RFC
HEADER_SIZE = 128
HEADER_FMT = ">4sIIIIIIII4sIIII18I"  # see _parse_header for field order
# Explicit big-endian ("forward byte, reverse bit" per RFC) struct format:
# magic(4s as 4 raw bytes, read as uint32 separately), samples, samples_offset,
# bases, bases_left_clip, bases_right_clip, bases_offset, comments_size,
# comments_offset, version(4s), sample_size, code_set, private_size,
# private_offset, spare[18]


@dataclass
class SCFHeader:
    samples: int
    samples_offset: int
    bases: int
    bases_offset: int
    comments_size: int
    comments_offset: int
    version: str
    sample_size: int
    code_set: int
    private_size: int
    private_offset: int


@dataclass
class SCFData:
    header: SCFHeader
    trace: np.ndarray          # (n_samples, 4) in A, C, G, T column order
    bases: str
    peak_indices: np.ndarray   # (n_bases,) int, index into `trace`
    probabilities: np.ndarray  # (n_bases, 4) uint8, per RFC prob_A..prob_T
    comments: dict


def _delta_decode(values: np.ndarray, modulus: int = 1 << 16) -> np.ndarray:
    """Reverse the SCF delta-delta encoding described in the RFC
    (two passes of cumulative sum, matching delta_samples2 with job=UN_DELTA_IT).

    Critically, the original C code operates on unsigned integers (uint16
    for 2-byte samples), so every addition implicitly wraps modulo 2^16.
    That wraparound must be replicated at each step here, or values that
    were encoded with a negative delta (stored wrapped-around) won't
    decode back correctly.
    """
    out = values.astype(np.int64).copy()
    for _ in range(2):
        running = 0
        for i in range(len(out)):
            out[i] = (out[i] + running) % modulus
            running = out[i]
    return out


def _delta_encode(values: np.ndarray, modulus: int = 1 << 16) -> np.ndarray:
    """Forward delta-delta encoding (job=DELTA_IT), inverse of _delta_decode.
    Deltas are taken modulo `modulus` to match the target unsigned integer
    width, exactly as the reference C implementation's implicit wraparound
    arithmetic would."""
    v = values.astype(np.int64).copy()
    n = len(v)
    out1 = np.empty(n, dtype=np.int64)
    prev = 0
    for i in range(n):
        cur = v[i]
        out1[i] = (cur - prev) % modulus
        prev = cur
    out2 = np.empty(n, dtype=np.int64)
    prev = 0
    for i in range(n):
        cur = out1[i]
        out2[i] = (cur - prev) % modulus
        prev = cur
    return out2


def read_scf(path: str) -> SCFData:
    """Parse an SCF v3.x file into an SCFData structure."""
    with open(path, "rb") as f:
        raw = f.read()

    if len(raw) < HEADER_SIZE:
        raise ValueError(f"{path}: too short to be a valid SCF file")

    magic = struct.unpack(">I", raw[0:4])[0]
    if magic != SCF_MAGIC:
        raise ValueError(
            f"{path}: bad SCF magic number {magic:#x} (expected {SCF_MAGIC:#x}); "
            "not an SCF file, or byte order is unexpected."
        )

    (
        samples, samples_offset, bases, bases_left_clip, bases_right_clip,
        bases_offset, comments_size, comments_offset,
    ) = struct.unpack(">8I", raw[4:36])
    version = raw[36:40].decode("ascii", errors="replace")
    sample_size, code_set, private_size, private_offset = struct.unpack(">4I", raw[40:56])
    # remaining 18 * 4 = 72 bytes are spare, ending at byte 128

    header = SCFHeader(
        samples=samples, samples_offset=samples_offset, bases=bases,
        bases_offset=bases_offset, comments_size=comments_size,
        comments_offset=comments_offset, version=version,
        sample_size=sample_size, code_set=code_set,
        private_size=private_size, private_offset=private_offset,
    )

    if version.startswith("3"):
        dtype = np.uint8 if sample_size == 1 else np.uint16
        itemsize = 1 if sample_size == 1 else 2

        # Version 3: four consecutive per-channel blocks (A, C, G, T), each
        # delta-delta encoded, each `samples` elements long.
        trace = np.zeros((samples, 4), dtype=np.int64)
        offset = samples_offset
        for ch in range(4):
            n_bytes = samples * itemsize
            chunk = raw[offset:offset + n_bytes]
            fmt = ">" + ("B" if sample_size == 1 else "H") * samples
            vals = np.array(struct.unpack(fmt, chunk), dtype=dtype)
            trace[:, ch] = _delta_decode(vals, modulus=1 << (8 * itemsize))
            offset += n_bytes

        # Bases section: peak_index[bases] (uint32), then probA/C/G/T[bases]
        # (uint8 each), then the called base chars, then 3*bases reserved bytes.
        off = bases_offset
        peak_indices = np.array(
            struct.unpack(f">{bases}I", raw[off:off + bases * 4]), dtype=np.int64
        )
        off += bases * 4

        probs = np.zeros((bases, 4), dtype=np.uint8)
        for ch in range(4):
            probs[:, ch] = np.frombuffer(raw[off:off + bases], dtype=np.uint8)
            off += bases

        base_chars = raw[off:off + bases].decode("ascii", errors="replace")

    elif version.startswith("1") or version.startswith("2"):
        # Versions 1/2: samples are interleaved (A,C,G,T) per point, NOT
        # delta-encoded (delta-delta encoding was introduced in v3). Sample
        # precision (sample_size) may be 1 or 2 bytes, per the header (or
        # assumed 1 byte if version < "2.00" -- not handled here since real
        # files in practice specify it).
        dtype = np.uint8 if sample_size == 1 else np.uint16
        itemsize = 1 if sample_size == 1 else 2
        n_vals = samples * 4
        fmt = ">" + ("B" if sample_size == 1 else "H") * n_vals
        n_bytes = n_vals * itemsize
        vals = np.array(
            struct.unpack(fmt, raw[samples_offset:samples_offset + n_bytes]),
            dtype=np.int64,
        )
        trace = vals.reshape(samples, 4)  # already in A, C, G, T column order

        # Bases section: `bases` flat 12-byte records:
        # peak_index(u32) prob_A prob_C prob_G prob_T(u8x4) base(char)
        # prob_sub prob_ins prob_del (u8x3)
        peak_indices = np.zeros(bases, dtype=np.int64)
        probs = np.zeros((bases, 4), dtype=np.uint8)
        base_chars_list = []
        off = bases_offset
        for i in range(bases):
            rec = raw[off:off + 12]
            peak_idx, pA, pC, pG, pT = struct.unpack(">IBBBB", rec[0:8])
            base_char = rec[8:9].decode("ascii", errors="replace")
            peak_indices[i] = peak_idx
            probs[i] = (pA, pC, pG, pT)
            base_chars_list.append(base_char)
            off += 12
        base_chars = "".join(base_chars_list)

    else:
        raise NotImplementedError(
            f"{path}: SCF version {version!r} detected; this reader implements "
            "versions 1.x, 2.x and 3.x."
        )

    comments_raw = raw[comments_offset:comments_offset + comments_size]
    comments = {}
    for line in comments_raw.decode("ascii", errors="replace").split("\n"):
        if "=" in line:
            k, _, v = line.partition("=")
            if k:
                comments[k] = v

    return SCFData(
        header=header, trace=trace.astype(np.float64), bases=base_chars,
        peak_indices=peak_indices, probabilities=probs, comments=comments,
    )


def write_scf(path: str, trace: np.ndarray, bases: str,
              peak_indices: np.ndarray, probabilities: np.ndarray,
              comments: dict | None = None, sample_size: int = 2) -> None:
    """Write an SCF v3.10 file. Mainly useful for round-trip testing this
    reader against the public spec, and for producing SCF files from this
    package's own basecall() output if you want a standard interchange
    format for downstream tools (Staden, BioPerl, sangerseqR, etc.)."""
    n_samples = trace.shape[0]
    n_bases = len(bases)
    comments = comments or {}

    dtype = np.uint8 if sample_size == 1 else np.uint16
    fmt_char = "B" if sample_size == 1 else "H"

    sample_blocks = b""
    for ch in range(4):
        encoded = _delta_encode(trace[:, ch].astype(np.int64))
        encoded = np.mod(encoded, 2 ** (8 * sample_size)).astype(dtype)
        sample_blocks += struct.pack(f">{n_samples}{fmt_char}", *encoded.tolist())

    bases_block = struct.pack(f">{n_bases}I", *peak_indices.astype(np.int64).tolist())
    for ch in range(4):
        bases_block += struct.pack(f">{n_bases}B", *probabilities[:, ch].astype(np.uint8).tolist())
    bases_block += bases.encode("ascii")
    bases_block += b"\x00" * (n_bases * 3)  # reserved sub/ins/del bytes

    comments_text = "\n".join(f"{k}={v}" for k, v in comments.items()) + "\n\0"
    comments_bytes = comments_text.encode("ascii")

    samples_offset = HEADER_SIZE
    bases_offset = samples_offset + len(sample_blocks)
    comments_offset = bases_offset + len(bases_block)

    header = struct.pack(">I", SCF_MAGIC)
    header += struct.pack(
        ">8I", n_samples, samples_offset, n_bases, 0, 0,
        bases_offset, len(comments_bytes), comments_offset,
    )
    header += b"3.10"
    header += struct.pack(">4I", sample_size, 0, 0, 0)
    header += struct.pack(">18I", *([0] * 18))
    assert len(header) == HEADER_SIZE, len(header)

    with open(path, "wb") as f:
        f.write(header)
        f.write(sample_blocks)
        f.write(bases_block)
        f.write(comments_bytes)
