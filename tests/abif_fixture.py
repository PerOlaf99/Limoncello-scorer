"""Write small ABIF (.fsa / .ab1) files for the loader and sizing tests.

There is no vendor .fsa in the repository, so the tests build their own to the
published ABIF layout (Applied Biosystems, "ABIF File Format Specification"):
a 128-byte header whose "tdir" entry points at a directory of 28-byte entries
(name, number, element type, element size, element count, data size, data
offset, handle); items of four bytes or fewer are stored inside the entry's
data-offset field.
"""
import struct

import numpy as np

SHORT, LONG, CHAR, PSTRING = 4, 5, 2, 18


def _item(name, number, etype, esize, payload, count):
    return (name.encode("ascii"), number, etype, esize, count, payload)


def short_array(name, number, values):
    vals = np.clip(np.round(np.asarray(values, float)), -32768, 32767).astype(">i2")
    return _item(name, number, SHORT, 2, vals.tobytes(), len(vals))


def pstring(name, number, text):
    b = text.encode("latin-1")
    return _item(name, number, PSTRING, 1, bytes([len(b)]) + b, len(b) + 1)


def chars(name, number, text):
    b = text.encode("latin-1")
    return _item(name, number, CHAR, 1, b, len(b))


def short(name, number, value):
    return _item(name, number, SHORT, 2, struct.pack(">h", value), 1)


def write_abif(path, items):
    """Write *items* (from the helpers above) as an ABIF file."""
    header_len = 128
    data = bytearray()
    entries = []
    for name, number, etype, esize, count, payload in items:
        if len(payload) <= 4:
            off_field = payload.ljust(4, b"\0")
            entries.append(struct.pack(">4sIHHII4sI", name, number, etype, esize,
                                       count, len(payload), off_field, 0))
        else:
            off = header_len + len(data)
            data += payload
            entries.append(struct.pack(">4sIHHIIII", name, number, etype, esize,
                                       count, len(payload), off, 0))
    dir_off = header_len + len(data)
    directory = b"".join(entries)
    head = b"ABIF" + struct.pack(">H", 101)
    head += struct.pack(">4sIHHIIII", b"tdir", 1, 1023, 28, len(entries),
                        len(directory), dir_off, 0)
    head = head.ljust(header_len, b"\0")
    with open(path, "wb") as fh:
        fh.write(head + bytes(data) + directory)


def gaussians(n, centres, height=1000.0, width=4.0, noise=3.0, seed=0):
    x = np.arange(n)
    y = np.random.default_rng(seed).normal(20.0, noise, n)
    for c in centres:
        y += height * np.exp(-0.5 * ((x - c) / width) ** 2)
    return y


def five_dye_fsa(path, sample_scans, ladder_scans, n=7500, well="B7",
                 dyes=("6-FAM", "VIC", "NED", "PET", "LIZ"), raw_too=True):
    """A 5-dye fragment run: the sample in channel 1 (6-FAM), the size
    standard in channel 5 (LIZ), the other channels near baseline."""
    analyzed = [gaussians(n, sample_scans, seed=1), gaussians(n, [], seed=2),
                gaussians(n, [], seed=3), gaussians(n, [], seed=4),
                gaussians(n, ladder_scans, height=900.0, seed=5)]
    items = []
    for k, num in enumerate((9, 10, 11, 12, 205)):
        items.append(short_array("DATA", num, analyzed[k]))
    if raw_too:   # raw data: same peaks, unseparated (a bit of everything)
        mix = np.mean(analyzed, axis=0)
        for k, num in enumerate((1, 2, 3, 4, 105)):
            items.append(short_array("DATA", num, analyzed[k] + 0.3 * mix))
    for k, dye in enumerate(dyes, start=1):
        items.append(pstring("DyeN", k, dye))
    items += [short("Dye#", 1, len(dyes)), chars("FWO_", 1, "GATC"),
              pstring("TUBE", 1, well), pstring("SMPL", 1, "sample 1")]
    write_abif(path, items)
