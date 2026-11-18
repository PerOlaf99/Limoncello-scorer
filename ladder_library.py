"""The size-standard (ladder) library: bundled ladders plus the user's own.

Two kinds of ladder:

* **Bundled** -- commercial standards whose band table comes from a primary
  source (kit insert, vendor software or datasheet).  They live in
  :data:`fragment_sizing.BUILTIN_LADDERS`, cite that source, and cannot be
  edited or deleted here.  A half-remembered table sizes every peak wrongly,
  so nothing is bundled without one.
* **User ladders** -- custom or imported ones, one JSON file each under
  ``<mark library>/ladders/``, beside the marks and preferences, so they move
  with the tool and survive updates::

      {"name": "ET550-R (lot 1234)", "vendor": "GE", "dye": "ET-ROX",
       "color": "#d62728", "lengths": [50, 75, ...], "exclude": [],
       "default_channel": 1, "source": "kit insert", "notes": ""}

Only ``lengths`` matters to sizing.  ``exclude`` lists bands that are in the
lane but must not anchor the curve; ``color`` is for display; and
``default_channel`` pre-fills the ladder channel (the actual channel is still
chosen per run -- the same kit can sit on another channel elsewhere).

Import (:func:`import_ladder_file`) reads a ladder JSON (one ladder, a list of
them, or ``{"ladders": [...]}``), and plain text / CSV / TSV -- including a
MegaBACE Genetic Profiler ``.tdb`` table: a row of four or more numbers is
taken as the length list, otherwise the first numeric column (a ``size`` /
``bp`` / ``length`` column when there is a header).

Headless: no Tk, so the CLI and the tests use it too.
"""
from __future__ import annotations

import csv
import io
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import fragment_sizing as fs
from fragment_sizing import Ladder

MIN_BANDS = 4                     # fewer cannot carry a sizing curve
MAX_CHANNELS = 8                  # detector channels a ladder may sit on

# display colours for common ladder dyes (matplotlib tab10-ish)
DYE_COLORS = {
    "rox": "#d62728", "et-rox": "#d62728", "liz": "#ff7f0e", "orange": "#ff7f0e",
    "fam": "#1f77b4", "6-fam": "#1f77b4", "hex": "#2ca02c", "vic": "#2ca02c",
    "joe": "#2ca02c", "tamra": "#7f7f7f", "ned": "#bcbd22", "tet": "#17becf",
    "cy5": "#9467bd", "d1": "#9467bd", "wellred d1": "#9467bd",
}


def ladder_dir() -> Path:
    """Where user ladders live (inside the mark library, so tests that point
    the library at a temporary directory move this with it)."""
    import mark_library
    return Path(mark_library.LIBRARY_DIR) / "ladders"


def key_of(name: str) -> str:
    """A file/lookup key from a ladder name: lower case, [a-z0-9_]."""
    key = re.sub(r"[^a-z0-9]+", "_", str(name).strip().lower()).strip("_")
    return key or "ladder"


def default_color(dye: str) -> str:
    return DYE_COLORS.get(str(dye or "").strip().lower(), "#555555")


def to_dict(ladder: Ladder) -> dict:
    def num(x):
        return int(x) if float(x).is_integer() else float(x)
    return {
        "name": ladder.name,
        "vendor": ladder.vendor,
        "dye": ladder.dye,
        "color": ladder.color or default_color(ladder.dye),
        "lengths": [num(b) for b in ladder.lengths],
        "exclude": [num(b) for b in ladder.exclude],
        "default_channel": ladder.default_channel,
        "source": ladder.source,
        "notes": ladder.notes,
    }


def validate(ladder: Ladder) -> List[str]:
    """Problems that make *ladder* unusable (empty list: fine)."""
    problems = []
    ls = list(ladder.lengths)
    if not str(ladder.name).strip():
        problems.append("the ladder needs a name")
    if len(ls) < MIN_BANDS:
        problems.append(f"a ladder needs at least {MIN_BANDS} bands "
                        f"(got {len(ls)})")
    if any(b <= 0 for b in ls):
        problems.append("band lengths must be positive")
    dups = sorted({b for b in ls if ls.count(b) > 1})
    if dups:
        problems.append("duplicate band lengths: "
                        + ", ".join(f"{b:g}" for b in dups))
    stray = [b for b in ladder.exclude if b not in ls]
    if stray:
        problems.append("excluded bands not in the ladder: "
                        + ", ".join(f"{b:g}" for b in stray))
    if ladder.default_channel is not None and not 1 <= ladder.default_channel <= MAX_CHANNELS:
        problems.append(f"default channel must be 1-{MAX_CHANNELS}")
    usable = len(ls) - len([b for b in ladder.exclude if b in ls])
    if len(ls) >= MIN_BANDS and usable < MIN_BANDS:
        problems.append(f"after exclusions fewer than {MIN_BANDS} bands remain")
    return problems


# --------------------------------------------------------------------------- #
# the library
# --------------------------------------------------------------------------- #
def builtin_ladders() -> Dict[str, Ladder]:
    return {k: fs.ladder_from_dict(v) for k, v in sorted(fs.BUILTIN_LADDERS.items())}


def user_ladders() -> Dict[str, Ladder]:
    """The user's ladders, by key.  A file that cannot be read is skipped
    rather than breaking every sizing run."""
    out: Dict[str, Ladder] = {}
    d = ladder_dir()
    if not d.is_dir():
        return out
    for path in sorted(d.glob("*.json")):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
            out[path.stem] = fs.ladder_from_dict(data, default_name=path.stem)
        except (OSError, ValueError, TypeError):
            continue
    return out


def all_ladders() -> List[dict]:
    """Every ladder for a picker: ``{"key", "ladder", "builtin"}``, bundled
    first."""
    rows = [{"key": k, "ladder": lad, "builtin": True}
            for k, lad in builtin_ladders().items()]
    rows += [{"key": k, "ladder": lad, "builtin": False}
             for k, lad in user_ladders().items()]
    return rows


def find(name_or_key: str) -> Optional[Ladder]:
    """A ladder by key or by (case-insensitive) name, bundled or user."""
    text = str(name_or_key).strip()
    key = key_of(text)
    for row in all_ladders():
        if row["key"] == text.lower() or row["key"] == key \
                or row["ladder"].name.strip().lower() == text.lower():
            return row["ladder"]
    return None


def save_user_ladder(ladder: Ladder, key: Optional[str] = None,
                     overwrite: bool = False) -> str:
    """Write *ladder* to the user library and return its key.

    Refuses an unusable ladder, a key that names a bundled ladder, and an
    existing key unless *overwrite* (an edit) is set."""
    problems = validate(ladder)
    if problems:
        raise ValueError("; ".join(problems))
    key = key or key_of(ladder.name)
    if key in fs.BUILTIN_LADDERS:
        raise ValueError(f"{ladder.name!r} would replace the bundled ladder "
                         f"{key!r}; give it another name")
    d = ladder_dir()
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{key}.json"
    if path.exists() and not overwrite:
        raise FileExistsError(f"a ladder called {ladder.name!r} is already in "
                              "the library")
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(to_dict(ladder), indent=2) + "\n", encoding="utf-8")
    tmp.replace(path)
    return key


def delete_user_ladder(key: str) -> None:
    if key in fs.BUILTIN_LADDERS:
        raise ValueError("bundled ladders cannot be deleted")
    path = ladder_dir() / f"{key}.json"
    if path.exists():
        path.unlink()


def unique_name(name: str) -> str:
    """*name*, or *name* (2), (3)... so it does not clash with the library."""
    taken = {key_of(r["ladder"].name) for r in all_ladders()} | \
        {r["key"] for r in all_ladders()}
    if key_of(name) not in taken:
        return name
    n = 2
    while key_of(f"{name} ({n})") in taken:
        n += 1
    return f"{name} ({n})"


def export_ladder(ladder: Ladder, path) -> None:
    Path(path).write_text(json.dumps(to_dict(ladder), indent=2) + "\n",
                          encoding="utf-8")


# --------------------------------------------------------------------------- #
# import
# --------------------------------------------------------------------------- #
_NUM_DOT = re.compile(r"[-+]?\d+(?:\.\d+)?")
_NUM_COMMA = re.compile(r"[-+]?\d+(?:,\d+)?")
_decimal_comma = False            # set per file by parse_lengths_text


def _numbers(text: str) -> List[float]:
    if _decimal_comma:
        return [float(t.replace(",", ".")) for t in _NUM_COMMA.findall(text)]
    return [float(t) for t in _NUM_DOT.findall(text)]


def parse_lengths_text(text: str) -> List[float]:
    """Band lengths from free text, a CSV/TSV or a Genetic Profiler table.

    A row holding four or more numbers is a horizontal list (the longest such
    row wins -- a GP ``test_plate.tdb`` keeps its ladder on one row).
    Otherwise the lengths run down a column: the ``size``/``bp``/``length``
    column when a header names one, else the first numeric field per row."""
    global _decimal_comma
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return []
    # A comma is a decimal mark only in a European-style file: fields split
    # by ";" or tabs and no "." anywhere.  Otherwise it separates values.
    _decimal_comma = (("\t" in text or ";" in text) and "." not in text
                      and "," in text)
    # one comma/semicolon/space separated row of numbers
    rows = [_numbers(ln) for ln in lines]
    wide = [r for ln, r in zip(lines, rows) if len(r) >= MIN_BANDS
            and not re.search(r"[A-Za-z]", ln)]
    if wide:
        return max(wide, key=len)
    # a column: find a header naming it
    dialect_delim = "\t" if any("\t" in ln for ln in lines) else (
        ";" if any(";" in ln for ln in lines) else ",")
    table = list(csv.reader(io.StringIO("\n".join(lines)), delimiter=dialect_delim))
    col = None
    head = [c.strip().lower() for c in table[0]] if table else []
    for i, name in enumerate(head):
        if re.fullmatch(r"(size|sizes|bp|length|lengths|fragment|fragment size"
                        r"|size \(bp\)|length \(bp\)|size_bp|length_bp)", name):
            col = i
            break
    out = []
    for r in table[1:] if col is not None else table:
        cells = [c.strip() for c in r]
        if col is not None:
            if col < len(cells) and _numbers(cells[col]):
                out.append(_numbers(cells[col])[0])
            continue
        for c in cells:
            nums = _numbers(c)
            if nums and re.fullmatch(r"[-+]?\d+(?:[.,]\d+)?", c) and len(nums) == 1:
                out.append(nums[0])
                break
    return out


def import_ladder_file(path) -> List[Ladder]:
    """Ladders read from *path* (JSON, or a text/CSV/TSV/.tdb length list).

    A text file yields one ladder named after the file, dye "custom": the
    caller (the GUI's editor) lets the user name it and pick the dye/colour
    before it is saved."""
    path = Path(path)
    raw = path.read_bytes()
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = raw.decode("latin-1")
    stripped = text.lstrip()
    if path.suffix.lower() == ".json" or stripped.startswith(("{", "[")):
        data = json.loads(text)
        if isinstance(data, dict) and isinstance(data.get("ladders"), list):
            items = data["ladders"]
        elif isinstance(data, list) and data and isinstance(data[0], dict):
            items = data
        else:
            items = [data]
        return [fs.ladder_from_dict(d, default_name=path.stem,
                                    default_source=f"imported from {path.name}")
                for d in items]
    lengths = parse_lengths_text(text)
    if not lengths:
        raise ValueError(f"no fragment lengths found in {path.name}")
    return [Ladder(name=path.stem, lengths=tuple(lengths), dye="custom",
                   source=f"imported from {path.name}")]
