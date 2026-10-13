"""Durable store for the operator's marks and the engine's calls.

The library is a record of *answers*: the four internal-standard scans an
operator placed on a well, and the call auto-genotyping read from them.  It
exists so a restart never costs the work of marking a plate -- the next
Auto-genotype reseeds itself from what is on disk -- and so the corrected
calls pile up as labelled rows for the machine-learning export.

One JSON file per run under :data:`LIBRARY_DIR`::

    {"run": "OY_...Run01", "updated": "...",
     "wells": {"A01": {"std": [[2135, "HOM1"], ...],
                       "std_source": "anchored",
                       "result": {"call": "hom-1", ...,
                                  "label_source": "engine"}}}}

``std`` is a list of ``[scan, name]`` pairs; ``result`` is one auto-genotype
row with a ``label_source`` saying whose call it is ("engine", "edited", ...).
A missing or unreadable library is never fatal: every reader degrades to "no
marks" so a locked disk cannot stop a well from opening.
"""
from __future__ import annotations

import csv
import datetime as _dt
import json
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

#: Where the marks live.  Beside the application, so the whole tool stays
#: portable (a USB stick keeps its library) and a test can point it at a
#: temporary directory without touching the operator's real plate history.
LIBRARY_DIR = Path(__file__).resolve().parent / "library"

_PREFS_FILE = "prefs.json"


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def _safe_name(run: Any) -> str:
    name = str(run).strip() or "run"
    for ch in ("/", "\\", ":", "\0"):
        name = name.replace(ch, "_")
    return name


def _run_path(run: Any) -> Path:
    return LIBRARY_DIR / f"{_safe_name(run)}.json"


def _clean(value: Any) -> Any:
    """Turn *value* into something ``json`` can write.

    Rows carry NumPy scalars and arrays once auto-genotyping has run; those
    are reduced to plain Python, and anything still unknown is stringified
    rather than allowed to abort a save (losing a whole plate's marks to one
    odd field would be far worse than an ugly cell).
    """
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, float):
        if value != value or value in (float("inf"), float("-inf")):
            return None
        return value
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_clean(v) for v in value]
    item = getattr(value, "item", None)
    if callable(item):
        try:
            return _clean(item())
        except Exception:                                    # noqa: BLE001
            pass
    tolist = getattr(value, "tolist", None)
    if callable(tolist):
        try:
            return _clean(tolist())
        except Exception:
            pass
    return str(value)


def _atomic_write(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(_clean(payload), ensure_ascii=False),
                   encoding="utf-8")
    tmp.replace(path)


def _read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _runs_on_disk() -> List[Tuple[str, Dict[str, Any], Path]]:
    if not LIBRARY_DIR.is_dir():
        return []
    out: List[Tuple[str, Dict[str, Any], Path]] = []
    for path in sorted(LIBRARY_DIR.glob("*.json")):
        if path.name == _PREFS_FILE:
            continue
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if not isinstance(data, dict):
            continue
        wells = data.get("wells")
        if not isinstance(wells, dict):
            continue
        run = data.get("run") or path.stem
        out.append((str(run), wells, path))
    return out


def _as_pairs(marks: Any) -> List[Tuple[int, str]]:
    pairs: List[Tuple[int, str]] = []
    for pair in marks or []:
        try:
            scan, name = pair[0], pair[1]
        except (TypeError, IndexError, KeyError):
            scan, name = pair, ""
        try:
            scan = int(round(float(scan)))
        except (TypeError, ValueError):
            continue
        pairs.append((scan, "" if name is None else str(name)))
    return pairs


# --------------------------------------------------------------------------- #
# runs: the marks on disk
# --------------------------------------------------------------------------- #
def load_run(run: Any) -> Dict[str, Dict[str, Any]]:
    """One run's saved entries, or ``{}`` if nothing (readable) is stored."""
    path = _run_path(run)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    if not isinstance(data, dict):
        return {}
    wells = data.get("wells")
    if not isinstance(wells, dict):
        return {}
    return wells


def save_run(run: Any, wells: Dict[str, Dict[str, Any]]) -> None:
    """Write a run's entries, replacing what was there.

    *wells* is expected to be the full set for the run (the callers load the
    current file, update one well and hand the whole mapping back), so a write
    is an overwrite rather than a merge.
    """
    payload = {
        "run": str(run),
        "updated": _dt.datetime.now().isoformat(timespec="seconds"),
        "wells": _clean(dict(wells or {})),
    }
    path = _run_path(run)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, ensure_ascii=False),
                   encoding="utf-8")
    tmp.replace(path)


def entry_from_picker(pk: Any) -> Dict[str, Any]:
    """The storable part of a picker: its standard mark.

    ``std_auto`` rides along as ``std_source`` so a plate-model mark can be
    told from an operator's pick after a restart -- the Auto path re-detects
    the former and must not treat it as a hand mark.
    """
    entry: Dict[str, Any] = {"std": [[s, n] for s, n in _as_pairs(
        getattr(pk, "std", None))]}
    # Remember the channel the standard sits on, so a restart can put the marks
    # back into the picker and not only the scan numbers.  Without it a restored
    # standard had no picked mains, so "add the peak that was missed" had
    # nothing to add to and the well was reported already marked and left alone.
    col = None
    for r in getattr(pk, "records", None) or []:
        try:
            if r.get("kind") == "main" and any(
                    abs(int(r.get("scan", -10 ** 9)) - s) <= 2
                    for s, _ in entry["std"]):
                col = int(r.get("col"))
                break
        except (TypeError, ValueError, KeyError):
            continue
    if col is not None:
        entry["std_col"] = col
    src = getattr(pk, "std_auto", None)
    if src:
        entry["std_source"] = str(src)
    length = getattr(pk, "length_bp", None)
    if length is not None:
        entry["length_bp"] = _clean(length)
    return entry


def apply_picker(pk: Any, entry: Dict[str, Any],
                 default_col: Any = None) -> None:
    """Give a freshly built picker back the standard it was saved with.

    The standard's marks are put back into the picker too, not just its scans:
    a restored standard that owned no picked mains could not be re-marked, so
    adding a peak that detection missed had nothing to add to and the well came
    back "already marked, left alone".  *default_col* is the standard channel to
    use when the entry predates ``std_col``.  ``std_auto`` is left falsy for a
    plain mark so auto-genotyping still routes it through the manual-seed path;
    a plate-model source is restored as-is.
    """
    if not isinstance(entry, dict):
        return
    pairs = _as_pairs(entry.get("std"))
    if not pairs:
        return
    col = entry.get("std_col", default_col)
    add = getattr(pk, "add_std_mark", None)
    if col is not None and callable(add):
        try:
            col = int(col)
        except (TypeError, ValueError):
            col = None
        if col is not None:
            for s, _n in pairs:
                try:
                    add(s, col)
                except Exception:                      # noqa: BLE001
                    pass
    pk.std = pairs
    src = entry.get("std_source")
    pk.std_auto = src if src else None
    if "length_bp" in entry:
        try:
            pk.length_bp = entry.get("length_bp")
        except Exception:                              # noqa: BLE001
            pass


def entry_from_row(row: Dict[str, Any], label_source: str = "engine") -> Dict[str, Any]:
    """The storable part of an auto-genotype row, under ``result``.

    The row keeps every measured field -- peak heights, SNRs, the standard's
    scans -- and gains ``label_source`` so the export can separate the
    engine's guesses from an operator's corrections.
    """
    data = _clean(dict(row)) if isinstance(row, dict) else {}
    if not isinstance(data, dict):
        data = {}
    data["label_source"] = str(label_source)
    return {"result": data}


# --------------------------------------------------------------------------- #
# the library as a whole: the export and the "what is on disk" view
# --------------------------------------------------------------------------- #
def summary() -> List[Tuple[str, int, int, Any]]:
    """``[(run, marked, results, updated), ...]``, newest last."""
    rows = []
    for run, wells, path in _runs_on_disk():
        marked = sum(1 for e in wells.values()
                     if isinstance(e, dict) and e.get("std"))
        results = sum(1 for e in wells.values()
                      if isinstance(e, dict) and e.get("result"))
        try:
            updated: Any = _dt.datetime.fromtimestamp(path.stat().st_mtime)
        except OSError:
            updated = ""
        rows.append((run, marked, results, updated))
    return rows


def runs() -> List[Dict[str, Any]]:
    """Every stored run as ``{"run": ..., "wells": {...}}``."""
    return [{"run": run, "wells": wells} for run, wells, _ in _runs_on_disk()]


# --------------------------------------------------------------------------- #
# models: the fitted internal-standard shapes
# --------------------------------------------------------------------------- #
def _models_dir() -> Path:
    """Where fitted models live: a subdirectory under the library.

    Kept under :data:`LIBRARY_DIR` (rather than a module constant) so a test
    that points the library at a temporary directory redirects the models with
    it, and so the ``*.json`` scan that lists runs never trips over one and
    reads a model as a (well-less) mark file.
    """
    return LIBRARY_DIR / "models"


def _model_path(run: Any) -> Path:
    return _models_dir() / f"{_safe_name(run)}.json"


def save_model(run: Any, model: Dict[str, Any]) -> None:
    """Persist one run's fitted standard under its name (overwrite).

    *model* is a ``PlateISModel.to_dict()`` payload.  Best-effort like every
    other write here: a model that cannot be saved costs a re-fit, never the
    plate, so a failure is logged and swallowed rather than raised.
    """
    data = {"run": str(run),
            "updated": _dt.datetime.now().isoformat(timespec="seconds"),
            "model": _clean(dict(model or {}))}
    try:
        _atomic_write(_model_path(run), data)
    except OSError as exc:
        print(f"library: could not save model {run}: {exc}", file=sys.stderr)


def load_model(run: Any) -> Dict[str, Any]:
    """One run's saved model payload, or ``{}`` when there is none/readable."""
    data = _read_json(_model_path(run))
    if not isinstance(data, dict):
        return {}
    model = data.get("model")
    return model if isinstance(model, dict) else {}


def saved_models() -> List[Tuple[str, Any, Dict[str, Any]]]:
    """``[(run, updated, model), ...]`` for every stored model, oldest first."""
    out: List[Tuple[str, Any, Dict[str, Any]]] = []
    models = _models_dir()
    if not models.is_dir():
        return out
    for path in sorted(models.glob("*.json")):
        data = _read_json(path)
        if not isinstance(data, dict):
            continue
        model = data.get("model")
        if not isinstance(model, dict):
            continue
        run = str(data.get("run") or path.stem)
        try:
            updated: Any = _dt.datetime.fromtimestamp(path.stat().st_mtime)
        except OSError:
            updated = ""
        out.append((run, updated, model))
    out.sort(key=lambda row: str(row[1]))
    return out


def ml_rows(wells: Dict[str, Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Flatten ``{run: {well: entry}}`` into one labelled row per genotyped well."""
    rows: List[Dict[str, Any]] = []
    if not isinstance(wells, dict):
        return rows
    for run in sorted(wells, key=str):
        per_run = wells.get(run)
        if not isinstance(per_run, dict):
            continue
        for well in sorted(per_run, key=str):
            entry = per_run.get(well)
            if not isinstance(entry, dict):
                continue
            result = entry.get("result")
            if not isinstance(result, dict):
                continue
            row = {"run": run, "well": well}
            for k, v in result.items():
                if k not in ("run", "well"):
                    row[k] = v
            rows.append(row)
    return rows


def write_ml(path: str | Path, rows: Iterable[Dict[str, Any]]) -> int:
    """Write labelled rows as CSV (default) or JSON Lines; return the count."""
    path = Path(path)
    rows = [r for r in rows if isinstance(r, dict)]
    if path.suffix.lower() in (".jsonl", ".ndjson"):
        with path.open("w", encoding="utf-8") as f:
            for row in rows:
                f.write(json.dumps(_clean(row), ensure_ascii=False) + "\n")
        return len(rows)

    fields: List[str] = []
    for row in rows:
        for k in row:
            if k not in fields:
                fields.append(k)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            flat = {}
            for k, v in row.items():
                if isinstance(v, (dict, list, tuple)):
                    flat[k] = json.dumps(_clean(v), ensure_ascii=False)
                else:
                    flat[k] = v
            writer.writerow(flat)
    return len(rows)


# --------------------------------------------------------------------------- #
# preferences
# --------------------------------------------------------------------------- #
def load_prefs() -> Dict[str, Any]:
    path = LIBRARY_DIR / _PREFS_FILE
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def save_prefs(prefs: Dict[str, Any]) -> None:
    try:
        LIBRARY_DIR.mkdir(parents=True, exist_ok=True)
    except OSError:
        return
    path = LIBRARY_DIR / _PREFS_FILE
    tmp = path.with_suffix(path.suffix + ".tmp")
    try:
        tmp.write_text(json.dumps(_clean(dict(prefs or {})), ensure_ascii=False),
                       encoding="utf-8")
        tmp.replace(path)
    except OSError:
        pass
