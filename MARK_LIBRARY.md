# The mark library

Status: implemented. Recorded in `mark_library.py`, `genotyping.PlateISModel`,
and the `_library_*` methods of `sequence_analyzer.py`.

The library is the answers file. Every internal-standard mark the operator
places, and every call auto-genotyping reads from it, is written beside the app
as it is made, so a restart never costs a marked plate, a fitted standard can
seed a later plate, and the corrected calls pile up as labelled rows for the
machine-learning export. It is one directory, `library/`, next to the
application (portable — a USB stick keeps its own history), and it is
git-ignored so an operator's plate history never lands in the repository.

Nothing in the library is load-bearing for opening a file. A missing, locked or
unreadable store degrades to "no marks": a well still opens and labels still
render. Losing the disk can cost a re-mark, never a run.

## What is stored, where

```
library/
  prefs.json                 UI preferences (e.g. the minimum het fraction)
  <run>.json                 one file per run: marks + calls
  models/<run>.json          one file per run: the fitted standard
```

A *run* is a data folder (the folder that holds the `.rsd` files); its library
name is that folder's name. One write per run per action, atomic (`.tmp` then
rename), so a crash mid-write cannot half-write a plate.

### `library/<run>.json` — marks and calls

```json
{"run": "OY_rs1695_T9_270910Run01",
 "updated": "2026-10-13T12:00:00",
 "wells": {
   "A01": {"std": [[2135, "HOM1"], [2214, "HOM2"], [2425, "HET1"], [2520, "HET2"]],
            "std_source": null,
            "std_col": 2,
            "result": {"call": "hom-1", "frac": 0.02, "flags": "",
                       "std_scans": "2135,2214,2425,2520",
                       "label_source": "engine"}}
 }}
```

- `std` — the standard mark: `[scan, name]` pairs, exactly what the picker's
  standard holds.
- `std_source` — how the mark came to be. `null` is an operator's hand mark;
  `anchored` / `anchored-merged` / `relaxed-…` say the plate model placed it.
  The distinction is load-bearing: Auto-genotype re-detects a model-placed mark
  through the model (that is what produces the source and confidence columns)
  and must not mistake it for a hand mark.
- `std_col` — the channel (ACGT index) the standard lives on, so a restart
  re-marks the same channel without re-asking.
- `result` — one auto-genotype row (well, call, frac, the four areas and
  significances, flags, reason), plus `label_source`:
  - `engine` — the automatic call;
  - `edited` — the operator corrected it in the results table.

Only wells that carry a standard are stored as marks. Peaks picked but not yet
named a standard stay in the session: the library is a record of answers, not of
half-finished clicks.

### `library/models/<run>.json` — the fitted standard

Written as the standard is learned (Semi automatic, or Auto-genotype on a
marked plate), and read back to seed a plate that has no marks of its own.

```json
{"run": "OY_rs1695_T9_270910Run01", "updated": "…",
 "model": {"geometry": {"d1": 81, "d2": 64, "d3": 10, "tol1": 12, …},
           "merged_geometry": null, "merged_seen": false,
           "centers": [2200, 2250, …], "n_marked": 6,
           "std_snr": [12.0, …], "std_heights": [900.0, …],
           "co_mig": [[-10, -9], [], [3], []], "sample_snr": [40, …],
           "scaled": null}}
```

The subdirectory matters: the scan that lists runs reads `library/*.json`, so a
model in `library/` itself would look like a well-less mark file. Keeping models
under `models/` keeps the two apart. `mark_library.runs()` and `summary()`
therefore never see a model.

## How it is used

**In a session.** Opening a well hands the picker back the mark it was saved
with (`_library_apply`, via `mark_library.apply_picker`), on the saved channel.

**Across a restart.** The same read is what reseeds a plate: nothing special,
the marks simply come from disk, and Auto-genotype fits the plate shape from
them as it would from hand marks. `_auto_learn_from_marks` counts
library-restored marks like hand marks — that is the point of saving them.

**Across plates.** Auto genotyping ▸ **Use a saved standard…** lists every
stored standard and hands one to the current sample window. It is used only for
a run that has **no marks of its own**: a plate the operator has already marked
always fits those first. Forget learned standard drops it again. The
transferable parts are the shape, the tolerances, the merged template and the
sample/standard co-migration; the band *positions* (`centers`) are plate
specific, so seeding works best plate-to-plate within the same assay. Because
`PlateISModel.find` also searches a stage without the position window, a drifted
centre is handled rather than fatal.

### The machine-learning export

**Auto genotyping ▸ Export library for ML…** gathers every run on disk (plus
anything still in the session) and flattens it to one row per genotyped well —
the measured fields plus the call and its `label_source`. Written as CSV by
default, JSON Lines if the name ends in `.jsonl`/`.ndjson`:

- every call the engine makes is one training row;
- every operator correction rides alongside it as `label_source: edited`, so a
  corrected call is never confused with the engine's guess.

The export does not train anything itself; it produces the file a model would be
trained on.

## Reusing a fitted standard

**Auto genotyping ▸ Use a saved standard…** lists the stored models newest
first, each with its one-line `describe_plate()` summary, and hands the chosen
one to the current sample window. It is applied by the next **Auto-genotype**,
and only to a run that has no marks to fit from — the operator's own marks
always win. **Forget learned standard** clears both the current session's models
and the loaded saved one.

The transferable part is the *shape*: `d1`/`d2`/`d3` and their tolerances, the
merged (three-band) template, and the sample/standard co-migration offsets
primed from the marks. The absolute band **positions** (`centers`) are
plate-specific, so when a saved model is applied to a different plate the
centres may be off; `PlateISModel.find()` also runs its stages without the
position window, so a drifted centre is still found and ranked on shape. Seeding
is therefore most reliable plate-to-plate within the same assay, and a plate the
operator has marked always fits its own first.

## Related

- `mark_library.py` — the store (read/write, models, prefs, ML flatten).
- `genotyping.PlateISModel.to_dict` / `from_dict` — the serializable form, plus
  `ISGeometry.to_dict` / `from_dict`.
- `tests/test_model_library.py` — the round-trip and store tests.
