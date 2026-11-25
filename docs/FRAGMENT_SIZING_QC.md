# Fragment sizing QC — automatic ladder detection (2026-11-25)

This file records how the automatic ladder detection in `fragment_sizing.py`
(commit 337feda) was checked against the MegaBACE archive. The sizing method
itself is described in `README.md` under *Fragment sizing → How it works*.

Raw outputs are kept outside the repository, in
`limoncello-data/fragment_sizing/qc_20261125/`:

| File | Contents |
|---|---|
| `auto_wells.csv`, `auto_bands.csv` | All 167 runs, one row per well / per ladder band |
| `rerun_wells.csv`, `rerun2_wells.csv` (+ `_bands`) | Rerun of the 97 problem runs with the final detection code |
| `problem_runs.txt` | The runs that were rerun |
| `summarize_rerun.py` | Before/after table per run |

To reproduce: `python3 tools/fragment_sizing_qc/validate_auto.py OUT.csv BANDS.csv [runlist.txt]`.

## What the engine does now

- **`--ladder auto` is the default.** For each run folder the engine tries every
  known kit (ET400-R, ET550-R, ET900-R, Promega ILS 500/600, GeneScan 400HD/500
  ROX, MapMarker 1000, Beckman SS400) in every channel. It works through wells
  spread over the plate until one kit/channel combination fits cleanly in 6
  wells. It tests at least 8 wells and at most 32.
- **A kit is named only if it fits at least 3 test wells.** If two kits fit
  almost equally well, the result says *kit uncertain*.
- **A fit must explain the tall peaks.** A ladder fit is accepted only if it
  leaves few unexplained tall peaks: at most 0.25 per band.

  | Trace | Unexplained tall peaks per band |
  |---|---|
  | True ladder | 0–0.09 |
  | Wrong kit | 0.19–0.41 |
  | Sample channel | 1.9–20 |

  This test stops kits being "fitted" to busy sample channels, and it stops
  ROX bleed into the neighbouring channel being taken for the ladder.
- **Plate rescue.** Once a run has at least 3 clean wells, the engine builds a
  template of where each band sits. It then re-fits weak wells to that
  template, allowing a stretch and shift plus a quadratic correction. A rescue
  is refused unless the rescued bands are at least 1.5× taller than the other
  peaks in their span.

## Validation on all 167 runs

| Measure | Result |
|---|---|
| Wells outside cell ID runs falsely passed as laddered | 326 before → 1 now |
| Clean cell ID wells | 4,213 |
| … of which rescued by the plate template | 408 |
| Leave-one-out sizing error, median | 0.61 bp |
| Leave-one-out sizing error, 90th percentile | 0.95 bp |

Runs without a ladder still report "no size standard recognised", and nothing
in them is sized. These are the CWOT, cap coating, control plate, PTEN, qPCR
and P1 ligation runs.

## Findings per run

**Cell ID plates that used to fail and now size**

| Run | Ladder | Clean wells |
|---|---|---|
| 191012_cell_id_rekjøring2 | ET550-R, **channel 3** | 43 |
| 231012_cellid_rekjøring3 | ET550-R ch1 | 62 |
| 250213_cell_id_p224_225_p229-234 | ET550-R ch1 | 51 |
| 251012_cell_id_rekjøring_4 Run01 / Run02 | ET550-R ch1 | 39 / 14 |
| 310311cellID | ET550-R ch1 | 36 |
| 060912_cell_id_NY | ET550-R, **channel 3** | 24 |
| 261112cellIDplate221112 | ET550-R ch1 | 20 |
| 151112_cell_id_rerun1 Run02 | ET550-R ch1 | 12 |

**The ladder is not always in channel 1.** 191012, 060912 and 210912 have
ET550-R in channel 3. Detection handles this by itself.

**080611 is not a false ladder.** ROX bleeds into channel 2 at about 10% of its
height. The real ladder is in channel 1, and detection now picks channel 1.

**050213: use Paulo's copy.** The MB2-backup copy of
`050213_cell_id_p320-328Run01` is truncated: its files are about 5 KB instead of
127 KB, and nothing in it can be sized. Paulo's copy sizes 43 clean wells.

**The kit is uncertain on some runs.**

| Run | ET550-R | Promega ILS 600 | Note |
|---|---|---|---|
| Cell_id_test_av_ILS600 | 2 clean wells | 9 clean wells | Both leave 0.05 unexplained. The run name and the clean-well count both point to ILS 600. |
| 050313_cell_id_genetisk_plate1 | 7 clean wells, 0.09 unexplained | 11 clean wells, 0.18 unexplained | Genuinely ambiguous. Check in the GUI before trusting the sizes. |

**Weak plates — do not trust the kit label.** The rerun found only 1–2 clean
test wells on these runs:

- 151112 Run01
- 170611cellid
- 23012012_MB2_Cell_ID_GFØ
- 220213 buffer test
- 220413_CWOT_PCR

The rerun named a kit for them, sometimes an exotic one (Beckman, ILS 500).
The rule that a kit needs at least 3 test wells was added after this rerun.
With it, these runs report no ladder or *kit uncertain*.

## Open points

- **Per-band calibration.** `auto_bands.csv` holds about 4,000 clean cell ID
  wells. It can show whether any ET550-R band runs consistently off its
  listed length.
- **Detection speed.** Detection takes about 2–5 minutes per 96-well run. The
  GUI (`sequence_analyzer.py`) still uses a fixed ladder and channel. To use
  auto detection there, the result should be cached per run.
