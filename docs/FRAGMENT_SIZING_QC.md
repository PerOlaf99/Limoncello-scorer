# Fragment sizing QC — automatic ladder detection (2026-11-25)

This file records how the automatic ladder detection in `fragment_sizing.py`
(commit 337feda; band correction added later) was checked against the MegaBACE archive. The sizing method
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

## ET550-R band lengths against the clean wells

**Method.** The check used 4,099 clean ET550-R cell ID wells from 82 runs. A
well counts if it shows at least 20 of the 22 bands. For each band, the engine
predicted its length from a local cubic in scan, fitted through the three bands
on each side. It then averaged the difference between that prediction and the
listed length.

A single band that sits off also pulls its neighbours' residuals the other
way. To separate that, the offsets were solved together, with a ridge penalty.
Smooth trends cannot be seen this way, because the per-well curve absorbs
them. Only a band that sits out of line with its neighbours shows up.

**250–550 bp: the listed lengths hold.** Every band is within ±0.26 bp in
every period, and the spread between runs is tight. The largest offsets are a
regular +0.18 at 450 and −0.26 at 475. Both are well below the 0.6 bp typical
sizing error, so no change is needed.

**100–175 bp: two ladder behaviours, split by date.**

| Runs | 125 bp | 150 bp | 175 bp |
|---|---|---|---|
| 2011 to June 2012 (28 runs) | 0.0 | −0.1 | −0.3 |
| September 2012 to 2013 (50 runs) | **+1.4** (raw +2.3) | −0.9 | −0.7 |

In the later runs the 125 bp band migrates about 2 bp long. In 130213, for
example, 100→125→150→175 sit about 180, 160 and 183 scans apart; in 2011
runs the gaps are an even 134, 140 and 143. The 125 band is also lower
relative to its neighbours: 0.87 against 1.18 of the well's median band.

The peaks are clean and single in both periods. Overlays are in
`band125.png` in the QC folder. Run settings do not explain the split:
almost all runs in both periods used 10 kV, 75 min and 60 °C. That points to
a change of ladder lot, or of polymer, in summer 2012.

The split is not perfectly clean:

- 220513 and 22032013 are 2013 runs but behave like the early runs;
- 090511 Run02 and 050311 are early runs with 125 running −2.5 bp;
- 280411 sits at +0.9.

**Looking for the cause (instrument, array, ladder lot).** The per-run
offsets are in `band125_by_run.csv` in the QC folder, together with each
run's backup folder and the size of its `CDErr.txt` and `TotalErr.txt`
logs. The lab ran two MegaBACEs, MB1 and MB2. Cell ID was run by Geir Frode
Øy ("gfø" in run names).

- **Last unshifted run:** `cell_id_12.06.2012` (−0.9 bp). **First shifted
  run:** `050912_cell_id_gfø` (+1.5 bp). There are no ET550-R runs in between,
  so the change happened between 12 June and 5 September 2012.
- **The instrument cannot be read from the data.** The `.rsd` files name only
  the dye set (GT DyeSet2: ET-ROX, FAM, NED, HEX). The backup folders are
  mixed: runs named "MB2" (`17022012_MB2_Cell_Id`, `22122011mb2gfø`) sit in
  `Data_MB1`. Early runs in `Data_MB1` and in `OLD DATA MB2` both sit near
  0 bp. All later runs (`MB2_ALL_DATA…/ALL RUNs Backup`, `MagaBace2_Paulo`)
  sit near +2 bp.
- **Capillary detection errors, September to November 2012.** From 050912
  (about 25,000 lines in `CDErr.txt`) to 211112, almost every run has
  10,000–460,000 capillary detection errors. Before that, runs had 2 to a few
  thousand. From December 2012 the logs are clean again, but the 125 band
  stays at +2 bp. So the errors do not cause the shift, but they point to a
  physical change on the instrument in late summer 2012: a new capillary
  array, polymer or service work.
- **Run settings** are the same before and after (10 kV, 75 min, 60 °C as
  set). The actual temperature is not logged. A single band moving while its
  neighbours stay put fits a change in how that fragment folds, which
  depends on the polymer, the denaturing conditions and the real
  temperature, or a ladder lot made differently. A plain temperature or
  voltage change would shift every band a little.

To settle it, the lab logbook should be checked for 12 June – 5 September
2012:

1. which MegaBACE (MB1 or MB2) ran cell ID before and after;
2. capillary array changes or service visits;
3. polymer (LPA) type or lot, and buffer;
4. the ET550-R lot.

The best test is the two 2013 runs without the shift, 22032013 (−0.5 bp) and
220513 (−0.8 bp). Both have clean error logs. If they ran on the other
instrument, or with an older ladder or polymer lot, that names the cause.

**Effect.** The sizing curve passes through every anchor. In the later runs it
is therefore bent by about 2 bp around 125 bp. A sample fragment near 125 bp
is sized up to about 2 bp short of where the neighbouring bands put it; at
110–140 bp the error is smaller. Sizes outside 100–175 bp are not affected.

**Fix: per-run band correction (`run_band_corrections`, on by default).**
The fix needs a run with at least 8 clean wells. In such a run, each band is
measured against its neighbours in every clean well, with at least two
neighbours on each side. A band is corrected when two things hold:

- its median offset is at least 0.75 bp and at most 4 bp;
- at least 75% of the wells agree on the direction.

Bands are corrected one at a time, the worst first, up to 3 per run. That
order matters, because one band off pulls its neighbours' offsets the other
way. The correction is kept only if it lowers the run's typical
(leave-one-out) sizing error.

The corrected lengths are used for the curve, the sizes and the fit quality.
The band table keeps `band_bp` (listed) and adds `band_used_bp` (used). The
run's note lists the corrections, e.g. `125 -> 127.17 bp`.

| Run | Bands corrected | Typical sizing error |
|---|---|---|
| 130213 (2013) | 125 → 127.17 | 0.64 → 0.28 bp |
| 151112 Run02 (2012) | 125 → 128.54 | — |
| 050311 (2011) | 150 → 152.97, 175 → 176.20, 125 → 124.09 | 0.91 → 0.33 bp |
| 090511 Run02 (2011) | 150 → 153.95, 175 → 176.46, 100 → 100.99 | 1.21 → 0.42 bp |
| 090611 (2011) | none | 0.25 bp |
| 220513 (2013, early pattern) | none | 0.48 bp |

In 130213, sample sizes between 100 and 180 bp move by up to 2.2 bp. The two
odd 2011 runs turned out to have the 150 bp band about 3–4 bp long, rather
than 125 short.

`size_files(..., band_correct=False)` gives the uncorrected sizes.

## The ladder engine, hardened (engine version 8)

The goal was ladder recognition that is solid for any kit and any user, not
only for this lab's ET550-R plates.

**Synthetic benchmark** (`tools/fragment_sizing_qc/benchmark_detection.py`).
The benchmark covers 6 plates of 12 wells for each of the 19 bundled kits,
plus 6 plates with no ladder. The plates vary:

- 4- and 5-channel instruments, with the ladder in any channel;
- capillary speed per well, curvature, and runs that stop before the kit's
  last bands;
- fading long bands, a dropped band, and 4–15% dye bleed;
- pull-up from the sample, and sparse, allele-pair or busy (45 peaks)
  sample channels.

| | First run (4 plates per kit) | Final (6 plates per kit) |
|---|---|---|
| Kit and channel right | 54 of 76 plates | **114 of 114** |
| No-ladder plates named "none" | 4 of 4 | 6 of 6 |
| "Kit uncertain" on a right answer | 13 | 0 |
| Time per 12-well plate | 45–70 s before the speed-up | 2.5–14 s |

In the first run, every wrong answer was a nested kit (ET550-R for ET900-R,
ILS 500 for ILS 600, SS400 vs MapMarker 400, MapMarker 1000 vs 1500).

What fixed it:

- **Speed.** The alignment step now runs in plain Python and gives identical
  results (checked on 20,000 random cases); it was 7× slower before. Kits
  with too few band-height peaks are skipped. Kits with identical tables are
  fitted once. Once one combination is clean in 2 wells, later wells only
  try the combinations that have fitted.
- **Nested kits.** ET550-R sits inside ET900-R, ILS 500 inside ILS 600,
  SS400 inside SS600, and MapMarker 1000 inside 1500. Tall peaks after the
  last matched band now count as unexplained. Among kits that explain the
  channel equally, the engine leaves out fits more than 3× looser than the
  tightest, then prefers the most bands matched (within one band), then
  the fewest kit bands. A kit with the same lengths over the bands seen is
  not reported as a doubt.
- **Runs that stop early.** A kit is judged on the bands the run reached,
  if that is at least 60% of the kit. The alignment no longer pays for
  bands past the end of the trace. Without this, an evenly spaced kit
  (GeneFlo 1000) slid by six bands, and ILS 600 was bent onto it.
- **Dyes.** A LIZ kit is only tried in channel 5 of a five-dye instrument.
  Between identical tables, the dye decides.
- **Weak plates.** When the leading answer has fewer than 6 clean test
  wells, detection reads the rest of the plate before deciding. This also
  tells the ladder's channel from its bleed on weak plates.
- **One engine.** A run's *ladder profile* is worked out once from all its
  wells: kit, channel, band template and band corrections
  (`find_run_ladder`, `RunLadder`). It is cached in
  `~/.cache/limoncello-scorer/ladders/`. The GUI's sizing now goes through
  the same `size_files` as `scorer.py size`. Its dialog offers "Automatic"
  for the kit and for the channel, and reports each run's profile.
- **All kits for all users.** "Automatic" tries every bundled kit and every
  ladder in the user's own library. The MegaBACE-only list survives as
  `MEGABACE_KITS`.

**The real runs.** This is the 160 runs that sized, compared with the
earlier engine (`tools/fragment_sizing_qc/compare_runs.py`, output in
`v8_wells.csv`):

| Measure | Result |
|---|---|
| Runs with the same kit, channel and clean wells | 147 |
| Runs that changed | 13 |
| Clean cell ID wells | 4,229 → **4,330** (486 of them rescued) |
| Leave-one-out sizing error, median | 0.30 bp |
| Leave-one-out sizing error, 90th percentile | 0.65 bp |
| Wells outside cell ID runs passed with a named ladder | **0** |
| Runs with band corrections | 54 |

The error was 0.61 and 0.95 bp before the band corrections.

ET550-R is found in 80 runs, 3 of them in channel 3. 76 runs report no size
standard.

What changed, run by run:

| Run | What happened |
|---|---|
| 220213 buffer test | ILS 600 with 1 clean well → ET550-R with 59 |
| Cell_id_test_av_ILS600 | ILS 500 → **ILS 600** (9 clean wells), as its name says |
| 050213 | 0 → 43 clean wells |
| 070411, 23012012, 151112 Run01, 220413_CWOT_PCR | → none: at most 2 clean wells on the whole plate |
| 251012 Run01, 280411, 060912_GFØ | Same kit; about 5 fewer clean wells each, because fits are now stricter |
| 21122011 | ILS 600 "kit uncertain", 18 wells; an earlier pass said ET550-R in channel 3 |
| 050313 | ILS 500 vs ILS 600 vs ET550-R; still unresolved |

21122011 and 050313 are weak plates where kits that share lengths fit about
equally. Both say so in their notes.

**Known weakness.** 170611cellid is a short run, about 3,970 scans, so the
ladder reaches only about 375 bp. ET550-R fits some of its wells exactly,
but on others its alignment fails, and ILS 500 fits instead. The run is
named ILS 500 and is most likely ET550-R. Short runs on kits that share
lengths remain the hardest case.

## Open points

- **Lot-level confirmation.** Lab records of when the ET550-R lot changed in
  2012 would confirm the cause.
- **Real data from other kits.** Every non-ET550-R result above is
  synthetic. LIZ, Beckman and MapMarker runs from other labs would test the
  engine on real ladders.
- **170611cellid** (see *Known weakness* above).
