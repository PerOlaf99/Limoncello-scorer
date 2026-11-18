# CTCE handoff (2026-11-18, night, from the home PC)

State of the CTCE species work (`ctce.py`, CTCE window in `sequence_analyzer.py`).
Everything is on `main`.

**Latest (home, 2026-11-18 night):** open item 1 below was done at home, committed in `eb1ff7a`.
- **First thing at work:** pull, then run the benchmark with `--plate "021123_mt_nucl_#19_20"`.
- **Check 021123's colder Homo1** (the operator's drag at Run01 A01 is 2307–2375).
- **Check the "temperatures disagree" flag:** see which wells it fires on, and whether those are the ones the operator had in mind.
- **Home PC files:**
  - The scratch scripts (plots, diagnostics) were saved to `/home/tv/limoncello_work/ctce/` before a restart.
  - `/home/tv/limoncello_work/ctce_flat_links/` holds symlinks to the run folders, usable as `--data`.

## Where it stands
The engine matches the hand drags on 16 marked runs: on average the engine
minus drags difference is +0.000, and the mean absolute difference is 0.012.
The marked runs cover:
- **Plates:** 010823 and 021123 #19_20.
- **Wells:** A01, A02, A03, A05 and A06.
- **Runs:** both runs of each plate.

| Plate | Expectation | Result |
|---|---|---|
| 021123 #19_20 (96 replicates; col 12 = 20 µl PCR) | one fraction everywhere | 92/96 wells, median 0.129, robust CV 0.21 (was 0.34) |
| 010823 #19_20 dilution curve (cols 1–12, rows = replicates) | flat across the dilution | 92/96 wells, column medians 0.093–0.125, median 0.109 |

Reproduce with:
```bash
python3 ctce_marks_benchmark.py --data "<folder with the run folders>" \
    --marks "docs/ctce_marks/More_measure_#19_20.csv" \
    --plate "021123_mt_nucl_#19_20" --plate "010823_mt_nucl_#19_20_dil_curve"
```
- `docs/ctce_marks/` holds a snapshot of the drags. The live file is
  `More_measure_#19_20.csv` in the repo root, and it is not committed.
- On this machine the data is in `/run/media/per/NVMe Disk/CTCE data/`, with 25 run folders:
  - #19_20: 010823, 021123, 031123 T5 LCM, 040723, 110523 PFU and 130323;
  - the OY rs1695 runs.
- At home the plates are under `/media/tv/Bilder og Data/tv/Nedlastinger/mtDNA …`
  (see the old notes). 021123 may not be there.

## Naming (from the operator)
- **Species:** primer peak, Homo1, Homo2, Hetero1, Hetero2 (formerly Hom2/Hom3/Het4/Het5).
- **"Hetero1+2":** the two heteroduplexes merged into one block. It counts once.

## Rules now in the engine
- **Block:** a species is its main peak plus its whole Taq +A stutter cluster.
  - The stutter cluster is about 15–22 scans after the main peak, is often a doublet, and can be as tall as the main peak.
  - `_back` takes the *tallest* peak near the stutter offset, not the first bump.
  - It then follows the stutter's tail through its own sub-peaks down to a real valley (below 20 % of the stutter peak).
- **Baseline:** every block has its own straight line, split at the valley, as dragged.
  - Blocks share a line only if the valley stays above half of the smaller peak (`MERGED_VALLEY`).
- **Heteroduplex candidates:**
  - They are grouped into doublets (gap ≤ 12 scans).
  - A run of 3+ peaks over 18+ scans is both heteroduplexes, split in the middle.
  - A cluster counts as stutter only if it lies one stutter offset after another cluster *and* is < 60 % of its height.
  - The two heteroduplexes are the pair with the most equal areas.
- **Presence (`present_species`):** Homo2 and each heteroduplex need ≥ 25 × noise and FWHM ≤ max(4 × Homo1's, 20 scans), so broad bumps are not species.
- **Fraction (`minor_fraction_named`):** p = (Homo2 + ΣHet/2) / total, from whatever is present.
  - **Homo2 absent** (low p, ~p²): p rests on the heteroduplexes.
  - **Few heteroduplexes** (dilute wells; the PCR is still running at the end, so there is little re-annealing): p rests on Homo2 / (Homo1 + Homo2).
  - **A lone narrow heteroduplex:** it stands in for its partner.
  - **A lone wide one** (≥ Homo1's block): it is "Hetero1+2".
  - **Neither Homo2 nor a heteroduplex:** **below detection, no number** (the operator's wish: no forced numbers).
- **Temperature:** measured at the temperature that resolves the species best.
  - Replicate runs there are averaged.
  - Averaging the other temperature in made it worse (robust CV 0.19 → 0.34): the heteroduplexes merge or shift at 56 °C.

## Open items, roughly in order
1. ~~**Homo1 at the colder temperature**~~ (worked on 2026-11-18 at home, on 040723 at 54/55 °C and 130323 at 53/54/55 °C; 021123 was not available there). Peaks move over one another as the temperature changes, which caused four kinds of error:
   - **Homo1 called static.** Between close temperatures Homo1 barely moves, so the runs paired it with itself and marked it "static". The engine then took a small peak before it (040723 A03: MF 0.202 → 0.146).
     - Fix: a static-marked candidate at least `HOMO1_OVER_STATIC` (3×) taller than the best other one is still Homo1.
   - **Homo2 lost to a small bump on its rise.** The footprint runs from 10 scans before a peak to 30 after, so a small bump on the rise of Homo2 borrowed its area (H01: a 234-high bump beat a 7433-high Homo2).
     - Fix: a Homo2 candidate with a peak more than twice as tall in its footprint is skipped. Outliers moved toward the plate value (F12 0.004 → 0.095, G12 0.005 → 0.123); 040723 column 12 went from 0.075 to 0.128.
   - **Wrong heteroduplex stutter offset.** The heteroduplex stutter offset is now taken from Homo2, not Homo1 (Homo1 +24 vs heteroduplexes +13–16 on 040723).
   - **Wrong cluster main peak.** A cluster's main peak is its first *big* peak, not a small bump before it. 040723 A01 is now 2527 / 2671 / 2818 / 2865, as the operator marked.
   - **New check: "temperatures disagree".** Every well is also measured at the plate's other temperatures, where they resolve the same species, each ≥ 1 % of the total. A well whose MF *relative to that temperature's plate median* moves by more than `TEMPS_DISAGREE` (×1.3) is flagged, naming the species whose share moved most. On 040723 4 of 29 checked wells are flagged; one is A01, at its Hetero2 (the extra bump the operator found hard to leave out). On 130323, at 53–54 °C a broad peak sits under Homo2. The MFs are stored in `WellResult.mf_by_temperature`.
   - **Benchmark at home** (010823 drags + plates):
     - drags: mean |diff| 0.0054 → 0.0045;
     - 010823: 95/96 wells, robust CV 0.28 → 0.25;
     - 040723: robust CV 0.25 → 0.22.
   - **Still to check at work:** whether 021123's colder Homo1 (operator's drag 2307–2375, 68 scans wide) is now right, and whether the flag fires on the wells the operator had in mind.
2. **021123 Run01 A03:** the heteroduplexes are ~40 % above the drags (17–20 k vs 12.5 k), which gives p 0.138 vs 0.121.
3. **021123 Run02 A03 / A05:** the operator placed 3 marks. Their Run02 A05 Homo1 of 39 631 looks like another peak. Ask before using them.
4. **Use the temperature shift between the two runs to *identify* species** (operator's hint), for example to confirm Homo2 in wells where it is weak. Do not use it for quantifying.
   - **Caveat:** on 010823 both runs say 55 °C in the files (Temperature 55, comment "CTCE (55-52)").
5. **Set the detection threshold** (`SPECIES_HET_MIN` = 25 × noise) on the dilute columns (010823 cols 7–12) once there are drags there.
6. **The CTCE window** (`sequence_analyzer.py`):
   - show the new names, "below detection" and "Hetero1+2" in the table and the graph labels;
   - show the primer peak label.
7. **rs1695 runs with the internal standard** (OY_rs1695_T1…T10 in "CTCE data") still need testing. These are 4-dye genotyping plates (ET-ROX, FAM, NED, HEX).
8. **Other plates to try:**
   - 040723 dilution (54/55 °C);
   - 110523 PFU (Pfu adds no A, so it is a stutter-free control);
   - 031123 T5 LCM (real samples).

## Tests
- `tests/test_ctce.py`: 14 tests (+ Homo1 marked static by the other temperature, a bump on Homo2's rise).
  - The synthetic species test now expects the drag convention for overlapping heteroduplexes (the sum reads low).
- The full suite passes.
