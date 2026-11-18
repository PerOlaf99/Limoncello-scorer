# CTCE handoff (2026-11-18, evening)

State of the CTCE species work (`ctce.py`, CTCE window in `sequence_analyzer.py`)
for continuing at home. Everything is on `main`.

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
1. **Homo1 at the colder temperature** covers an extra peak, so its area there is slightly too large (operator's note). Check whether the measuring temperature's Homo1 block includes it, and trim it.
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
- `tests/test_ctce.py`: 12 tests.
  - The synthetic species test now expects the drag convention for overlapping heteroduplexes (the sum reads low).
- The full suite passes.
