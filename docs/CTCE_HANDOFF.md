# CTCE handoff (2026-11-19, work PC)

State of the CTCE species work (`ctce.py`, CTCE window in `sequence_analyzer.py`).
Everything is on `main`.

## 2026-11-23: answer keys in place, and the plan for a general engine

**Data** (NVMe disk, `limoncello-data/`; the rules are in `answer_keys/LEGEND_DECISIONS.md`):
- **Catalogue rebuilt:** 7,770 runs, now including `from_mybook` (348) and `from_sids` (69), which were skipped before.
  - Script: `tools/answer_keys/merge_catalog.py`.
  - The content fingerprint reads only the first 64 KB, so it can't tell a renamed copy from the original. For example, `120618_#1_P7F9_STD#2_Arstad Run02` (MyBook) is `120618_#2…`'s data under a new name. **The sample is #2.**
- **Rulings on the paper sheets (operator):**
  - red "0.0" = WT;
  - "Dårlig PCR" = no PCR, the block is not scored;
  - "OBS template" = the call is kept;
  - a strike-out = re-evaluation (crossed out with nothing instead = WT).
- **Lavage 2015 + Christian's Excel** (`add_lavage_excel.py`):
  - The paper is the primary record; an empty Excel cell = not entered.
  - Excel code 3 "neg" = no PCR product.
  - The key now carries a sample ID per well.
  - **P1F5:** Excel put plate 2's failed PCR under plate 1, and the runs confirm the paper.
- **Expert quantitative keys: fractions from peak heights picked by hand** (the operator and Paulo, every electropherogram):
  - `testis5_3d_p2f1/p2f2_scoring.csv`, ~3,100 spots (M1 = WT, M2… = mutants, up to 3 per spot);
  - `tumor_3d_scoring.csv`, ~4,900 spots: T170 breast tumour + N170 node + L170 liver met (P1F1), T231 bladder + LN231 node
    (= M231/Met231; P2F6, P5F6), Liver13 (P1F15). Excel "Fractions" = the WT share.
  - Builders: `tools/answer_keys/build_testis_3d.py`, `build_tumor_3d.py`.
  - LN231 P2F6 S41/S46 and P5F6 S46 share ~30 identical values, probably copy-paste; they are flagged.
  - The workbooks space slices evenly; sections are 12 µm apart (`z_by_section`).

**Height vs area (operator):** a correctly drawn area = the molecules that passed the detector, and is more accurate than height.
Plan:
1. The engine measures **both** on the same peaks.
2. Check its *height* fraction against the hand heights. This tests the peak picking.
3. Then measure area vs height, to get a correction per fragment.

**Design for the general engine** (agreed direction; no end-to-end ML):
1. **Signal:** baseline, sample channel, alignment on the STD (or on each well's own main peak; without a STD, peaks drift by up to 300 scans between capillaries).
2. **Product present?** Main-peak height and cleanliness against the plate.
3. **WT shape from the plate itself:** the median of its most ordinary wells. This works when every sample carries a *different* mutation (mutation scans), which a fixed list of peak positions cannot handle. Fall back to a library WT shape when most wells are mutant (LCM tumour plates).
4. **Extra peaks beyond the WT shape:**
   - 2–3 sharp peaks after the main peak = mutant;
   - heteroduplexes only = low MF;
   - stutter and the +A shoulder must not count.
5. **Fit and quantify:** areas and heights → p.
6. **Call:** WT / LOW / MUT / no PCR, from a small classifier on those numbers, trained on the paper keys. Validate leaving out a whole study at a time.
7. **GUI:** a review queue sorted by confidence; corrections feed the library and the training set.
8. **First:** a harness that scores the engine against every answer key.

**Prototype** `ctce_consensus_prototype.py` (steps 2–4 only; Lavage plates, against the paper key):
- **No PCR vs WT** by main-peak height: AUC 0.83–0.99 on 8 of 9 plates.
- **MUT/LOW vs WT:** AUC 0.57–0.96 (P1F13 p1 0.96, P2F2 p2 0.87, P2F6 p2 0.80).
- **Why it misses:** messy no-PCR wells crowd out the mutants, so the product test must come first. Stutter-aware extra-peak scoring is the main thing left.
- **P1F13 plate 2 is inverted** (the paper's no-PCR wells have *taller* product than WT, AUC 0.008), so a plate mix-up is suspected. Check it in the GUI.
- **The 090215 P2F1/P2F2 Lavage runs** may not be the runs that were scored: the sheets are dated 19.03.15.

**Latest (work PC, 2026-11-19 evening): towards one general engine.** The operator's requirements:
- The engine must work for any fragment.
- It reports a **directional variant-2 fraction from 0 to 1**. LCM wells (~30 cells) span the whole heteroplasmy range, and SNP plates are hom1/het/hom2.
- A lone homoduplex gives 0 or 1.
- Runs with the internal standard (STD) need no marks.

Findings:
- **Today's species search is #19_20-specific.** On rs1695 T9 it found 48 of the 306 species the standard marks. It assumes:
  - Homo1 is the largest peak;
  - the heteroduplexes lie 1.5–4 Homo1→Homo2 distances after Homo1;
  - a band around the primer peak.
- **The block measurement carries over.** Placed by the standard, it gives MF within 0.04 (mean) of `tests/data/rs1695_measured.csv`.
- **A general fit works.** It learns the plate pattern as spacing ratios between the four species, fits it per well on prominent clusters, snaps to real candidates and checks chemistry (the larger homoduplex is never much smaller than a heteroduplex).
  - On T9 it gives MF mean |diff| 0.012, with 96 % of wells within 0.05.
  - On T5 the plate-learned pattern matched T9's, and MF came within 0.05 in 81 % of wells.
  - T2 did badly in that test, but its pattern (≈1.088 / 1.108) was probably right: the operator's marks at 55 °C give heteroduplexes at 1.076 / 1.086, against 1.134 / 1.174 on T9 at 54 °C, so the spacing moves a lot with temperature. The failure came from wrong simulated "marks", which put Hetero1 and Hetero2 at the same place.
  - **Agreed design:** the operator clicks Homo1 → Homo2 → Hetero1 → Hetero2 on ~6 wells per fragment and temperature (heterozygous wells on SNP plates). Learn the pattern from those, fit every well.
  - Scripts: `NVMe Disk/limoncello-data/ctce/scripts_20261119/` (`gfit.py`, `rs_IS.py`, `rs_fitscore.py`).
- **Standard (STD) detection, scored against the operator's GUI clicks** (197 wells, 8 plates): 155 exact (≤ 3 scans), against 38 for the old single-channel search.
  - **Rule:** exclude the sample's own channel (the strongest; heterozygous samples also make four peaks), search every other channel (no dye matrix is applied, so the standard bleeds and the operator marked whichever channel showed it best), skip clipped copies, take the strongest clean quartet.
  - **Still wrong:**
    - IL10, where the heteroduplex bands are 5–7 scans apart and the 9-point smoothing in `genotyping._is_candidates` merges them;
    - a weak heteroduplex band, which pulls in a later stray peak (CYBA N1).
  - **Needs the operator's OK:** the fix goes into `genotyping.py`, which the genotyping auto-caller shares.
  - **Fragments:** `oy_standard_survey.csv` ranks them; rs1695, ABCC2, CYBA and GSTM1 have the cleanest standards.
- **Answer key:** only the operator's clicks in Limoncello's Genotyping tool.
  - Here: `~/Dokumenter/IS.csv`, `RS1695_T9_IS_Peaks.csv`, `CYBA_N1_A01_C04_IS_Peaks.csv`.
  - At home: `library/*.json`, including a 96-well rs1695 plate. `library/` is gitignored, so bring it over.
  - **Not usable:** `~/CE/90_genotyping_tools/is_peaks_*.csv`, an older tool that writes the same positions in every well.
- **Two-temperature test:** T1 is not usable (no standard at 56 °C, smeared traces). Use `OY_rs1695gradient_150910` (57 °C / 59 °C).
- **LCM failed PCRs** (primer and primer-dimer only) need a rule. On `080419_MIT97_P1F14_96spots_LCM` the peak positions show no clean gap between the primer region and the product, so ask the operator to mark a few failed and a few amplified wells first.
- **More answer keys:**
  - `excel_peaks/C276_C528/` (Excel: nine named slots per well, with positions, heights and areas); the raw runs aren't on the disk yet.
  - `Chris_96` (279 runs): real tumour mutations scored by eye, no standard.
- **New in the program (uncommitted, tests pass):** `analyzer_core.temperature_cycling` — a cycled run makes the current ripple — shown as "Temperature cycling" in Run info.
- **Data:** see `NVMe Disk/limoncello-data/README_DATA.md` (catalogue, by-kind index, 110 GB of verified duplicates removed).

**Earlier (home, 2026-11-18 night):** open item 1 below was done at home, committed in `eb1ff7a`.
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
  - Heteroduplexes that have begun to merge (valley above 20 % of the smaller peak, `MERGED_VALLEY`, was 50 %) are **one area under one line**, split at the valley only for display (operator, 2026-11-19: each heteroduplex counts half towards the homoduplexes, so the split does not matter).
    - Robust CV: 021123 0.23 → 0.21, 010823 0.25 → 0.23, 040723 0.22 → 0.20.
    - The drags were made as two areas, so the engine now sits a little above them: mean |diff| 0.0116 → 0.0156, mean +0.005.
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
2. **Heteroduplex pairing** (fixed 2026-11-19 at work).
   - **The shape:** each heteroduplex is a doublet, its main peak plus its +A stutter 8–16 scans on. The gap inside a doublet overlaps the gap between the two doublets (11–16 scans), so the old "≤ 12 scans is one doublet" rule failed:
     - 021123 Run01 A05 paired Hetero1 with its own stutter;
     - 021123 Run02 A01 lost Hetero2.
   - **Pairing fix:** peaks are now paired in order (`HET_DOUBLET_SPAN` = 18), which matches every operator mark.
   - **Width fix:** `_fwhm` now stops at a valley, so a merged peak is no longer taken for a broad bump.
   - **021123 robust CV:** 0.21 → **0.15**, 95/96 wells.
   - **010823:** 0.23 → 0.26. Columns 2–3 rose to about 0.137 now that both heteroduplexes are found; the dilute columns 7–9 stay at about 0.08–0.09 (item 6).
   - **040723:** 0.20 → 0.21.
   - **Block ends:** the operator ends heteroduplex drags at the knee where the steep fall meets the slow tail. The engine goes on to the floor.
     - Knee ends were tried (fall < 2–7 % of the peak per scan). They made every plate less even (010823 0.26 → 0.27–0.34), so the full extent stays.
     - Because of this the engine reads above the drags (mean +0.008, mean |diff| 0.0175).
   - **Colder Homo1** (021123 Run01 A01) is no longer too large: the block is 2311–2354 with 318.7 k, against the drag's 2307–2375 with 360.6 k.
   - **"Temperatures disagree"** on 021123 (before the pairing fix) fired on:
     - B01 (Hetero1);
     - C02 and C06 (Homo2 at 56 °C).
     - Ask the operator whether these are the wells meant.
2b. **Baseline correction** (2026-11-19, operator's wish for an objective area):
   - **The tick box:** "Baseline" in the view bar (also View → Baseline correction, with View → Baseline window…).
     - It takes the slow background off every trace shown, using `analyzer_core.background_baseline`, the same morphological opening as the engine's `hump_baseline`, 101 scans by default.
     - Drags measure the corrected trace and are measured again when the box is toggled.
     - The CSV records it in its `baseline_correction` column.
   - **Engine:** measuring on hump-corrected traces made no difference: robust CV within ±0.01 at windows 101/201/301. The block lines already follow the background, so the engine is unchanged.
   - **Next:** the operator re-drags 021123 A01–A05 with the box ticked, then the benchmark is rerun against those drags.
3. **021123 Run01 A03:** the heteroduplexes are ~40 % above the drags (17–20 k vs 12.5 k), which gives p 0.138 vs 0.121. The pairing is right; the difference is the block-end convention (item 2): the engine ends at 2751 against the drag's 2735, and its line runs lower. Decide with the operator which convention is right.
4. **021123 Run02 A03 / A05:** the operator placed 3 marks. Their Run02 A05 Homo1 of 39 631 looks like another peak. Ask before using them.
5. **Use the temperature shift between the two runs to *identify* species** (operator's hint), for example to confirm Homo2 in wells where it is weak. Do not use it for quantifying.
   - **Caveat:** on 010823 both runs say 55 °C in the files (Temperature 55, comment "CTCE (55-52)").
6. **Set the detection threshold** (`SPECIES_HET_MIN` = 25 × noise) on the dilute columns (010823 cols 7–12) once there are drags there. 010823 A05: the operator marked a Hetero1 doublet at 2745–2769 whose peaks (232 and 182) are below the 265 threshold, so the engine finds only Hetero2.
7. **The CTCE window** (`sequence_analyzer.py`):
   - show the new names, "below detection" and "Hetero1+2" in the table and the graph labels;
   - show the primer peak label.
8. **rs1695 runs with the internal standard** (OY_rs1695_T1…T10 in "CTCE data") still need testing. These are 4-dye genotyping plates (ET-ROX, FAM, NED, HEX).
9. **Other plates to try:**
   - 040723 dilution (54/55 °C);
   - 110523 PFU (Pfu adds no A, so it is a stutter-free control);
   - 031123 T5 LCM (real samples).

## Tests
- `tests/test_ctce.py`: 17 tests (+ Homo1 marked static by the other temperature, a bump on Homo2's rise, merging heteroduplexes as one area, doublets paired in order, width stopping at a valley).
  - The synthetic species test now expects the drag convention for overlapping heteroduplexes (the sum reads low).
- The full suite passes.
