# CTCE handoff (2026-11-18)

State of the CTCE species work (`ctce.py`, CTCE window in `sequence_analyzer.py`),
for continuing on another machine. Branch: `ctce-baseline-hets`, not merged to main.

## What changed on this branch
- **Baseline.** A broad hump (500–2000 counts) sits under the species, mostly Hom2's long tail.
  - `hump_baseline()` removes it (a morphological opening 101 scans wide) before peaks are found and ranked (`species_trace()`).
  - Blocks are measured on `quant_trace()` (flat floor plus the clipped-peak rebuild), above a **straight line between the block ends**. This is what the GUI's area drag does.
- **Blocks** (`species_blocks()`):
  - A block starts at the foot of the main peak's front (`_front`).
  - It ends after its *first* stutter peak, at the valley or where the trace goes flat (`_back`).
  - Neighbours that truly overlap share one line, split by a drop at the valley.
  - The stutter offset is measured per well (`stutter_of`): +24 scans on 040723, about +15 on 010823.
- **Hom3 position ratio** (r3) = median of the dominant peak after Hom2's cluster in the strong wells.
- **Heteroduplexes** are picked as the most *equal-area* pair. They must lie 1.5–4 Hom2→Hom3 distances after Hom2; on #19_20 they sit at about 2.0 and 2.3. A peak just after a comparable one is that peak's stutter, never a het. One het is used only when no pair exists.
- **GUI:**
  - "Clear shown wells" button in the picking/area bar.
  - After the CTCE analysis, a shown well draws the hump (blue band, "not counted") and each block's baseline (red line, block shaded). This works in the normal view and in area-drag mode.

## Reference measurements (user's drags, 040723_mt_nucl_#19_20_dil_curveRun01, A01, Ch2)
| species | drag | area | engine block | engine area |
|---|---|---|---|---|
| Hom2 | 2520–2565 | 644 672 (clipped in the trace) | 2518–2569 | 686 792 (de-clipped) |
| Hom3 | 2665–2698 | 60 419 | 2664–2699 | 61 768 |
| Het4 | 2812–2842 | 27 012 | 2811–2846 | 27 793 |
| Het5 | 2859–2882 | 19 712 (stops before the extra bump at 2885) | 2858–2882 | 18 566 |

The user also confirmed 010823 #19_20 Run01 A01 positions: Hom2 2552, Hom3 2679, Het4 2802, Het5 2831. The engine now hits all four exactly.

## Domain rules from the user
- Analyse one channel, the sample channel. An internal standard on another channel is only used to align runs.
- Peaks of interest change with temperature. Duplexes have a sharp front and a tail; ssDNA is symmetric. Taq gives a stutter peak after each species.
- If a peak is out of scale (flat top), measure it in another channel. This is done by rebuilding the clipped part from crosstalk × fitted gain.
- **Heteroduplexes:** the two hets come from Watson(variant 1)/Crick(variant 2) and Crick(v1)/Watson(v2), so in theory they have equal area. They form at the end of the PCR or with a boil-and-reanneal step after it.
- **Mass action:** with minor-variant fraction p:
  - major hom (1−p)², each het p(1−p), minor hom p²
  - so the major hom is larger than each het, which is larger than the minor hom
  - hets merged into one peak can exceed the major hom
  - at low p the minor hom disappears, leaving 3 peaks (the major hom and 2 hets)
  - minor fraction = (A_minorHom + ΣA_het/2) / A_total

## Open questions and issues
1. **Baseline convention.** The hand-integrated `data(1).csv` (010823 plates) *includes* the hump; the user's new drags exclude it. With the line baseline the engine is about 2.7× below the CSV on Hom3. Proposal: re-baseline the benchmark on fresh drags of 5–10 wells per plate. **Waiting for the user's decision.**
2. **Mass action is violated in the data.** The peak called the second homoduplex is about 30× larger than p² predicts from the het areas, and larger than each het in most wells. Asked the user: was there a final denature/reanneal step? Is that peak really the minor homoduplex? Until this is answered, the ordering is *not* enforced, and 3-peak wells (no minor hom) are not modelled.
3. **The hets are only roughly equal in the data.** Strong wells 040723 A03/B03 have Het4 ≈ 2× Het5. On the strongest half of the wells, Het5/Het4 has a median of 1.02 on 040723 (83 % within ×1.5) and 1.27 on 010823 (60 %; cause unknown).
4. **Faint dilution wells often get only one het** (40 of 92 on 040723): the second is below the 25σ threshold. Idea: a lower threshold for the second het.
5. **010823 A03, A04 and C06 are mislabelled.** For example, A03 picks Hom3 at 1.03 × Hom2 instead of about 1.049; these traces are unusual.
6. **#8 fragment** (010823 #8 Run01/Run02): weaker, with no clear het pattern. More reference drags are needed.

## Data used (home machine)
- mtDNA run zips: `/media/tv/Bilder og Data/tv/Nedlastinger/mtDNA copy-20230831T174221Z-001/`. Plates 010823 / 040723 / 130323, #19_20 and #8. The sample channel is Ch2 (G).
- Hand reference: `/media/tv/Bilder og Data/tv/Nedlastinger/data(1).csv` (semicolon-separated: plate; well; Primer Peak 1; Hom2; Hom3; Het4; Het5).
- The user has more CTCE data at work. New plates with hand drags of a few wells would be the best benchmark.

## How to check
- `python3 -m pytest -q tests/test_ctce.py`: 12 tests on synthetic plates, covering species blocks, the pattern, static exclusion and clipped peaks.
- Benchmark idea:
  1. run `ctce.analyse_plate([run folders])`;
  2. compare `res.areas` and the labelled peaks (`p.label in ctce.SPECIES_NAMES`) to the user's drags;
  3. check the Het5/Het4 equality and the mass-action consistency.
