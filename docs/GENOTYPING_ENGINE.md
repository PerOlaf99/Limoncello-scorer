# SNP genotyping engine (`snp_genotyping.py`)

Started 2026-11-27. The engine calls SNP genotypes from CTCE runs that have a
four-peak internal standard (STD), with no marks.

## What a well looks like

The standard is the SNP fragment as a heterozygote, run in its own dye. It
shows four bands, always in this order:

    Homo1 → Homo2 → (long gap) → Hetero1 → Hetero2

The two heteroduplexes can run so close together (5–12 scans on ABCC2 and
IL10) that they merge into one band.

The sample is the same fragment in another dye. It runs 6–15 scans ahead of
the standard. The calls are:

- **hom1:** a band at Homo1 only.
- **hom2:** a band at Homo2 only.
- **het:** both homoduplexes, at any ratio (allelic imbalance is flagged
  `ai`). When the heteroduplexes resolve, their band(s) are present too.

## How the engine works (per plate)

1. **Channels.**
   - The sample is the plate's strongest dye.
   - The standard is looked for first in the channel least like the sample,
     by the plate's median correlation with the sample channel. This matters
     because the sample bleeds strongly into Channel 1, and in a het that
     bleed copy has four bands too.
2. **Spacing, by plate vote.**
   - Each well proposes ordered band sets from its 14 strongest
     standard-channel peaks. The spacing most wells reproduce wins.
   - A proposal only counts if:
     - Homo1→Homo2 is 45–200 scans;
     - Homo2→Hetero1 is 120–330 scans and at least 1.3 × the first gap
       (Per's marks: 49–140 and 135–291 over six genes);
     - its weakest band is at least 30% of its strongest (one fragment, one
       amount);
     - no band except Hetero2 is a satellite, i.e. has a peak over twice its
       height within 30 scans before it;
     - the sample has a peak beside Homo1 or Homo2 that is at least half the
       biggest sample peak in the standard's stretch.
   - The four-band pattern wins over the three-band (merged) one when it has
     60% of the wells and either 75% of the vote strength, or it explains the
     three-band pattern as "Hetero1 skipped".
   - A four-band proposal whose fourth band is just the satellite every band
     has is rejected.
3. **Validation.**
   - A candidate standard is kept only if the sample sits beside it, not on
     top of it. The plate's sample offset must be at least 3 scans; an offset
     of about 0 means a dye front or bleed-through.
   - A plate whose PCR failed has no product to anchor on. The standard is
     then looked for without the anchor, and the plate note says so.
4. **Placement in each well.**
   - Spacing varies between capillaries (on T9, Homo1→Homo2 runs 57–100
     scans), so every (H1, H2, Het) set is scored:
     - by its weakest and mean band height;
     - minus a soft penalty for each gap's log-ratio to the plate's;
     - plus bonuses when the sample sits beside a homoduplex, and when the
       partner heteroduplex is where the plate puts it.
   - The partner is searched on both sides, because the band placed as
     Hetero1 may really be Hetero2.
5. **Sample.**
   - The standard's crosstalk is subtracted from the sample channel. The
     bleed fraction is the smallest sample/standard ratio at the standard's
     band scans.
   - The plate's sample offset is learned.
   - Peak heights are read at Homo1, Homo2 and the heteroduplex band(s), in
     units of the noise.
6. **Call** (`call_well`):

   | Condition | Call |
   |---|---|
   | stronger homoduplex below 40× noise | no-call ("no product") |
   | minor/major homoduplex ≥ 0.40, **or** ≥ 0.05 with the heteroduplex band ≥ 0.20 × the major | het (`ai` flag if minor/major < 0.33) |
   | otherwise | hom1 or hom2, by the stronger homoduplex |

   Each row also has `frac2`, the Homo2 share of the homoduplexes. It runs
   from 0 to 1 and is directional, as the LCM work needs.

## Alignment on the standard

`align_plate(arrays, rows, method)` resamples every well onto the plate's
reference axis, using the standard bands the engine placed in each well. The
reference is the plate's median position of each band. `align_map` gives a
single well's mapping. The methods are:

- **piecewise** (default): exact through every band, straight lines between
  them, and the affine slope beyond them.
- **affine**: a shift plus a stretch, fitted to the bands.
- **scale**: Per's 2005 CEQ 8000 normalisation (`KINGSTON/Report.pdf`,
  "Calculations"): each time divided by the sum of the four standard peak
  times, so one factor per well.

The check uses things the alignment doesn't fit: how tightly the sample's
main peak lines up across the hom1 and hom2 wells, and how well three bands
predict the fourth. Means over the 12 scored plates:

| Method | Sample peak spread (MAD, scans) | 4th band from 3 (median error, scans) |
|---|---|---|
| none | 84.9 | – |
| shift (Homo1 only) | 5.4 | – |
| scale (2005) | 9.5 | 10.4 |
| affine | 5.8 | 7.0 |
| **piecewise** | **2.3** | **6.6** |

A single scale assumes migration time starts at zero. On the MegaBACE plates
the capillary-to-capillary drift is mostly an offset plus a stretch, so even
a plain shift beats it. Piecewise wins.

`tools/genotyping_qc/overlay_plate.py RUN OUT.png` draws a plate before and
after alignment, coloured by call. On T9 the 96 wells go from an 800-scan
spread to one stack per species.

## Answer keys

Two files in `limoncello-data/answer_keys/`:

- **`snp_genotypes.csv`** (Per's hand calls, 12 plates, 1,138 wells).
  - Built by `tools/genotyping_qc/build_snp_key.py` from three sources:
    - `~/Dokumenter/Genotype_*.ods` (ABCC2 N11, CYBA N2 at 59 and 60 °C, rs1695 N1);
    - the old GUI's `Genotyping.xlsx` (rs1695 T1 and T5, IL10 N2, ABCC2 N1,
      CYBA N1 Run01/Run02, GSTA1 N1; code 4 = het);
    - the T9 plate key.
  - The run names were matched to folders. "Run1" is Run01; ABCC2 N11 is the
    010910 run.
- **`snp_std_marks.json`** (Per's standard marks, 17 runs, 1,111 wells).
  - From `~/Dokumenter/IS.csv` and the mark library on
    `78B0C7DE1FA7081C/electropherogram/library/`. The library marks name the
    species.
  - Three marks mean Per saw the heteroduplexes as merged.
  - The `result` calls stored in that library were made by the engine, not
    by Per, so they are not used as a key.

## Results

Run with `python tools/genotyping_qc/benchmark_snp.py --placement`. The raw
runs are read from `/run/media/per/78B0C7DE1FA7081C/OY`.

| Engine | Calls agree with Per | Standard placed as marked |
|---|---|---|
| `genotyping.auto_genotype`, no marks | 238/1,138 = 21% | – |
| old plate model seeded by a plate vote | 644/1,138 = 57% | – |
| **`snp_genotyping`** | **1,078/1,138 = 94.7%** | **963/1,111 = 86.7%** |

The call thresholds were chosen with each plate left out in turn, and came out
the same for 10 of the 12 plates. That leave-one-plate-out score is 94.3%.

Calls per plate:

| Plate | Agree |
|---|---|
| rs1695 T5 | 96/96 |
| IL10 N2 | 95/96 |
| ABCC2 N1 | 94/96 |
| CYBA N1 Run01 | 87/89 |
| CYBA N1 Run02 | 88/91 |
| rs1695 T9 | 92/96 |
| GSTA1 N1 | 90/94 |
| rs1695 N1 | 91/96 |
| ABCC2 N11 | 88/96 |
| CYBA N2 59 °C | 87/96 |
| rs1695 T1 | 87/96 |
| CYBA N2 60 °C | 83/96 |

Most of the disagreements left are of four kinds:

- **Very low AI hets** (rs1695 N1 D03 and E06, T9 G09 and G10): the minor
  homoduplex is 1–3% of the major, which is at the noise.
- **Wells Per left uncalled** because the standard was missing or the run
  bad ("no STD", "0"). The engine doesn't check the standard's own strength
  yet.
- **rs1695 T1:** a weak plate, smeared at 56 °C.
- **CYBA homozygotes with a large stray band at the heteroduplex position.**

## Open

- **Failed-PCR plate, fixed (2026-11-27).** ABCC2 N2 160910 had 16 false
  calls: in those wells the standard was placed on the primer front, where
  the sample's own primer peaks sit. Now any well whose Homo1 lands more
  than 450 scans from the plate's median Homo1 is placed again within that
  window (`POSITION_WINDOW`; Per's marks show real standards up to ~490
  scans off, the fronts ~800). The plate now has 1 call left (E09, a large
  sample peak right at Homo2, possibly a well that worked), and the standard
  is placed as marked in 84 of 93 wells (was 75). A stricter per-well test
  (product beside the standard, not on it) was tried and dropped: it also
  no-called 47 of rs1695 T1's hand-called wells, where the standard is weak
  at 56 °C.
- **No call when the standard is weak or absent in a well.** Per marks such
  wells "no STD".
- **Areas next to heights.** Per's rule is that a correctly drawn area is
  the better measure.
- **Dye-overlap matrix.** `genotype_plate(..., matrix="peaks")` removes the
  plate's own dye overlap first (`spectral.matrix_from_peaks`).
  - **Channel choice (2026-11-27).** Every channel but the sample's is voted
    on, and the strongest vote that passes the checks wins. The checks are
    that enough wells agree and that the sample sits beside the standard
    (offset ≥ 3 scans, so not a bleed copy or a front). The old rule, "the
    channel least like the sample first", picked an empty channel on the
    CYBA plates once the bleed was removed.
  - **On raw traces** (the default) this gives 93.8% (was 94.3%). On CYBA it
    now reads the standard's copy in Channel 1, where the standard's dye
    reads at about 43%.
  - **On corrected traces** it gives 92.7% (was 79–81%), and it finds the
    standard on its own channel, Ch3, on every plate.
  - **Still worse with the matrix:** rs1695 T1 (71/96 against 83/96), and
    standard placement on several ABCC2 plates (72% overall against 86%).
    So raw traces stay the default.
- **GUI:** an "Automatic" mode in Auto genotyping that uses this engine.
- **More keys:** the old GUI's `genotypes_OY_*.csv` cover about 10 genes but
  are not usable as a key (same positions in every well). More of Per's hand
  sheets would widen the test.
