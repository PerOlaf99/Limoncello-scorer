# Clean-room provenance log: Cimarron 3.12 clone (`cimarron_basecaller/cimarron_clone.py`)

Goal: reproduce Cimarron 3.12's base calls (MegaBACE vendor caller) from raw
traces, learned only from Cimarron's *outputs* (.esd result files) on raw
traces (.rsd) plus public literature. No Cimarron/MegaBACE binary was run,
disassembled or decompiled.

## Sources consulted

| Source | What was used |
|---|---|
| This repository (worktree): `cimarron_basecaller/spacing_caller.py`, `simple_caller.py`, `rsd_io.py`, `configs.py` | Our own independent caller (`track_bases`, GOLDEN preset) and readers; the clone re-uses its signal path. Its "mb.cxx" references point at source printed in the public patent EP0944739A1. |
| `cleanroom/data/m13/{rsd,esd}`, `cleanroom/data/mb4000/{rsd,esd}` | Raw traces and Cimarron 3.12's result files for the same wells (only `--split even` used for development). |
| `cleanroom/data/m13/m13_M77815.1.fa` | M13mp18 reference: used once to check whether junk mb4000 calls (ours and Cimarron's) are M13 sequence at all. |
| `cleanroom/tools/esd_read.py`, `cleanroom/tools/evaluate.py` | .esd tag reader, dye order; scoring tool. |
| `.esd` processed-trace section (decoded from the bytes by us) | 20-byte records from file start to the first tag: int32 + 4 x float32. |
| Patent EP0944739A1, https://patents.google.com/patent/EP0944739A1/en (first 100k chars, via a web fetch summary) | Mobility-shift search (shift vector with one lane fixed, maximise summed envelope); begin/end point heuristics (zones, primer peak, "End set back 350 samples from end"); cubic-spline up-sampling for low samples/band; 2048-scan processing segments with 148-scan overlap. |

Nothing forbidden was opened, listed or run. No `cimarrontv*`, `DLL_FUNCTION_MAP*`,
`CIMARRON_MASTER*` or `decomp*` file was encountered.

Development scripts (not part of the repo) live in the session scratchpad
`work/` folder next to `cleanroom/`.

## Observations (black box) and the decisions resting on them

### .esd layout
- The processed-trace section has exactly `END TRACE PROCESSING - START TRACE PROCESSING`
  records, each `int32 + 4 float32`. The four floats are in **instrument channel
  order** (the run's dye order, e.g. TGCA), not ACGT: re-ordering by `dye_order()`
  makes its peaks line up with Cimarron's called bases (plot of m13 A02 4000-4300).
- The int32 is the raw **Current** column of the .rsd: `marker[i] == current[START TRACE PROCESSING - 1 + i]`
  for 100% of records on every m13 well tested, so processed index *i* is raw
  scan `START TRACE PROCESSING - 1 + i`, and a called base's PEAK POSITION is
  the processed-trace argmax index + `START TRACE PROCESSING + 1` (offset -1 in
  ~95% of bases).

### Coordinate convention (positions)
- Cross-correlating each processed channel with the raw (baseline-subtracted)
  channel: the processed signal is the raw signal delayed by a **per-well**
  amount (G channel, first 2048-scan segment: 5..15 scans on m13), with a fixed
  relative per-dye pattern (A ~ G, C ~ +1, T ~ +5 relative). The delay changes
  in integer steps at roughly the 2048/148 segment boundaries and grows by 2-4
  scans over the read. Consistent with the patent's per-segment mobility
  alignment (one lane fixed), but the absolute per-well part is unexplained.
- Decision: positions = our scan + 7.5 + 0.0003*(scan-2000). The intercept and
  drift are the pooled medians of Cimarron-minus-ours over matched bases (m13
  even: 7 -> 9 along the read; mb4000 normal wells 9 -> 11; one compromise for
  both plates). The per-well spread (+-3) remains: pos<=1 tops out near 45%.
- Wells where Cimarron's positions run at slope 0.5 or 0.25 against raw scans
  (mb4000 broad-peak wells): Cimarron **decimates** the trace by 2 or 4 there and
  reports `raw/k + START SIGNAL`-like coordinates (START SIGNAL ~ START RAW/k).
  Slopes 1.15-1.7 (B22, D02, L24, P22...) match the patent's spline up-sampling.
  Not emulated yet (phase 1).

### Read span
- End: `END SEQUENCE = END TRACE PROCESSING - 50..61` on every well (an edge
  margin). On m13 most capillaries end with a large broad "end-of-run blob" (all
  channels rise ~4x then fall to the pre-run floor); END RAW DATA sits either
  ~35 scans after the blob peak or ~95 before it (two modes, not separable by
  any feature found). Without a blob, END RAW DATA = trace end - 53..120 (m13)
  or - 166..214 (mb4000); no signal feature found at that point (aligned
  averages over 100 wells show none).
  Decision: stop at min(trace end - 240, blob peak - 60) (grid-searched on the
  even wells of both plates: end median 47 / 32 scans, was 146 / 266).
- Start: Cimarron calls the faint early bands from the first sustained rise of
  the signal above the pre-run noise (m13: ~120 scans before track_bases' 10%
  region). Decision: start tracking at the first scan where the 11-scan-smoothed
  max-channel z-score (vs. scans 200-1200) stays > 8 for 10 scans (start median
  error 120 -> 10 scans on m13). On mb4000 Cimarron skips the first ~55 scans
  after START RAW DATA (an isolated leading peak); a start skip/taper was tried and
  did not move agreement/extra, so it is not used.
- Empty capillaries (6 mb4000 wells, flat traces): Cimarron still emits ~100 calls
  at positions beyond the trace end. Decision: return no calls when the 90th
  percentile of the max-channel pre-run z-score is < 12 (empty wells 0.3-7.7,
  all others >= 20). Removes ~21,800 extra bases on mb4000/even.

### Tracker
- Removing track_bases' quality trim and starting earlier exposed spacing
  collapse (B10: the running spacing locked onto half the band spacing; 1,536
  insertions). Clamping the running spacing to [0.7, 1.5] x the read's global
  spacing (GOLDEN uses [0.5, 2.0]) fixed it and cut insertions read-wide:
  m13/even extra 0.143 -> 0.057, agreement 0.928 -> 0.946.
- A second pass anchored to a quadratic expected-spacing curve (patent idea) made
  things worse (the curve's extrapolation at the read end is too low); dropped.
- Broad-peak reads (mb4000): Cimarron's band spacing is ~1.35x track_bases'
  peak-picking spacing estimate (it counts sub-peaks), both on wells Cimarron
  processed at full rate and on wells it decimated; our tracker then called ~2
  bases per band (D02, M10, O02, O14 at ~0 agreement). The autocorrelation of
  the channel-summed processed signal, searched in [0.8, 2.2] x that estimate,
  matches Cimarron's spacing on clean and moderately broad wells (D02 14/14,
  M10 13/13, G10 17/17, H10 18.5/18) but reads ~5-20% long on clean m13 reads
  (9.7 vs 9.2). Decision: reference spacing = max(peak estimate, autocorr/1.2),
  or 1.4 x the peak estimate when no periodicity is found. A hard switch
  (replace when autocorr > 1.25 x estimate) was rejected: the m13 ratios reach
  1.23, too close. mb4000/even 0.792/0.385 -> 0.839/0.170; m13/even
  0.945/0.055 -> 0.945/0.047.
- Tried and dropped: a start skip / linear start taper (no effect), a hard
  "broad read" clamp floor by spacing threshold (works, 0.827/0.200 on mb4000,
  but the threshold sat on the boundary of clean reads).
