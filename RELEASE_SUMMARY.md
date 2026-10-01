# Summary: One-Click Automated Genotyping for Limoncello

## Objective
- Implement the user-approved one-click batch auto-genotyping flow: select wells/run, press one action, automatically call each well, show a results table with progress/failure reasons, and export.
- Preserve the current inline manual peak-picking and drag-area modes; no channel ticking for the automated path.

## Important Details
- User selected the recommended option: "Build it: select wells, one action, results table". Trace-derived ``cut`` work is deferred.
- On this plate's ``ACTG`` dye order (Ch1=A, Ch2=C, Ch3=T, Ch4=G), the standard is **Ch3 (T)** and the sample is **Ch2 (C)**. This pair is exactly ``auto_mark_std``'s ``channel=3`` default; both constants stay explicit because swapping them does not fail loudly — it silently scores the sample's own peaks as the standard and returns confident nonsense. Every result row records the channels it was measured on.
- The ground truth ``rs1695_measured.csv`` is reproduced to 95/96 wells; the only divergence is the documented intentional one for G09 (heteroduplex at the detection floor, called hom‑2 on purpose).
- ``PeakPicker.auto_mark_std(channel=3, cut=1900)`` sets the standard quartet; after marking, four areas and their significances are measured at the standard's peak positions on the sample channel and passed to ``scorer.t9_call(hom1,hom2,het1,het2,snrs)``, which returns ``(call, fraction, flags)``.
- ``find_is_quartet(trace, cut=1900)`` on the raw ``doc.acgt`` column reliably finds the standard quartet in 96/96 wells. The 125‑scan position spread and ~2.6× area ratio make the T/G quartet unambiguous once the dye order is fixed.
- ``noise_sigma`` = 1.4826 × MAD of a 9‑scan Savitzky‑Golay residual; ``t9_call`` uses sigmas = (hom1,hom2,het1,het2) above this noise when given; omitted sigmas let the call rest on areas alone.
- ``SATELLITE_AMBIGUOUS_FRAC = 0.15`` is uncalibrated against known minor alleles; flags it raises are advisory (“look at this”), not confident calls. G09 remains a borderline H2.
- Do not modify ``/home/tv/electropherogram`` without a separate request.
- Preserve rollback refs: ``origin/260929-feat-scorer-packaging`` at ``7098cbe``; unmerged ``feat-abd-qval2-limoncello-theme`` contains ``2fa99d8d`` and ``d06ca51f``.
- Rebuild and clean‑extract ``Limoncello scorer.zip`` (51 files, 128 tests passing), then commit/push without force‑pushing.
- GUI launched with T9 trial folder, PID 42339, log clean.

## Work State
### Completed
- T9 engine integrated (`c5eb9ef`); `find_is_quartet()`, `acgt_index_for_channel()`, `PeakPicker.auto_mark_std()`, satellite ambiguity handling, and scorer updates.
- Inline drag-area mode completed (`4534266`); removed `ManualAreaDialog`/`ManualAreaEditor`, fixed motion guard, made area preview work.
- Mutually exclusive mode checkbuttons and manual‑help coverage completed (`0d9bcca`); archive follow‑up `ff14941`.
- Channel bar corrected (`a7deb3f`); labels now read Ch1→Ch4 in physical order; Run info distinguishes “Dye order (channels)” from “Column layout.”
- Redundant `Signal: Volts` toolbar label removed (`ee8536a`); existing figure‑level `Volt` Y‑label remains.
- Drag channel‑selection help and all‑channels‑off status message completed (`6a8b6fc`).
- Latest validation on the merged line: **168 tests pass, 7 skipped**, compile and whitespace checks pass, and the release zip clean‑extracts with the full suite passing (see *Mass-action MF* below for what was added after the 128‑test snapshot).
- GUI is running with the T9 trial folder (`/tmp/opencode/t9_trial/OY_rs1695_T9_270910Run01`), log clean.

### Active
- Batch auto-genotyping engine implemented in `genotyping.py:auto_genotype()`. 95/96 wells on the T9 plate reproduce the manual ground truth; only G09 diverges (documented intentional divergence). The engine raises a human‑readable *reason* when a call cannot be made (no quartet, weak sample, etc.) and records the physical channels used.
- Planned UI: a **Genotyping ▸ Auto‑genotype selected wells** menu item that loops over the current selection, builds a per‑well result table (well, call, frac, hom1‑hom2‑het1‑het2 areas, snrs, flags, reason), shows a progress dialog, and offers an **Export** button (CSV/JSON/XLSX). The table uses the same plumbing as the existing save‑peaks dialog (`_gen_save`, `save_table`) so rows feed directly into CSV/JSON/XLSX export without overwriting the manual pick/area state.
- The engine satisfies the test suite: 17/17 tests pass.
- Need: wire the menu action, results table/progress, export; then rebuild zip, commit/push, relaunch GUI.

### Blocked
- No hard blocker for the chosen basic scope; the only open design question is whether auto results share the existing pick table or use a separate mode — the current plan uses a separate table/export that does not interfere with manual rows.

## Handover — verified state at commit 17b6892
The engine is complete and validated, but **it is not reachable from the GUI yet**.
There is no `Auto‑genotype` menu item, so nothing to click; it is callable from
Python only. Do not assume the feature is usable in the app.

Verified headlessly against the 96‑well T9 plate (dye order ACTG):

```
wells scored        : 96
matched manual score: 95
disagreements       : 1  [('G09', 'het', 'hom-2', '')]
rows with a reason  : 1  [('H01', 'weakest sample duplex is only 32x the noise on Ch2')]
elapsed             : 0.3s total, 0.003s per well
het ai-flagged      : 2
```

- G09 is the one **documented intentional** divergence: its heteroduplex sits at
  the detection floor, so it is called hom‑2 on purpose (see `KNOWN_DIVERGENCE`
  in `tests/test_rs1695.py`).
- H01 has no sample and correctly reports *why* rather than guessing.
- The measured sigmas track `tests/data/rs1695_measured.csv` to ~1%
  (A12: 74.2 / 3300.5 / 73.0 / 85.7 measured vs 74.2 / 3300.5 / 73.0 / 87.9 fixture).

### Reproduce without a display
```bash
cd Limoncello scorer          # the app directory of this repository
python3 -m pytest tests/test_auto_genotype.py -q          # 17 passed
python3 -m pytest -q                                      # full suite: 168 passed, 7 skipped
python3 -m pytest tests/test_massaction.py -q            # 27 passed, incl. real A01 well
python3 -m pytest tests/test_duplex_ui.py -q             # 3 passed (needs DISPLAY)
```

### Pitfall to avoid on the next machine
`/home/tv/t9_rs1695_analysis/data/t9raw.npz` stores wells in the plate's
**physical channel order**, while `doc.acgt` is always **A,C,G,T**. On this plate
the npz→acgt mapping is `[2, 3, 1, 0]`. Mixing these up silently points the
detector at the sample channel; it was the cause of an hour of false leads here.
Two earlier "detections" were also an artefact of passing an already-smoothed
trace into `find_is_quartet()`, which smooths internally a second time.

## Mass-action MF (CTCE)
Ported onto this line from the pre-V3 branch, so it sits beside the
auto-genotyping engine rather than being lost by the switch.

- `PeakPicker.mark_duplex()` tags **one** allelic position (the last-picked
  main's own cluster on its channel, so two positions in one well stay apart)
  with the same `HOM1`/`HOM2`/`HET1`/`HET2` names the internal standard gets.
- `PeakPicker.mass_action(rec)` returns
  `MF = (A_MUT + ½ × A_HET) / (A_WT + A_MUT + A_HET)` plus `ai`, `a_wt`, `a_mut`,
  `a_het`, `n_homoduplex`. The ½ term is what makes a clean heterozygote read
  0.5 rather than 0.25, and it carries the whole fraction below ≈5 % MF, where
  the mutant strands have all re-annealed and no mutant homoduplex is visible.
- Peak count decides the split: 4 = 2 homoduplexes + 2 heteroduplexes,
  **3 = 1 + 2** (a homozygote), 2 = 2 + 0. Naming a 3-peak position as a
  truncated standard would call a heteroduplex a homoduplex and halve a low MF.
- GUI: **Genotyping ▸ Tag duplex species for MF…**, an **MF** column in the pick
  table, and the MF formula in the Help text. Export gains `mf` and `ai`
  columns; `fraction` keeps its plain meaning and both stay blank until a
  position is tagged.
- Tests: `tests/test_massaction.py` (27, incl. the real OY A01 well reading
  MF 0.1223) and `tests/test_duplex_ui.py` (3, menu → status → table → CSV).

## Auto-genotyping in the GUI
Done. **Genotyping ▸ Auto-genotyping** now carries three commands:

- **Auto-genotype selected wells…** — scores every selected well, progress in
  the status bar every 8 wells, and the rows land in their own results table
  below the plot: well, call, frac, hom1/hom2/het1/het2, snr1‑4, flags and the
  reason for a no-call. A well that fails to load gets a `load-error` row
  rather than aborting the batch — one bad file must not take the other 95
  down. The table is separate from the manual pick/area table, and
  `self.pick_tree` is cleared so `_sync_pick_table` cannot paint over it.
- **Channel roles (standard / sample)…** — a small dialog for the two roles,
  each labelled with the base that channel holds under the run's dye order.
  This is the *next* item below, done.
- **Save auto-genotype table…** — writes the rows through the same
  `genotyping.save_table` the manual pick table uses, so both land in the same
  CSV/XLSX/JSON shapes.

Tests: `tests/test_auto_genotype_ui.py` (6) drives all of it through the real
window, including asserting the menu really offers the cascade, that the
same-channel pair is refused before any well is scored, and that the A-row
calls match `tests/data/rs1695_expected.csv`.

## Channel roles, first step
The standard/sample pair is now **asked for** rather than assumed, which is
the first half of making the engine dye-agnostic. `DEFAULT_IS_CHANNEL = 3` /
`DEFAULT_SAMPLE_CHANNEL = 2` are still the *initial* values in the dialog, and
`auto_genotype()` still takes them as keyword defaults — so the API is
unchanged and the 95/96 result is untouched — but a run on any other kit can
now be pointed at the right channels from the GUI, and setting the same
channel as both is refused before scoring (swapping them silently scores the
sample's own peaks as the standard).

What is **not** done yet, and is the reason the engine is still not general:
- The roles are one pair for the whole run. The four-dye case wants a duplex
  standard on *every* channel at once, which is a list of per-channel sets,
  not two integers.
- `PeakPicker.std` is a single flat `[(scan, name)]` set and `duplex_of()`
  matches on scan alone, so tagging more than one channel can cross-label
  peaks from different channels into one position. `length_bp` is likewise a
  single scalar. This is the shared root of both the manual MF path and
  `auto_genotype`, and it should become per-channel before either is trusted
  on a multiplexed plate.
- `DEFAULT_IS_CUT = 1900` is the T9 kit's injection-front offset.

## Next Move
1. Make `PeakPicker` hold a duplex set **per channel** (`std` keyed by column,
   `duplex_of` matching column as well as scan, `length_bp` per channel), then
   let `Tag duplex species for MF…` pick which channel it is tagging. This
   unblocks the four-dye question and the multiplex `duplex_of` collision.
2. Extend `auto_genotype()` to take a per-channel role list so one call can
   call every channel's assay, and add a `std_scans`/reason per channel.
3. Re‑run the 96‑well check, the full suite, rebuild `Limoncello scorer.zip`
   with a clean‑extract test, then commit and push **without force‑pushing**.

## Known limitations to carry forward
- `DEFAULT_IS_CUT = 1900` is T9‑specific; the trace‑derived `cut` was deferred by choice.
- The standard/sample channel pair is a property of the **kit**, not of the run: standard Ch3 (T), sample Ch2 (C) on ACTG. It cannot be inferred from a trace, so it stays explicit in every result row.
- `SATELLITE_AMBIGUOUS_FRAC = 0.15` is still uncalibrated — no known minor‑allele plate. Treat its flags as "look at this", not as calls.

## Relevant Files
- `/home/tv/Limoncello-scorer/genotyping.py`: `auto_genotype()` (new), `find_is_quartet()`, `acgt_index_for_channel()`, constants `DEFAULT_IS_CHANNEL = 3`, `DEFAULT_SAMPLE_CHANNEL = 2`, `DEFAULT_IS_CUT = 1900`, `_noise_sigma()`, `_quartet_segments()`, `T9_MIN_DOMINANT_SIGMA = 40.0`.
- `/home/tv/Limoncello-scorer/tests/test_auto_genotype.py`: 17 tests, all passing — synthetic traces, channel defaults, ground‑truth 96‑well validation, failure paths.
- `/home/tv/Limoncello-scorer/scorer.py`: `t9_call()`, `t9_allele_fraction()`, `T9_MIN_DOMINANT_SIGMA = 40.0`.
- `/home/tv/Limoncello-scorer/sequence_analyzer.py`: Genotyping menu (`genotyping_m` at line 437‑480), `_gen_save()`, `_build_pick_table()`, `_sync_pick_table()`, existing mode handling.
- `/home/tv/Limoncello-scorer/Limoncello scorer.zip`: rebuilt release archive; clean extraction passes 128 tests.
- `/tmp/opencode/t9_trial/OY_rs1695_T9_270910Run01`: extracted 96‑well GUI trial data.
- `/home/tv/Nedlastinger/OY_rs1695_T9_270910Run01.zip`: authoritative source archive.
- `/home/tv/t9_rs1695_analysis/data/t9raw.npz`: cached T9 traces for validation.