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
- Latest validation: 128 tests pass, compile and whitespace checks pass, 51‑file release zip clean‑extracts with 128 tests passing.
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
cd /home/tv/Limoncello-scorer
python3 -m pytest tests/test_auto_genotype.py -q          # 17 passed
python3 -m pytest -q                                      # full suite
```

### Pitfall to avoid on the next machine
`/home/tv/t9_rs1695_analysis/data/t9raw.npz` stores wells in the plate's
**physical channel order**, while `doc.acgt` is always **A,C,G,T**. On this plate
the npz→acgt mapping is `[2, 3, 1, 0]`. Mixing these up silently points the
detector at the sample channel; it was the cause of an hour of false leads here.
Two earlier "detections" were also an artefact of passing an already-smoothed
trace into `find_is_quartet()`, which smooths internally a second time.

## Next Move
1. Add the **Genotyping ▸ Auto‑genotype selected wells** menu command in `sequence_analyzer.py` (near `genotyping_m`, line ~437‑480), calling `g.auto_genotype()` over `self.selected` wells.
2. Add the results table and progress/status updates; columns are well, call, frac, hom1, hom2, het1, het2, snr1‑4, flags, reason, is_channel, sample_channel, std_scans. Reuse the existing `save_table` export path (`_gen_save`, `genotyping.save_table`).
3. Keep auto results separate from the manual pick/area table so the two do not overwrite each other.
4. Re‑run the 96‑well check, the full suite, rebuild `Limoncello scorer.zip` with a clean‑extract test, then commit and push **without force‑pushing**.

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