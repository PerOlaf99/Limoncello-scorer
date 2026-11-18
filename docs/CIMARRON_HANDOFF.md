# Cimarron-like base caller: handoff (2026-11-18)

A pure-Python caller that reproduces Cimarron 3.12's base calls from raw
MegaBACE `.rsd` traces, so Limoncello can offer "Cimarron-like" calls without
Wine or the vendor DLL. This branch (`cimarron-clean-rewrite` on GitHub; locally `worktree-agent-a3d76a20c7a371f50`) is
**not merged into main**.

## Clean-room rules (keep them)
- **Black box only.** The targets are the DLL's outputs (`.esd`) on the operator's own plates.
- **Never read or port:**
  - the Ghidra decompile, `DLL_FUNCTION_MAP.md` or `CIMARRON_MASTER.md`;
  - CE's `cimarrontv.py`;
  - the source printed in patent EP0944739A1.
- **This branch is the clean rewrite ("option B").** Its implementer saw only two prose specs (`SPEC_signal.md`, `SPEC_calling.md`) plus black-box measurements. `CLEANROOM_rewrite.md` is the provenance log.
- **The phase-2 branches saw the patent's printed code and are never published:**
  - `worktree-agent-af94c3e4a66233f34` (signal);
  - `worktree-agent-abe6c3d36f0335285` (calling).
- Phase 1 (`worktree-agent-a77afb0535a44d93c`, black box only) is merged into this branch.

## Code
- `cimarron_basecaller/cimarron_signal.py`: raw trace → processed trace. This covers:
  - span;
  - baseline;
  - the separation matrix, estimated from the run;
  - cepstral deconvolution;
  - block delays;
  - flattening.
- `cimarron_basecaller/cimarron_calling.py`: processed trace → sequence and qualities. This covers:
  - windows and band detection;
  - fuzzy omission;
  - the gap splitter and the N rule;
  - stitching;
  - good-sequence limits.
- `cimarron_basecaller/cimarron_clone.py`: the entry points `call(...)` and `call_rsd(path)`.
- Tests: `tests/test_cimarron_{signal,calling,clone}.py`.

## Results (from `CLEANROOM_rewrite.md`)
"Agreement" is the share of Cimarron's bases reproduced. "Extra" is the share of extra bases we call.

| | Agreement | Extra | Position ≤ 1 scan |
|---|---|---|---|
| MB1000 M13, even wells (development), end to end | 0.953 | 0.031 | 0.53 |
| MB1000 M13, odd wells (held out, scored once), end to end | 0.947 | 0.034 | 0.55 |
| MB4000 (GRYPHON2, look only), end to end | 0.744 | 0.917 | 0.29 |
| Calling stage on Cimarron's own processed trace (m13 even / odd) | 0.982 / 0.983 | 0.015 / 0.019 | 0.98 |

**The calling stage is essentially done; the remaining gap is the signal stage.**

## Next steps, in order (the operator wants MB1000 ≥ 99 % first, then MB4000)
1. **Late-read sharpness.** About 75 % of the errors are in the last 30 % of the read, where the spacing falls to about 6 scans and our trace is less resolved than Cimarron's.
   - Every global sharpening tweak lay on one trade-off line, with agreement − extra ≈ 0.92.
   - The mechanism: a calling window that spans a sharp and a blurred segment fits one expected width, so the blurred bands get split.
   - Try width handling that varies along the read *inside* a calling window.
2. **Per-well common delay.** An oracle test with the delays taken from Cimarron's trace raised "position ≤ 1 scan" from 0.53 to 0.95. Find a predictor for it.
3. **Read start on the odd wells.** The median start error is 104 scans, against 38 on the even wells, and the odd wells miss 482 bases in the first 10 % of the read.
4. **MB4000.** Cimarron decimates many broad-peak wells there; emulate that. The optics differ from the MB1000's, so key the separation matrix by MACHINE ID.
5. **Wire into the GUI** as an optional "Cimarron-like" caller, only after review and a decision to merge.

## Data and tools (not in git; saved outside the repo)
Everything lives in `/home/tv/limoncello_work/cleanroom/` (it was moved from the session scratch folder before a restart):
- `data/m13/{rsd,esd}`: the 96-well MB1000 plate. Develop on the even columns only; the odd ones are held out.
- `data/mb4000/{rsd,esd}`: the second instrument.
- `specs/SPEC_signal.md`, `specs/SPEC_calling.md`: the specs the rewrite was written from.
- `tools/esd_read.py`, `evaluate.py`, `e2e.py`, `e2e_rsd.py`, `oracle.py`: the `.esd` reader and the scorers.
- `work_rewrite/`: the rewrite's scratch scripts.
  - Do **not** open `work_signal/` or `work_calling/`; they belong to the phase-2 agents that saw the patent.

The original targets also exist at:
- `/home/tv/electropherogram/MB1000_M13_DT/MB1000_M13_DT_Cp312_MD1/*.esd`;
- `99_archive/MB4000_DEMO_DATA`.

The specs and data are on this PC only. To continue on another PC, copy `/home/tv/limoncello_work/cleanroom/` across, apart from `work_signal/` and `work_calling/`.
