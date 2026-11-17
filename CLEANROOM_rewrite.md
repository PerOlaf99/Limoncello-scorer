# Clean-room provenance log: Cimarron 3.12 clone, rewrite (signal + calling chain)

Modules: `cimarron_basecaller/cimarron_signal.py`, `cimarron_basecaller/cimarron_calling.py`,
`cimarron_basecaller/cimarron_clone.py` (entry points `call`, `call_rsd`);
tests `tests/test_cimarron_signal.py`, `tests/test_cimarron_calling.py`,
`tests/test_cimarron_clone.py`. Phase 1 (span heuristics, empty-capillary rule)
is logged in `CLEANROOM.md`.

The implementer of this rewrite never saw the patent EP0944739A1 or any
source code printed in it. Everything here was written from two functional
specifications (prose, tables, maths) and from black-box measurements of
Cimarron's outputs.

## Sources consulted

| Source | Use |
|---|---|
| `cleanroom/specs/SPEC_signal.md` | Design of the signal stage (span, baseline, separation table and run estimate, cepstral deconvolution, block delays, flattening). |
| `cleanroom/specs/SPEC_calling.md` | Design of the calling stage (windows, band detection, expected curves, fuzzy omission, gap splitter, N rule, qualities, stitching with fallback, end rule, good-sequence limits). |
| `cleanroom/data/m13/{rsd,esd}` | Raw traces and Cimarron 3.12 results. Development on the 48 EVEN-column wells only; the 48 odd wells were evaluated once at the end. |
| `cleanroom/data/mb4000/{rsd,esd}` | Final look only (different instrument, GRYPHON2); no tuning. |
| `cleanroom/tools/esd_read.py`, `evaluate.py`, `e2e.py`, `oracle.py` | Read only; `esd_read.py` copied to the scratch folder, `e2e.py` adapted there to call `cimarron_clone.call_rsd`. |
| This repository (`main`) and the phase-1 branch `worktree-agent-a77afb0535a44d93c` (merged) | `rsd_io.read_rsd`; phase-1 empty-capillary rule `has_signal`; the clone's interface. |
| General literature | Homomorphic (cepstral) deconvolution and liftering; Mamdani-style fuzzy inference with centroid defuzzification; least-squares polynomial fits. |

Not consulted: the patent in any form, the branches/worktrees named as
forbidden, `cleanroom/work_signal/`, `cleanroom/work_calling/`, any other
scratchpad folder, any binary, disassembly or decompilation. No forbidden
material was encountered. (One development command that would have read the
odd wells' quality values was blocked before it ran; it was rerun on even
wells only.)

Scratch scripts (not part of the repo) live in `cleanroom/work_rewrite/`:
`common.py`, `sig_eval.py`, `flat_test*.py`, `call_eval.py`, `cache.py`,
`classify.py`, `feat_cmp.py`, `sp_bias.py`, `spec_ratio.py`, `xtalk.py`,
`delay_pred.py`, `e2e.py`.

## Design decisions (one line each: source or measurement)

### Signal stage (`cimarron_signal.py`)
- Frame: output row i = raw row start_tp - 1 + i, end_tp - start_tp rows — SPEC_signal §1.
- Span (start/end rules) — SPEC_signal §2 (same onset rule as phase 1).
- Baseline: moving mean 101 of moving min 101, span only, edges repeated — SPEC_signal §3.
- Separation applied as M @ x (rows = outputs); PLAYER1 table — SPEC_signal §4.
- **Deviation: the separation matrix is estimated from the run for known instruments too** (table only as seed; `ESTIMATE_KNOWN`). Measured end to end, m13 even, Cimarron's span: run estimate 0.9403 / extra 0.0199 vs PLAYER1 table 0.9314 / 0.0226 vs each well's own Cimarron matrix (read from its .esd) 0.9338 / 0.0213; read-middle errors drop from ~250 to ~30.
- Per-segment estimation of the matrix: tried, worse at the read end (0.927); not used.
- Deconvolution: 2048-row segments, step 1900, end-aligned last one, ownership by overlap midpoints, cepstral lifter 12..22, Gaussian low-pass with width W — SPEC_signal §5.2-5.3.
- W = 1210.8 / spacing, clipped to [40, 400] — SPEC_signal §5.4.
- **Deviation: the spacing-estimation pass uses W = 170 instead of 135.** Measured: with 135 the band-count spacing is 7 % / -3 % / 5 % / 23 % too long in segments 0 / 1 / 2 / last (merged bands late in the read), so the last segment is too blurred; with 170 the bias is +2 % / -5 % / -3 % / +6 %. End to end (Cimarron's span) 0.9421 / 0.0182 -> 0.9555 / 0.0336.
- Block delays and PLAYER1 delay table — SPEC_signal §6.
- Flattened output: lambda = Gaussian(100) of the 15-row running max of the envelope, 0.60 * internal / lambda, channel gains/offsets, clip — SPEC_signal §7.
- Outputs `internal`, `trace` (flattened), `divisor` (lambda), `gain`, `offset`, `separated`, `estimated` — SPEC_signal §8 / task.
- `unflatten` inverts `flatten` exactly where neither clip was active (test found the negative-offset case: a clipped inner value still gives flattened = -offset > 0).
- Footer metadata and dye order parsed from the .rsd footer's tag records (same tag layout as the .esd, observed in the files).

### Settling the flattening model (the two specs disagree)
Target: Cimarron's stored trace in the .esd. Measurements on m13 even, with our
internal trace and Cimarron's span:
- Per-peak amplitude: log(stored peak / our internal peak, local max within +-4 rows), de-meaned per well, regressed on log divisor. 257-tap box mean of the envelope (zero outside the trace, /129): residual variance 0.0441, free slope -1.08 (-1.00 in the read middle). Smoothed running maximum: 0.0478, slope -0.93 (-0.81 middle). Box with repeated edges: 0.0490. -> Cimarron's flattening divisor is the calling spec's 257-tap box mean of the envelope.
- Stored-envelope invariants: CV of the 257-scan mean 0.026 vs CV of the smoothed running max 0.036 (box mean is the flatter, i.e. normalised, quantity).
- The fixed-iteration inversion (SPEC_calling §3.3) is stable on Cimarron's trace but localises on our running-max flattened trace (inverted envelope 6x the mean in one region; agreement 0.58): only a trace whose 257-scan envelope mean is constant inverts sensibly — again the box model.
- Whole-trace correlation with the stored trace, per-channel gain/offset fitted: running max 0.801, box 257 0.783 (box 129: 0.793). With our imperfect internal trace the running max reproduces the stored trace better, because the box mean depends on our (too broad) peak shapes.
- Decision: the model of what Cimarron does is the 257-tap box mean (documented); the emitted `flattened` trace keeps SPEC_signal §7's running-max divisor because it matches the stored trace better (acceptance metric), and the calling stage is fed the unflattened trace, so no estimate of the flattening is needed in the chain.

### Integration
- **The caller is fed `caller_trace` = flattened * lambda / 0.60**: the exact inverse of the lambda division, keeping the channel gains, the offsets that suppress small cross-talk bumps (A offset 0.11) and the clipping. Measured: internal (clipped) 0.9403 / 0.0199, caller_trace 0.9421 / 0.0182; our flattened trace through the §3.3 inversion 0.58. Cross-talk at band peaks (other/dominant channel) in our internal trace is 2-4x Cimarron's in the G and A channels; Cimarron's stored trace has those bumps suppressed and the calling rules were fitted on it.
- Amplitude trend is not critical: internal / box-flattened-then-inverted / oracle trend from the stored trace all within 0.930-0.932.
- Scale: the caller multiplies its input by 0.5 / mean envelope — SPEC_calling §2.

### Calling stage (`cimarron_calling.py`)
- All rules, breakpoints and constants as SPEC_calling §4-15; fuzzy centroids on a fine grid (trapezoid rule).
- **Spec inconsistency: quality centroid domain.** §13.3 says "centroid ... over the domain 0...1" but gives rule 9 alone = 99.0 and rule 8 alone = 87.167, which are the centroids of the untruncated triangles; with the domain 0..1 rule 9 alone gives 95.8. Using the triangles' full extent (-0.1..1.105): exact-quality agreement on Cimarron's own trace 0.043 -> 0.782 (spec: 0.789).
- Band rebuild ("make touching"): ties move the left band's right boundary; width = r - l (raw list) and r - l + 1 after a rebuild — §7.
- Gap splitter for the first pair: r_prevgap uses the first band's (median) gap as 0 — §11.1 ambiguity, read as "0 for the first band".
- GC richness uses the up-to-5 bands before band j and their left gaps — §11.1 ambiguity.
- An inserted band that lands on an existing peak position is dropped; positions are forced strictly increasing after stitching (needed for the interface).
- Stitching with the no-match fallback — §4.2 (required); never discards the rest of the read.
- On Cimarron's own trace (flattened=True): m13 even agreement 0.9824 / extra 0.0150 / pos<=1 0.982 / quality exact 0.782 — the spec's §16 reference (0.9823 / 0.0149 / 0.982 / 0.789) is reproduced. Good-sequence limits from Cimarron's qualities: 47/48 even wells exact.

### Clone (`cimarron_clone.py`)
- `call(trace, base_order)` keeps the evaluator interface; the columns are put into the PLAYER1 dye order (TGCA) and processed without metadata (separation is estimated anyway, PLAYER1 delays/gains as dye-chemistry default — SPEC_signal §9).
- N is returned as "N" (Cimarron's sequence contains N); the phase-1 test that forbade letters outside base_order was relaxed accordingly, and the phase-1 test pinning the old GOLDEN-based internals was removed with those internals.
- Empty capillaries return no calls — phase 1 (`CLEANROOM.md`).
- `call_rsd(path)` reads the .rsd, its machine id and channel order.

## Results
(see the final section, filled in at the end of the work)
