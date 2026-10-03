# Sample peak co-migration: the missing check

Status: diagnosed, **fixes not yet implemented**. Recorded `89cddaee`.

## The observation

The internal standard is added to the same PCR as the sample, so the sample's
own product co-migrates with it. Measured on `OY_ABCC2_N10_170910Run01`, the
distance between a sample peak apex and the corresponding standard band is
**about 9-12 scans** (operator-confirmed ballpark). Well-behaved wells measure
-15 to -9; none of the correct homoduplexes sit more than ~15 scans off.

| well  | W1 apex vs IS1 | W2 apex vs IS2 | W1 snr | truth      | current call |
|-------|---------------|----------------|--------|------------|--------------|
| A01   | -14           | +129           | 1679   | hom-1      | hom-1        |
| D07   | -15           | +46            | 571    | hom-1      | hom-1        |
| H11   | -15           | +81            | 54     | hom-1      | hom-1        |
| F12   | **+12**       | -44            | 39     | no-call    | hom-1        |
| H12   | -13           | +58            | 8      | no-call    | hom-1        |
| H01   | -9            | **+82**        | 143    | hom-1      | **hom-2**    |

## Why this matters

`auto_genotype()` integrates whatever peak lies inside each window, and
`_quartet_segments()` sets those windows from midpoint to midpoint. Nothing
asks whether the peak is at the position the standard says it should be. With
`d1 = 81`, the H2 window spans ~170 scans while the expected H2 position is
fixed by the standard, so a strong unrelated peak falls inside it and is
measured as the second allele.

That is the whole of H01: its H1 apex co-migrates correctly (-9), but a
94.6-sigma peak at +82 from IS2 was read as H2, giving a false hom-2. The
operator reads H01 as a weak hom-1.

## Three separate defects, not one

1. **H01 - no apex tolerance.** A hom call needs the measured apex within
   tolerance of the corresponding IS band, not merely inside the window.
   H01's +82 is far outside the 9-12 ballpark.

2. **H12 - `t9_call` lets a heteroduplex stand in for the homoduplex.**
   In `scorer.t9_call`, `dom = max(sigmas)` (scorer.py:422) is taken over all
   four sigmas. H12's homoduplex reads 8.5 sigma -- below the 40 sigma dominant
   floor -- but het1 reads 44, so `dom` clears the floor and a hom-1 is
   returned. For a *hom* call the homoduplex itself has to clear the floor. A
   strong heteroduplex with a missing homoduplex is the signature of a bulge or
   a nonspecific product, not a genotype. H12 is the only marked well in this
   run that currently calls hom with a sub-threshold homoduplex.

3. **H11 - unflagged post-IS spike.** The call (hom-1) is correct. There is a
   16-sigma peak at scan 3139, +469 past IS H1, which is noise after the IS and
   deserves a flag. Informational only, does not change the call.

## F12: still needs a decision

F12 is genuinely ambiguous and is **not** resolved by the co-migration check
alone. It has a real peak at 2652 (W1 apex +12 from IS1, so it *does*
co-migrate) which the operator describes as "very small". The 1025-sigma peak at
scan 2322 is not the product -- confirmed -- but no current rule looks at
off-IS peaks at all.

So F12 turns on a criterion the co-migration check does not supply: whether a
single co-migrating fragment that is small, with nothing at H2 and nothing at
either heteroduplex, is a no-call. Two readings give different code:

- **Reject on co-migration alone**: one co-migrating fragment plus a
  dominant off-IS peak means the well is not the expected product -> no-call.
- **Reject on amplitude**: F12's co-migrating peak is small in absolute terms
  relative to plate, independent of position.

F12 at 39 sigma is also just under `T9_MIN_DOMINANT_SIGMA = 40`, so a floor
change would move it -- but changing that constant would affect every well on
every plate and should not be done to accommodate one well. **Left untouched
pending the operator's call.**

## Not yet done

- No code change has been made for any of the three defects above.
- `tests/test_co_migration.py` does not exist yet.
- Once implemented, re-run the 96-well check and the full suite; ABCC2_N10 is
  currently 92/96 and the expectation is 94/96 with F12 still open.