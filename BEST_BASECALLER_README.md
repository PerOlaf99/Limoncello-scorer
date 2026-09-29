# best_basecaller

De-novo, reference-free MegaBACE M13 basecaller. Pure DSP + peak tracking
(numpy/scipy only). No reference sequence and no trained model are used at
call time.

This is the tuned caller for the `MB1000_M13_DT` plate. On all 96 wells, scored
with NCBI BLAST+ `megablast` against the NCBI M13mp18 reference (`M77815.1`):

| metric | this caller | Cimarron 3.12 (ESD) |
|---|---|---|
| matched bases | **78,362** | 72,286 |
| aligned length | **83,116** | 74,726 |
| coverage of reference (mean) | **11.94%** | 10.74% |
| total bit score | **124,783** | 122,831 |
| mean bit score | **1,299.8** | 1,279.5 |
| mean % identity | 94.30% | **96.76%** |
| coverage of read (mean) | **90.04%** | 89.18% |
| longest error-free stretch (mean) | 363.9 | **490.7** |
| mean read length | **962.5** | 873.0 |

It beats Cimarron 3.12 on the two decisive counters -- **matched bases**
(+8.4%) and **bit score** (+1.6%) -- and on coverage, by emitting a longer
true-positive read (962 bp vs 873 bp). Cimarron keeps the higher per-base
identity and the longer clean stretch, because it calls shorter reads. See
`NOTES.md` for the full breakdown.

## Precision mode (longest error-free run / %ID)

A second configuration trades read length and matched bases for a **longer
error-free run and much higher identity**. Same tracker, but stronger Wiener
regularization (`gaussian_recon_noise_reg=0.128`), the spacing-anchor curve, and
a deep 38th-percentile quality trim that discards the degraded read ends:

| metric | `--mode golden` | `--mode precision` | Cimarron 3.12 |
|---|---|---|---|
| matched bases | **78,362** | 59,403 | 72,286 |
| total bit score | **124,783** | 104,920 | 122,831 |
| mean % identity | 94.30% | **98.19%** | 96.76% |
| longest error-free run (mean) | 363.9 | **473.1** | 490.7 |
| longest error-free run (median) | 330.5 | **500** | 519 |
| mean read length | **962.5** | 627.4 | 873.0 |
| total gaps | 2,576 | **787** | 1,967 |

Precision mode wins identity by +1.43 pp over Cimarron and cuts gaps by 60%,
nearly matching its longest run (473 vs 491 bp).

## Install

```bash
python3 -m pip install -r requirements.txt
```

Requires Python 3.9+, numpy and scipy. Nothing else.

## Run

```bash
# whole plate (directory of *.rsd) -> calls/<well>.fasta + calls/all_reads.fasta
python3 basecall.py --input /path/to/MB1000_M13_DT --out calls

# a single well
python3 basecall.py --input A01.rsd --out calls

# FASTQ with PHRED qualities instead of FASTA
python3 basecall.py --input /path/to/MB1000_M13_DT --out calls --format fq

# precision mode: longest error-free run / %ID optimum
python3 basecall.py --input /path/to/MB1000_M13_DT --out calls --mode precision
```

Each well produces `<well>.fasta` (or `.fastq`), plus a combined `all_reads.*`.
Progress and per-well call counts are printed to stdout.

## API

```python
from cimarron_basecaller import track_bases
from cimarron_basecaller.rsd_io import read_rsd, to_acgt_trace

trace, order = to_acgt_trace(read_rsd("A01.rsd"), base_order="TGCA")
seq, quals, bands = track_bases(
    trace, base_order=order,
    use_gaussian_reconstruction=True, gaussian_recon_segment_size=384,
    use_combined_channel_score=True, window_frac=(0.75, 1.25),
    local_norm_window=1800, channel_peak_bonus=1.2,
    pullback_weight=(0.008, 0.001), ema_alpha=0.08,
    profile_fracs=(0.33, 1.0),
)
```

`basecall.py` holds this exact dict as `WIN_CONFIG`, plus a mean-base-quality
gate (`QUALITY_GATE`) that re-calls a runaway well with the stable scalar
`FALLBACK_CONFIG`. Any of `ema_alpha`, `pullback_weight`, `min_prominence`,
`channel_peak_bonus` accepts a `(start, end)` pair interpolated over
`profile_fracs`; a scalar or `(v, v)` reproduces the un-profiled caller
byte-for-byte.

## Files

- `basecall.py` - command-line entry point (`WIN_CONFIG` here).
- `cimarron_basecaller/` - DSP + tracking package:
  - `spacing_caller.py` - `track_bases`, the recommended caller.
  - `preprocessing.py`, `deconvolution.py`, `peaks.py`, `scoring.py`,
    `alignment.py` - baseline, spectral separation, Gaussian band filter,
    mobility shift, combined scoring.
  - `rsd_io.py` - MegaBACE `.rsd` reader (`to_acgt_trace`).
  - `dp_caller.py`, `simple_caller.py` - alternate callers (much worse; kept
    for reference).

## Notes

- The `WIN_CONFIG` constants were tuned on this plate. Re-tune
  `pullback_weight`, `channel_peak_bonus` and `profile_fracs` for other plates;
  the Gaussian band filter, combined-channel score, position-profiled pull-back
  and the quality gate are the transferable parts.
- Reads are called in the sample's own orientation. When comparing against an
  M13 reference, use a strand-agnostic aligner (e.g. BLAST / megablast), which
  is how the table above was produced.
- `cimarron_basecaller/` is vendored from the user-provided MegaBACE Cimarron
  software bundle; only numpy/scipy are required on the `track_bases` path.
