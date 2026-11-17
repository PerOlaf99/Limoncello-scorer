#!/usr/bin/env python3
"""De-novo MegaBACE M13 basecaller (best tuned configuration).

Reference-free: pure DSP + peak tracking. No reference sequence and no trained
model are used to call bases; the M13 genome is only used for scoring.

Given one .rsd file or a directory of .rsd files, writes one FASTA per well and
a combined FASTA of all reads.

Usage:
    python basecall.py --input /path/to/MB1000_M13_DT --out calls
    python basecall.py --input A01.rsd --out calls
    python basecall.py --input plate_dir --out calls --format fq   # FASTQ with calibrated Phred
"""
import argparse
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from cimarron_basecaller import track_bases
from cimarron_basecaller.phred import phred_qualities
from cimarron_basecaller.rsd_io import read_rsd, to_acgt_trace

# Best configuration on MB1000_M13_DT. Tuned against the project GOLDEN
# STANDARD (NCBI BLAST+ megablast vs M13mp18 M77815.1): beats Cimarron 3.12 on
# matched bases (+8.4%) and bit score (+1.6%). The position-profiled pull-back
# (0.008 -> 0.001 over the last 2/3 of the read) recovers the degraded 3' tail,
# and retuning the Wiener band filter (sigma_scale 1.05, noise_reg 0.06) cuts
# substitutions and lengthens the error-free runs without losing matched bases.
from configs import GOLDEN as WIN_CONFIG  # noqa: E402

# Fallback for the rare well whose profiled read runs away (tracker follows
# noise past the true end; those reads are implausibly long and have low mean
# base quality). Same stable scalar config the profiled caller is compared
# against; keeps every well alignable.
FALLBACK_CONFIG = dict(
    use_gaussian_reconstruction=True,
    gaussian_recon_segment_size=384,
    use_combined_channel_score=True,
    window_frac=(0.75, 1.25),
    local_norm_window=1800,
    channel_peak_bonus=1.0,
    pullback_weight=0.019,
    ema_alpha=0.10,
)
QUALITY_GATE = 2.0

# Precision configuration: optimises the longest error-free run and per-base
# identity instead of matched bases, accepting a shorter read. Same tracker as
# WIN_CONFIG with (a) strong Wiener regularization (noise_reg 0.128), (b) the
# spacing-anchor curve, and (c) a deep 38th-percentile quality trim. On the
# 96-well plate this gives 98.19% identity / 473 bp longest run vs Cimarron
# 3.12's 96.76% / 490.7 bp (golden mode: 94.30% / 363.9 bp, but more matched
# bases). Use `--mode precision` to select it.
from configs import PRECISION as WIN_CONFIG_PRECISION  # noqa: E402

PHRED_OFFSET = 33


def basecall_well(path, config=WIN_CONFIG, base_order="TGCA"):
    """Return (sequence, phred_quals) for one .rsd file.

    The qualities are calibrated Phred scores (see cimarron_basecaller.phred),
    priced with the model fitted for *config*; a read that drops to the
    fallback config is priced with the golden model.
    """
    rsd = read_rsd(path)
    trace, order = to_acgt_trace(rsd, base_order=base_order)
    model = "precision" if config == WIN_CONFIG_PRECISION else "golden"
    seq, heights, bands = track_bases(trace, base_order=order, **config)
    # QUALITY_GATE is on the tracker's peak heights, not on Phred.
    if len(heights) == 0 or float(np.mean(heights)) < QUALITY_GATE:
        seq, heights, bands = track_bases(trace, base_order=order, **FALLBACK_CONFIG)
        model = "golden"
    return seq, phred_qualities(bands, model)


def _write_fasta(fh, name, seq, width=60):
    fh.write(">%s len=%d\n" % (name, len(seq)))
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def _write_fastq(fh, name, seq, quals):
    fh.write("@%s len=%d\n%s\n+\n%s\n" % (
        name, len(seq), seq,
        "".join(chr(min(q, 93) + PHRED_OFFSET) for q in quals)))


def _collect(path):
    if os.path.isdir(path):
        return sorted(os.path.join(path, f) for f in os.listdir(path)
                      if f.lower().endswith(".rsd"))
    if path.lower().endswith(".rsd"):
        return [path]
    raise SystemExit("input must be a .rsd file or a directory of .rsd files")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help=".rsd file or directory")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--format", choices=("fa", "fq"), default="fa",
                    help="fa = FASTA (default), fq = FASTQ with calibrated Phred qualities")
    ap.add_argument("--mode", choices=("golden", "precision"), default="golden",
                    help="golden = matched-bases/bit-score optimum (default); "
                         "precision = longest-error-free-run/%%ID optimum")
    ap.add_argument("--base-order", default="TGCA",
                    help="plate dye order (default TGCA, the DYEnamic order)")
    args = ap.parse_args()
    config = WIN_CONFIG if args.mode == "golden" else WIN_CONFIG_PRECISION

    wells = _collect(args.input)
    if not wells:
        raise SystemExit("no .rsd files found in %s" % args.input)
    os.makedirs(args.out, exist_ok=True)

    ext = ".fasta" if args.format == "fa" else ".fastq"
    combined = os.path.join(args.out, "all_reads" + ext)
    total = 0
    with open(combined, "w") as cf:
        for path in wells:
            name = os.path.splitext(os.path.basename(path))[0]
            seq, quals = basecall_well(path, config, args.base_order)
            total += len(seq)
            if args.format == "fa":
                with open(os.path.join(args.out, name + ext), "w") as fh:
                    _write_fasta(fh, name, seq)
                _write_fasta(cf, name, seq)
            else:
                with open(os.path.join(args.out, name + ext), "w") as fh:
                    _write_fastq(fh, name, seq, quals)
                _write_fastq(cf, name, seq, quals)
            print("%s  %d calls" % (name, len(seq)))
    print("wrote %d wells (%d bases) to %s" % (len(wells), total, args.out))


if __name__ == "__main__":
    main()
