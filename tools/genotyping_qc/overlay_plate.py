"""Overlay a plate's wells before and after aligning on the standard.

    python overlay_plate.py RUN_FOLDER OUT.png [method] [base_order]

Top: the sample channel of every well as recorded.  Bottom: the same wells
on the plate's axis (snp_genotyping.align_plate), coloured by call, with the
reference standard bands marked."""
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from analyzer_core import load_trace  # noqa: E402
import snp_genotyping as sg  # noqa: E402

COLORS = {"hom1": "tab:blue", "hom2": "tab:red", "het": "tab:green", "no-call": "0.75"}


def main(folder, out, method="piecewise", order="ACTG"):
    folder = Path(folder)
    docs = {f.stem: load_trace(f, base_order=order) for f in sorted(folder.glob("*.rsd"))}
    rows = sg.genotype_plate(docs, base_order=order, run_name=folder.name)
    arrs = {w: np.asarray(d.acgt, float) for w, d in docs.items()}
    aligned, ref = sg.align_plate(arrs, rows, method)
    ps_note = next(iter(rows.values()))["plate_note"]
    samp = None
    for ch in range(1, 5):
        if f"sample on Ch{ch}" in ps_note:
            samp = sg.acgt_index_for_channel(order, ch)
    samp = 1 if samp is None else samp
    fig, ax = plt.subplots(2, 1, figsize=(15, 9))
    lo, hi = int(ref[0] - 150), int(ref[-1] + 150)
    for w, a in arrs.items():
        c = COLORS.get(rows[w]["call"], "0.75")
        ax[0].plot(sg.flatten(a[:, samp]), color=c, lw=0.5, alpha=0.6)
        if w in aligned:
            ax[1].plot(np.arange(lo, hi), sg.flatten(aligned[w][:, samp])[lo:hi], color=c, lw=0.5, alpha=0.6)
    for r in ref:
        ax[1].axvline(r, color="k", ls=":", lw=0.8)
    allstd = [s for v in sg.std_from_rows(rows).values() for s in v]
    ax[0].set_xlim(min(allstd) - 150, max(allstd) + 150)
    ax[1].set_xlim(lo, hi)
    for a_ in ax:
        a_.set_ylim(0, None)
    ax[0].set_title(f"{folder.name}: sample channel as recorded ({len(arrs)} wells)", fontsize=10)
    ax[1].set_title(f"aligned on the standard ({method}); dotted: reference bands. "
                    "blue hom1, red hom2, green het, grey no-call", fontsize=10)
    plt.tight_layout()
    plt.savefig(out, dpi=70)


if __name__ == "__main__":
    main(*sys.argv[1:])
