"""Leftover crosstalk after a dye-overlap correction: for every clean
single-dye peak (rank-one window) of a run, the strongest other channel
relative to the dye's own after correction.  Lower is better; a perfect
matrix leaves ~0.  Usage: python matrix_check.py RUN_FOLDER [base_order]"""
import sys, warnings
from pathlib import Path
import numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from analyzer_core import load_trace  # noqa: E402
import spectral  # noqa: E402


def leftover(chans, m, cut=1000):
    vals = []
    for x in chans.values():
        y = m.correct(x) if m is not None else x
        for v, _w in spectral.peak_signatures(x, cut):
            pass
        sig = spectral.peak_signatures(y, cut)
        for v, _w in sig:
            o = np.sort(v)[-2]          # strongest channel after its own
            vals.append(o)
    return np.median(vals), np.percentile(vals, 90), len(vals)


def main(folder, order=None):
    folder = Path(folder)
    chans = {f.stem: np.asarray((load_trace(f, base_order=order) if order else load_trace(f)).channels, float)
             for f in sorted(folder.glob("*.rsd"))}
    mats = {"none": None, "estimated": spectral.matrix_from_peaks(chans, cut=1000)}
    for name, m in spectral.run_matrices(folder).items():
        mats["run " + name] = m
    for extra in sys.argv[3:]:
        for name, m in spectral.read_xtalk_tdb(extra).items():
            mats["file " + name] = m
    for name, m in mats.items():
        med, p90, n = leftover(chans, m)
        print(f"{folder.name:30s} {name:16s} leftover median {med:.3f}  p90 {p90:.3f}  ({n} peaks)")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 and sys.argv[2] != "-" else None)
