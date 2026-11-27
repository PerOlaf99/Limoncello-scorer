"""Fragment sizing with no dye-overlap matrix vs the automatic one.

Per run: the kit/channel found, clean wells and median LOO error, both ways.
Runs: the cell ID runs that sized in the v8 validation (every 6th) plus a
spread of other kinds (where a matrix must not create a false ladder).

    python matrix_sizing_compare.py OUT.csv [n_other]
"""
import csv
import glob
import multiprocessing as mp
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, "/home/per/Limoncello-scorer")
CAT = "/run/media/per/NVMe Disk/limoncello-data/catalog/runs_current.csv"
V8 = "/run/media/per/NVMe Disk/limoncello-data/fragment_sizing/qc_20261125/v8_wells.csv"


def one(job):
    import fragment_sizing as fs
    name, kind, path = job
    files = sorted(glob.glob(path + "/*.rsd"))
    out = {"run": name, "kind": kind}
    for tag, matrix in (("none", "none"), ("auto", "auto")):
        try:
            res = fs.size_files(files, "auto", sample_channels=[1], matrix=matrix,
                                cache_dir=None)
        except Exception as e:      # noqa: BLE001
            out[tag] = f"error {str(e)[:60]}"
            continue
        clean, loo = 0, []
        for r in res:
            nb = len(r.ladder.lengths) - len(r.ladder.exclude)
            if fs.ladder_fit_ok(r.quality, nb) and "no size standard" not in r.detection:
                clean += 1
                if r.quality.get("rms_error_bp") is not None:
                    loo.append(r.quality["rms_error_bp"])
        det = res[0].detection if res else ""
        found = "none" if "no size standard" in det else f"{res[0].ladder.name} ch{res[0].ladder_channel}"
        out[f"{tag}_found"] = found
        out[f"{tag}_clean"] = clean
        out[f"{tag}_loo"] = round(float(np.median(loo)), 3) if loo else ""
        out[f"{tag}_matrix"] = res[0].matrix if res else ""
        out["wells"] = len(res)
    return out


def main():
    out_path = sys.argv[1]
    n_other = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    sized = []
    seen = set()
    for r in csv.DictReader(open(V8)):
        if r["run"] not in seen and "no size standard" not in r["detection"]:
            seen.add(r["run"])
            sized.append((r["run"], "cell_id", r["path"].strip()))
    jobs = sized[::6]
    cat = [r for r in csv.DictReader(open(CAT))
           if r["kind"] not in ("fragment_or_genotyping",) and int(r["n_rsd"] or 0) >= 12]
    rng = np.random.default_rng(1)
    pick = rng.choice(len(cat), size=min(n_other, len(cat)), replace=False)
    jobs += [(cat[i]["name"], cat[i]["kind"], cat[i]["path"]) for i in sorted(pick)]
    with mp.Pool(6) as pool:
        rows = []
        for k, row in enumerate(pool.imap_unordered(one, jobs)):
            rows.append(row)
            print(f"{k + 1}/{len(jobs)} {row['run'][:40]:40s} none: {row.get('none_found')} "
                  f"{row.get('none_clean')} {row.get('none_loo')} | auto: {row.get('auto_found')} "
                  f"{row.get('auto_clean')} {row.get('auto_loo')}", flush=True)
    keys = sorted({k for r in rows for k in r})
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
