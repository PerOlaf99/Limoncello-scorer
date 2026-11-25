"""Step 2: size every fragment run with ladder="auto" + plate rescue.
Per well: chosen kit/channel, clean fit, rescued, LOO RMS.  Per band: scan, height, residual (for table calibration)."""
import sys, csv, glob, os, warnings, multiprocessing as mp
warnings.filterwarnings("ignore")
sys.path.insert(0, "/home/per/Limoncello-scorer")
CAT = "/run/media/per/NVMe Disk/limoncello-data/catalog/runs_current.csv"
OUT, BANDS = sys.argv[1], sys.argv[2]

def one(row):
    import fragment_sizing as fs
    files = sorted(glob.glob(row["path"] + "/*.rsd"))
    W, B = [], []
    try:
        res = fs.size_files(files, "auto", sample_channels=[1])   # ladder only matters; one cheap sample channel
    except Exception as e:
        return [dict(run=row["name"], kind=row["kind"], well="", detection="run error: " + str(e)[:100], path=row["path"])], []
    for r in res:
        nb = len(r.ladder.lengths) - len(r.ladder.exclude)
        ok = fs.ladder_fit_ok(r.quality, nb)
        W.append(dict(run=row["name"], kind=row["kind"], well=r.well, ladder=r.ladder.name, channel=r.ladder_channel,
                      detection=r.detection, clean=int(ok), rescued=int(any("plate's other" in w for w in r.warnings)),
                      anchors=r.quality.get("anchors"), n_bands=nb, rms_bp=r.quality.get("rms_error_bp"),
                      rms_spacing=r.quality.get("rms_error_spacing"), path=row["path"]))
        if ok:
            for s in r.standards:
                B.append(dict(run=row["name"], well=r.well, ladder=r.ladder.name, band_bp=s["band_bp"], scan=s["scan"],
                              height=s["height"], width=s["width"], residual_bp=s["residual_bp"]))
    return W, B

if __name__ == "__main__":
    rows = [x for x in csv.DictReader(open(CAT)) if x["kind"].startswith(("fragment_sizing", "fragment_or"))]
    if len(sys.argv) > 3:                      # optional: only the runs named in this file
        keep = set(open(sys.argv[3]).read().split("\n")) - {""}
        rows = [x for x in rows if x["name"] in keep]
    wk = ["run", "kind", "well", "ladder", "channel", "detection", "clean", "rescued", "anchors", "n_bands", "rms_bp", "rms_spacing", "path"]
    bk = ["run", "well", "ladder", "band_bp", "scan", "height", "width", "residual_bp"]
    with open(OUT, "w", newline="") as f1, open(BANDS, "w", newline="") as f2, mp.Pool(int(os.environ.get("NP", 5))) as pool:
        w1 = csv.DictWriter(f1, fieldnames=wk); w1.writeheader()
        w2 = csv.DictWriter(f2, fieldnames=bk); w2.writeheader()
        for i, (W, B) in enumerate(pool.imap_unordered(one, rows), 1):
            w1.writerows(W); w2.writerows(B); f1.flush(); f2.flush()
            print(f"{i}/{len(rows)}", flush=True)
