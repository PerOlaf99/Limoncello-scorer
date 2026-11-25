"""Step 1: per run, which (ladder channel, ladder kit) fits best.  12 wells spread over the plate,
every physical channel x every ROX-type kit, ladder fit only (no sample sizing)."""
import sys, csv, glob, os, warnings, multiprocessing as mp
warnings.filterwarnings("ignore")
sys.path.insert(0, "/home/per/Limoncello-scorer")
from pathlib import Path
CAT = "/run/media/per/NVMe Disk/limoncello-data/catalog/runs_current.csv"
OUT = sys.argv[1]
KITS = ["et400r", "et550r", "et900r", "promega_ils600", "promega_ils500", "genescan500_rox",
        "genescan400hd_rox", "mapmarker1000", "beckman_ss400"]

def fit(fs, y, lad):
    pk = fs.detect_peaks(y, min_prominence_frac=fs.DEFAULT_LADDER_PROMINENCE_FRAC, min_height_frac=fs.DEFAULT_LADDER_HEIGHT_FRAC)
    noise = fs._baseline_noise(y); top = max((p.height for p in pk), default=0)
    if not pk or noise <= 0 or top < fs.MIN_LADDER_SNR * noise:
        return dict(anchors=0, top=top, noise=noise)
    cand = fs.ladder_candidates(pk, len(lad.lengths))
    m = fs.align_peaks([p.scan for p in cand], lad.lengths)
    q = fs.ladder_fit_quality([s for s, _ in m], [b for _, b in m])
    hs = {p.scan: p.height for p in cand}
    mh = [hs[s] for s, _ in m]
    return dict(anchors=len(m), rms=q["rms_error_bp"], sp=q["rms_error_spacing"], top=top, noise=noise,
                band_med=float(np.median(mh)) if mh else 0)

def one(row):
    import fragment_sizing as fs
    from analyzer_core import load_trace, n_channels
    global np; import numpy as np
    files = sorted(glob.glob(row["path"] + "/*.rsd"))
    if not files: return []
    pick = files[::max(1, len(files) // 8)][:8]
    try:
        m = fs.plate_matrix_for([Path(f) for f in files], prefer="auto")
    except Exception:
        m = None
    out = []
    lads = {k: fs.load_ladder(k) for k in KITS}
    for f in pick:
        try:
            d = load_trace(Path(f))
            raw, corr = fs._channels_for_sizing(d, d.base_order or "ACTG", m)
        except Exception as e:
            continue
        for c in range(corr.shape[1]):
            y = corr[:, c]
            for k, lad in lads.items():
                r = fit(fs, y, lad)
                out.append(dict(run=row["name"], kind=row["kind"], well=Path(f).stem, channel=c + 1, ladder=k,
                                n_bands=len(lad.lengths), path=row["path"], **r))
    return out

if __name__ == "__main__":
    rows = [x for x in csv.DictReader(open(CAT)) if x["kind"].startswith(("fragment_sizing", "fragment_or"))]
    keys = ["run", "kind", "well", "channel", "ladder", "n_bands", "anchors", "rms", "sp", "top", "noise", "band_med", "path"]
    with open(OUT, "w", newline="") as fh, mp.Pool(5) as pool:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for i, res in enumerate(pool.imap_unordered(one, rows), 1):
            w.writerows(res); fh.flush(); print(f"{i}/{len(rows)}", flush=True)
