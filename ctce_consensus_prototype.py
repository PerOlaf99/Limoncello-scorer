"""PROTOTYPE (2026-11-23): WT / mutant / no-PCR from peak pattern, with no settings per fragment.
Each well is lined up on its own main peak, the plate's WT shape is the median of its most ordinary wells,
and a well's score is its extra peaks beyond that shape.  The no-PCR test is main-peak height against the plate.
Usage: ctce_consensus_prototype.py <run folder> <fragment> <plate>   (scored against the Lavage paper key)
See docs/CTCE_HANDOFF.md, 2026-11-23."""
import sys, glob, csv, numpy as np
sys.path.insert(0, "/home/per/Limoncello-scorer")
from pathlib import Path
import analyzer_core as ac
from scipy.ndimage import uniform_filter1d, median_filter
PRE, POST = 120, 350

def plate(folder):
    W = {Path(f).stem[:3]: np.asarray(ac.load_trace(Path(f)).acgt, float)[:, :4] for f in sorted(glob.glob(folder + "/*.rsd"))}
    ch = int(np.argmax(np.median([w.max(0) for w in W.values()], 0)))
    T = {}
    for k, w in W.items():
        s = w[:, ch] - median_filter(w[:, ch], 301)
        a = int(np.argmax(s[1000:])) + 1000                         # the well's own main peak
        if a - PRE < 0 or a + POST > len(s): continue
        T[k] = (s[a - PRE:a + POST], s[a])
    return T

def score(folder):
    T = plate(folder)
    H = np.array([h for _, h in T.values()])
    N = {k: x / max(h, 1) for k, (x, h) in T.items()}
    R = np.median(list(N.values()), 0)
    for _ in range(2):
        sc = {k: np.clip(uniform_filter1d(n - R, 5), 0, None)[PRE + 8:].sum() for k, n in N.items()}
        cut = np.percentile(list(sc.values()), 50)
        R = np.median([N[k] for k in N if sc[k] <= cut], 0)
    out = {}
    for k, n in N.items():
        extra = np.clip(uniform_filter1d(n - R, 5), 0, None)[PRE + 8:]
        out[k] = dict(excess=float(extra.sum()), maxex=float(extra.max()), rel_height=float(T[k][1] / np.median(H)))
    return out

def auc(pos, neg):
    return np.mean([p > q for p in pos for q in neg]) if pos and neg else float("nan")

if __name__ == "__main__":
    folder, frag, pl = sys.argv[1:4]
    key = {x["well"]: x["code"] for x in csv.DictReader(open(
        "/run/media/per/NVMe Disk/limoncello-data/answer_keys/arstad_lavage_scoring.csv"))
        if x["fragment"] == frag and (x["plate"] or "1") == pl}
    S = score(folder)
    g = lambda codes, f: [f(S[w]) for w in S if key.get(w) in codes]
    m = g(("MUT", "LOW", "VERYLOW"), lambda s: s["maxex"]); wt = g(("WT",), lambda s: s["maxex"])
    print(f"{frag} p{pl}: MUT/LOW vs WT  AUC {auc(m, wt):.3f} (n {len(m)}/{len(wt)})   "
          f"WT vs NOPCR by main-peak height AUC {auc(g(('WT',), lambda s: s['rel_height']), g(('NOPCR',), lambda s: s['rel_height'])):.3f}")
    top = sorted(S, key=lambda w: -S[w]["maxex"])[:15]
    print("   top 15 by extra peaks:", " ".join(f"{w}:{key.get(w,'?')}" for w in top))
