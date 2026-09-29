"""
Alignment-based scoring against ground truth, giving numbers analogous to
BLAST's "Query Cover" / "Percent Identity" so accuracy is measured
honestly rather than via loose substring heuristics.
"""
from __future__ import annotations
from Bio import Align


def score_against_truth(called: str, truth: str) -> dict:
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1

    alignments = aligner.align(called, truth)
    best = alignments[0]

    aligned_called, aligned_truth = str(best[0]), str(best[1])
    matches = sum(1 for a, b in zip(aligned_called, aligned_truth) if a == b and a != "-")
    aligned_len = len(aligned_called)
    identity = matches / aligned_len * 100 if aligned_len else 0.0
    query_cover = aligned_len / len(truth) * 100 if len(truth) else 0.0

    return {
        "identity_pct": identity,
        "query_cover_pct": query_cover,
        "aligned_len": aligned_len,
        "matches": matches,
        "score": best.score,
        "called_len": len(called),
        "truth_len": len(truth),
    }
