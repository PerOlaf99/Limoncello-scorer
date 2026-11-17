"""
cimarron_basecaller: a from-scratch base-calling implementation for
MegaBACE 1000 CE raw trace data, informed by (but not limited to) the
publicly disclosed algorithm in EP0944739A1 (Cimarron base caller
patent, application withdrawn -- no patent was ever granted).

Entry point: `track_bases` (empirically tuned, validated against real
Cimarron 3.12 ground truth -- see BEST_BASECALLER_README.md). The earlier
patent-literal pipeline and the DP / peak-picking callers scored far worse
and were removed; they remain in the git history.
"""

from .spacing_caller import track_bases, TrackedBase

__all__ = ["track_bases", "TrackedBase"]
__version__ = "0.3.0"
