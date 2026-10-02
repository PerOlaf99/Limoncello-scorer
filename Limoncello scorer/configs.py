"""Named basecaller presets for the Limoncello CE Analyzer.

These are the plate-validated configurations from the bundled best_basecaller
release (see ``basecall.py`` and ``BEST_BASECALLER_README.md``).  They are
consumed by :meth:`analyzer_core.AnalysisSettings.to_track_kwargs`, which passes
the dict straight to ``cimarron_basecaller.track_bases``.

Every key below is a *real* ``track_bases`` keyword, so a preset is authoritative
for that caller; the GUI only overrides individual knobs the user has actually
changed from their defaults.

The two decisive configurations (96-well MB1000_M13_DT, NCBI BLAST+ megablast vs
M13mp18 M77815.1):

* ``mb1000_length``  -- the matched-bases / bit-score optimum (longest true read).
* ``mb1000_accuracy`` -- the longest-error-free-run / %ID optimum (deep Q-trim).

The remaining names are kept so saved settings and menus stay valid:

* ``pos_bonus07`` / ``pos_profile`` are the legacy aliases of the two above.
* ``hz_soften`` is the length preset with milder Wiener regularization.
* ``mb4000_*`` currently alias the MB1000 presets (no MegaBACE 4000 spectral
  CHM is bundled with this release -- see the warning in the README).
"""

# Matched-bases / bit-score optimum (a.k.a. "golden" mode). The position-
# profiled pull-back recovers the degraded 3' tail; the retuned Wiener band
# filter (sigma_scale 1.05, noise_reg 0.06) cuts substitutions without losing
# matched bases.
GOLDEN = dict(
    use_gaussian_reconstruction=True,
    gaussian_recon_segment_size=384,
    gaussian_recon_noise_reg=0.06,
    gaussian_recon_sigma_scale=1.05,
    use_combined_channel_score=True,
    window_frac=(0.75, 1.25),
    local_norm_window=1800,
    channel_peak_bonus=1.2,
    pullback_weight=(0.008, 0.001),
    ema_alpha=0.08,
    profile_fracs=(0.33, 1.0),
)

# Longest-error-free-run / %ID optimum (a.k.a. "precision" mode): stronger
# Wiener regularization, the spacing-anchor curve and a deep 38th-percentile
# quality trim.
PRECISION = dict(GOLDEN)
PRECISION.update(
    gaussian_recon_noise_reg=0.128,
    use_spacing_anchor_curve=True,
    trim_quality_percentile=38.0,
)

# Length preset with a gentler band filter (historically the "mid hard-zone"
# softening).  The hard-zone deconvolution stage is not part of this caller
# build, so this is a conservative variant of GOLDEN rather than a distinct
# pipeline.
SOFT = dict(GOLDEN)
SOFT.update(gaussian_recon_noise_reg=0.10)

CONFIGS = {
    # unified instrument x mode presets
    "mb1000_accuracy": PRECISION,
    "mb1000_length": GOLDEN,
    "mb4000_accuracy": PRECISION,   # alias -- no MB4000 CHM bundled
    "mb4000_length": GOLDEN,        # alias -- no MB4000 CHM bundled
    # legacy aliases (kept so saved settings resolve)
    "pos_bonus07": PRECISION,
    "pos_profile": GOLDEN,
    "hz_soften": SOFT,
    # "raw_peaks" is handled directly in run_basecall (no track_bases call)
}
