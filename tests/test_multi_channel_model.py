"""One standard, several sample channels.

Some runs carry more than one sample dye per well against a single internal
standard.  The standard is one physical object, so the plate's shape is the same
for every dye; only where the sample product sits relative to the standard
changes with the fluorophore.  ``per_channel_models`` fits the shape once and
gives each sample channel its own copy primed with its own product offsets.
These tests pin that split.
"""
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import genotyping  # noqa: E402

# base_order "ACTG": Ch1 -> col 0 (A), Ch2 -> col 1 (C), Ch3 -> col 2 (T),
# Ch4 -> col 3 (G).  The standard lives on Ch3, the samples on Ch1/Ch2 -- two
# different dyes, so two different shifts.
IS_COL = genotyping.acgt_index_for_channel("ACTG", 3)
CH1_COL = genotyping.acgt_index_for_channel("ACTG", 1)
CH2_COL = genotyping.acgt_index_for_channel("ACTG", 2)

SEED_SCANS = [2600, 2681, 2745, 2755]


def _gauss(y, centre, height, width=2.0):
    """Add one Gaussian peak to a channel, in place."""
    lo, hi = max(0, int(centre) - 6), min(len(y), int(centre) + 7)
    for i in range(lo, hi):
        y[i] += height * np.exp(-0.5 * ((i - centre) / width) ** 2)


def _wells(count=10, ch1_shift=25, ch2_shift=-11):
    """``count`` wells: a quartet on the IS channel, one product band per dye.

    The product is a single strong band on H1, so it is the loudest band in its
    segment and ``prime_plate_co_migration`` records its offset as the dye's
    shift.  Ch1 and Ch2 are given different shifts -- that difference is the
    whole point.
    """
    traces = {}
    for i in range(count):
        rng = np.random.default_rng(i + 1)
        acgt = rng.normal(0.0, 1.0, (3400, 4))
        for scan in SEED_SCANS:                       # the standard, on Ch3
            _gauss(acgt[:, IS_COL], scan, 800.0)
        _gauss(acgt[:, CH1_COL], 2600 + ch1_shift, 800.0)
        _gauss(acgt[:, CH2_COL], 2600 + ch2_shift, 800.0)
        traces[f"A{i:02d}"] = acgt
    return traces


def _fit():
    traces = _wells()
    seeds = {w: list(SEED_SCANS) for w in traces}
    models, found = genotyping.per_channel_models(
        traces, seeds, is_col=IS_COL, sample_channels=[1, 2], base_order="ACTG")
    return traces, seeds, models, found


def test_returns_one_model_per_sample_column():
    _t, _s, models, found = _fit()
    assert set(models) == {CH1_COL, CH2_COL}
    assert len(found) == 10                           # the standard was placed


def test_the_shape_is_shared_across_channels():
    _t, _s, models, _f = _fit()
    a, b = models[CH1_COL], models[CH2_COL]
    assert a.geometry is b.geometry                   # one fit, two dyes
    assert a.centers == b.centers
    assert a.merged_seen == b.merged_seen
    assert a.std_snr == b.std_snr
    assert a.std_heights == b.std_heights


def test_each_channel_learns_its_own_product_offset():
    _t, _s, models, _f = _fit()
    # H1 (duplex 0) is the loudest product band on both channels, so its learned
    # offset is that dye's shift: +25 on Ch1, -11 on Ch2.
    assert models[CH1_COL].co_migration_median(0) == 25.0
    assert models[CH2_COL].co_migration_median(0) == -11.0


def test_channels_do_not_share_their_offset_lists():
    _t, _s, models, _f = _fit()
    models[CH1_COL].co_mig[0].append(999.0)
    assert 999.0 not in models[CH2_COL].co_mig[0]


def test_copy_detaches_the_mutable_references():
    traces, seeds, _models, _found = _fit()
    base = genotyping.semi_auto_plate_model(
        traces, seeds, is_col=IS_COL, samp_col=None)[0]
    clone = base.copy()
    assert clone.geometry is base.geometry            # shape may be shared
    assert clone.co_mig == [[], [], [], []]           # offsets start empty
    clone.std_snr.append(1.0)
    clone.co_mig[0].append(5.0)
    assert clone.std_snr != base.std_snr              # a real copy
    assert base.co_mig[0] == []                       # untouched


def test_the_standard_channel_cannot_be_a_sample():
    traces = _wells()
    seeds = {w: list(SEED_SCANS) for w in traces}
    try:
        genotyping.per_channel_models(
            traces, seeds, is_col=IS_COL, sample_channels=[3], base_order="ACTG")
    except ValueError as e:
        assert "standard channel" in str(e)
    else:
        raise AssertionError("the standard channel was accepted as a sample")
