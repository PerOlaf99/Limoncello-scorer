"""Where the sample product sits relative to its internal standard.

The internal standard is spiked into the same PCR as the sample, so a real
product co-migrates with it -- landing a short, repeatable distance from the
band it belongs to, not exactly on it.  On ABCC2_N10 that distance has a median
of 11 scans *before* the fitted centre, and every genuine hom-1 there falls
between -20 and -7.

That offset is load-bearing, because a duplex's segment is far wider than the
product.  With ``d1`` at 81 scans the H2 window spans ~170, so taking the
segment maximum hands the measurement to whatever else happens to sit inside
it.  These tests pin the two places that shows up as a wrong call.
"""
import sys
import pathlib

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import genotyping  # noqa: E402
import scorer  # noqa: E402


def _trace(peaks, sigma=1.0, n=3400):
    """One channel of ``acgt``: flat noise plus Gaussian peaks at given scans."""
    rng = np.random.default_rng(7)
    y = rng.normal(0.0, sigma, n)
    for scan, height in peaks:
        i = int(scan)
        y[i] += height * sigma
    return y


def _quartoctet(centres):
    """A standard quartet at *centres*, the way find_is_quartet reports it."""
    return [int(round(c)) for c in centres]


def _model_from(offsets, base=2600.0, **kw):
    """A model primed from wells whose product sits at *offsets* from H1."""
    model = genotyping.PlateISModel.from_marks([_quartoctet([base, base + 81, base + 145, base + 155])
                                                for _ in range(6)])
    wins = []
    for o in offsets:
        acgt = np.zeros((3400, 4))
        acgt[:, genotyping.acgt_index_for_channel("ACTG", 3)] = _trace(
            [(base, 900), (base + 81, 900), (base + 145, 900), (base + 155, 900)])
        # The sample's product lands at the offset the model is learning.
        acgt[:, genotyping.acgt_index_for_channel("ACTG", 2)] = _trace(
            [(base + o, 800)])
        wins.append((acgt, genotyping.acgt_index_for_channel("ACTG", 2),
                     _quartoctet([base, base + 81, base + 145, base + 155])))
    genotyping.prime_plate_co_migration(model, wins)
    return model


class TestCoMigrationWindow:
    def test_learns_the_plate_offset(self):
        model = _model_from([-12] * 8 + [-9, -15, -11, -13])
        win = model.co_migration_window(0)
        assert win is not None
        assert win[0] < -12 < win[1]

    def test_too_few_references_installs_no_gate(self):
        # Refusing to gate is the safe direction: the ungated search is exactly
        # the pre-existing behaviour, so a missing gate can never regress.
        model = _model_from([-12, -11, -10, -13, -9, -12, -11])
        assert model.co_migration_window(0) is None

    def test_minority_of_bad_references_cannot_widen_the_gate(self):
        # median/MAD over a mode is unmovable by a minority of junk, which
        # is the whole reason for the robust pair: ordinary noise in the H2
        # neighbourhood on a hom-1 plate must not redefine H1.
        model = _model_from([-12] * 8 + [45, -48, 42, 47, -45, 38])
        win = model.co_migration_window(0)
        assert win is not None
        assert win[0] <= -12 <= win[1]
        assert win[1] < 20


    def test_gate_wider_than_a_product_is_refused(self):
        # No mode at all -- a reference set that is just scattered.  A gate
        # wider than any real spread is a contaminated set, and installing it
        # would re-admit the off-target peak it exists to exclude.
        model = _model_from([-45, -35, -25, -15, -5, 5, 15, 25, 35, 45])
        assert model.co_migration_window(0) is None


class TestGateRejectsOffProduct:
    def _well(self, product_offset, decoy_offset, decoy_height):
        base = 2600.0
        acgt = np.zeros((3400, 4))
        acgt[:, 3] = _trace([(base, 900), (base + 81, 900),
                             (base + 145, 900), (base + 155, 900)])
        acgt[:, 1] = _trace([(base + product_offset, 800),
                             (base + decoy_offset, decoy_height)])
        return acgt

    def _row(self, acgt, model, base=2600.0):
        class _Doc:
            pass
        d = _Doc()
        d.acgt = acgt
        d.well = "T01"
        return genotyping.auto_genotype(
            d, is_channel=3, sample_channel=2, base_order="ACTG",
            std_scans_manual=_quartoctet([base, base + 81, base + 145, base + 155]),
            is_model=model)

    def test_peak_on_the_wrong_side_of_the_standard_is_not_product(self):
        """F12: only a peak at +12, so there is no co-migrating product."""
        model = _model_from([-12] * 8 + [-9, -15, -11, -13])
        # A decoy well past H2's centre, not anywhere near a product position.
        row = self._row(self._well(400, 12, 44), model)
        assert row["call"] == "no-call"

    def test_co_migrating_peak_is_still_called(self):
        """The gate must not throw away the peak it is meant to protect."""
        model = _model_from([-12] * 8 + [-9, -15, -11, -13])
        row = self._row(self._well(-12, 400, 44), model)
        assert row["call"] == "hom-1"


class TestHomRequiresItsOwnSignal:
    """A hom allele has to be there, not merely out-areaed by a rival.

    ``t9_call`` gates on ``max(sigmas)``, which proves only that *some* band is
    strong.  A heteroduplex will carry a homozygous call that its own
    homoduplex does not support, and an off-target peak integrated across a wide
    segment will out-area a real allele by accident of geometry.  Both of those
    are ABCC2 wells.
    """

    def test_strong_heteroduplex_cannot_carry_a_hom_call(self):
        # H12's shape: hom1 barely present, het1 four times over the floor.
        call, _frac, _flags = scorer.t9_call(925, 0, 46340, 0, (8.5, 2.4, 44.0, 9.2))
        assert call == "no-call"

    def test_off_target_area_cannot_outrank_a_real_allele(self):
        # H01's shape: a genuine hom1, and a second "allele" that is large in
        # area but only 16.9 sigma at the band it is being read from.
        call, _frac, _flags = scorer.t9_call(14332, 25282, 972, 189,
                                             (144.6, 16.9, 0.6, 3.8))
        assert call == "hom-1"

    def test_genuine_homozygote_still_calls(self):
        call, _frac, _flags = scorer.t9_call(231448, 0, 40700, 0,
                                             (1717.7, -3.3, 7.4, 16.9))
        assert call == "hom-1"

    def test_real_hom_2_still_calls(self):
        call, _frac, _flags = scorer.t9_call(1000, 200000, 500, 500,
                                             (5.0, 900.0, 2.0, 2.0))
        assert call == "hom-2"

    def test_area_alone_still_calls_when_no_significance_is_given(self):
        # A hand-built peak table carries areas and nothing else.  The new
        # requirement must not silently no-call a table that never had sigmas.
        call, _frac, _flags = scorer.t9_call(14332, 25282, 0, 0)
        assert call == "hom-2"

    def test_heterozygote_unaffected(self):
        call, _frac, flags = scorer.t9_call(190067, 39272, 10722, 3808,
                                            (2714.3, 185.7, 100.9, 107.1))
        assert call == "het"
        assert "ai" in flags