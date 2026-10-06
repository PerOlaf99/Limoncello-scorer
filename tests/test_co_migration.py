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


class TestModelBuildingPrimesTheGate:
    """The gate has to be learned where the model is actually built.

    ``prime_plate_co_migration`` is only useful if the production path calls
    it: ``from_wells`` hands the marked wells straight to it, so a model built
    the way the GUI builds one carries the calibration, while a caller that
    knows only the standard column (and cannot name the sample column) is left
    exactly as before -- ungated, fixed-radius search.
    """

    def _model(self, offsets, with_samp=True):
        base = 2600.0
        wells = []
        for o in offsets:
            acgt = np.zeros((3400, 4))
            acgt[:, 3] = _trace([(base, 900), (base + 81, 900),
                                 (base + 145, 900), (base + 155, 900)])
            acgt[:, 1] = _trace([(base + o, 800)])
            wells.append((acgt, 3, _quartoctet(
                [base, base + 81, base + 145, base + 155])))
        if with_samp:
            return genotyping.PlateISModel.from_wells(wells, samp_col=1)
        return genotyping.PlateISModel.from_wells(wells)

    def test_from_wells_primes_the_gate(self):
        model = self._model([-12] * 8 + [-9, -15, -11, -13])
        win = model.co_migration_window(0)
        assert win is not None
        assert win[0] < -12 < win[1]

    def test_omitting_samp_col_installs_no_gate(self):
        # The GUI path that knows only is_col keeps the old behaviour -- no
        # gate -- instead of silently gating to a column it does not own.
        model = self._model([-12] * 12, with_samp=False)
        assert model.co_migration_window(0) is None

    def test_built_model_rejects_the_off_product_peak(self):
        # The F12 shape end to end: the model built from from_wells, then the
        # well carrying only a decoy peak, must read no-call.
        model = self._model([-12] * 8 + [-9, -15, -11, -13])
        row = TestGateRejectsOffProduct()._row(
            TestGateRejectsOffProduct()._well(400, 12, 44), model)
        assert row["call"] == "no-call"

    def test_built_model_still_calls_the_co_migrating_peak(self):
        # And the same model must not throw away the real product.
        model = self._model([-12] * 8 + [-9, -15, -11, -13])
        row = TestGateRejectsOffProduct()._row(
            TestGateRejectsOffProduct()._well(-12, 400, 44), model)
        assert row["call"] == "hom-1"


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


class TestWellSelfConsistency:
    """The verdict reads a well from inside out, not through the plate gate.

    A gate is the plate's *average* drift, and a genuine run that shifted every
    delta together must not be read through an average that no longer applies
    (F05's CYBA het at +8/+8/+17).  So the well's dominant band anchors the
    call, the plate only referees the cases that cannot decide themselves.
    """

    def test_lone_band_off_the_plate_home_is_an_orphan(self):
        # F12: only one band at all, +12 on a plate whose products live at
        # -12.  Nothing else ran, so the plate is the only referee -- and the
        # peak contradicts it.
        model = _model_from([-12] * 8 + [-9, -15, -11, -13])
        cands = {0: (12.0, 585.0)}
        assert genotyping._co_migration_verdict(
            cands, [43.9, 3.3, 2.5, 2.7], model) == []

    def test_lone_band_on_the_plate_home_is_a_product(self):
        model = _model_from([-12] * 8 + [-9, -15, -11, -13])
        cands = {0: (-9.0, 585.0)}
        assert genotyping._co_migration_verdict(
            cands, [120.0, 0.1, 0.0, 0.0], model) == [0]

    def test_dominant_beats_an_internal_pair_that_contradicts_it(self):
        # E11: two residual peaks in the heteroduplex region agree with each
        # other (both ~-39) but agree with nothing else -- and the well's real
        # allele, its dominant band, sits at the plate's -13.  The pair must
        # not outvote the dominant.
        model = _model_from([-12] * 8 + [-9, -15, -11, -13])
        cands = {0: (36.0, 1209.0), 1: (-13.0, 56588.0),
                 2: (-45.0, 2312.0), 3: (-33.0, 1555.0)}
        assert genotyping._co_migration_verdict(
            cands, [19.2, 3785.9, 25.9, 24.4], model) == [1]

    def test_drifted_hom_dominant_kept_when_other_signal_present(self):
        # D07: a real T9 hom-2 whose product drifted to -33, 22 scans off the
        # plate's -11.  A second band is present, so the well genuinely ran,
        # and the plate's average no longer speaks for it: the dominant stands.
        model = _model_from([-11] * 8 + [-9, -13, -12, -10])
        cands = {1: (-33.0, 17638.0), 3: (16.0, 3702.0)}
        assert genotyping._co_migration_verdict(
            cands, [8.7, 223.0, -23.9, 23.2], model) == [1]

    def test_drifted_set_is_all_kept(self):
        # F05: every band shifted right together, so the whole set -- not the
        # plate average -- is the well talking.
        model = _model_from([-10] * 12)
        cands = {0: (21.0, 1429.0), 1: (8.0, 1534.0), 2: (17.0, 1574.0)}
        assert genotyping._co_migration_verdict(
            cands, [75.4, 116.5, 51.5], model) == [0, 1, 2]


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

class TestHetRequiresTheSecondAllele:
    """A heterozygote makes both homoduplexes, so one of them must be there.

    Taq A splits CYBA's product into two sub-peaks a near-constant distance
    apart.  The second one is a conformer of the same allele, not a second
    allele, but it sits in the heteroduplex neighbourhood and reads as one.
    Five wells were called ``het`` this way -- CYBA_N1 B01 plus CYBA_N2 59C
    A04, C05, D03 and F03 -- every one of them with the second homoduplex
    reading exactly zero area while a reference genotype says homozygous.

    Across every het call carrying a reference genotype, 0 of 31 genuine
    heterozygotes had an empty second homoduplex against 5 of 11 false ones, so
    requiring that area costs no true call.

    The floor on that second allele is significance, not area, and the second
    homoduplex is held to it too -- see
    test_a_second_allele_must_be_visible_not_merely_present.
    """

    def test_taq_conformer_is_not_a_second_allele(self):
        # CYBA_N2 59C A04: one strong homoduplex, a merged heteroduplex band
        # well above the floor, and no second homoduplex at all.
        call, _frac, _flags = scorer.t9_call(200546, 0, 46340, None,
                                             (1180.0, 0.0, 96.4))
        assert call == "hom-1"

    def test_split_heteroduplex_still_needs_the_second_homoduplex(self):
        # The separated form of the same trap: both heteroduplexes present and
        # comparable, but still only one allele's homoduplex.
        call, _frac, _flags = scorer.t9_call(145988, 0, 30011, 31002,
                                             (742.1, 0.0, 88.0, 91.4))
        assert call == "hom-1"

    def test_genuine_heterozygote_still_calls(self):
        # Both homoduplexes present, heteroduplexes merged: unchanged.
        call, _frac, _flags = scorer.t9_call(190067, 39272, 10722, None,
                                             (2714.3, 185.7, 100.9))
        assert call == "het"

    def test_requirement_is_presence_not_magnitude(self):
        # The second allele does not have to match the first in size.  A band
        # carrying a thousandth of the dominant's area is still a band, and a
        # het whose minor allele is faint is a normal het.
        call, _frac, _flags = scorer.t9_call(190067, 1200, 10722, None,
                                             (2714.3, 45.0, 100.9))
        assert call == "het"

    def test_a_second_allele_must_be_visible_not_merely_present(self):
        # RS1695_N2 A04, which called het on this: heteroduplex peaks clear of
        # the noise at 40 and 56 sigma, a dominant homoduplex at 6447, and a
        # first homoduplex at 7.5 sigma -- 0.7% of the real peak.  Real hets on
        # the same plate clear T9_MIN_DOMINANT_SIGMA on both homoduplexes
        # (tightest: G10 at 46.0), so this is what separates the two without
        # costing a true call.
        call, _frac, _flags = scorer.t9_call(611754, 6130, 0, 4465,
                                             (7.5, 6447.3, 40.3, 56.2))
        assert call == "hom-2"

    def test_a_minor_allele_pinned_to_the_window_edge_is_not_an_allele(self):
        # RS1695_N2 A01, which called het on this: the sample carries one peak,
        # in the hom2 window at 2334, and the only thing in the hom1 window is
        # that peak's own rising flank reaching back past IS band 1. So the
        # hom1 reading of 34 sigma is the same molecule counted twice.
        call, _frac, _flags = scorer.t9_call(390123, 76560, 779, 7461,
                                             (34.2, 7045.2, 89.1, 95.3))
        assert call == "hom-2"
