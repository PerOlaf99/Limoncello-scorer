"""Mass-action MF (CTCE): MF = (A_MUT + 1/2 A_HET) / (A_WT + A_MUT + A_HET).

Pure-unit checks against hand-computed areas, plus the real OY A01 well to
prove the wiring on actual picks.
"""
import unittest
from pathlib import Path

import pytest

pytest.importorskip("tkinter")
pytest.importorskip("matplotlib")

import genotyping as g  # noqa: E402

A01 = Path("/media/per/78B0C7DE1FA7081C/OY/OY_rs1695_T9_270910Run01/A01.rsd")


def mkrec(scan, col, kind, area, base="A", gid=None):
    return {"file": "synthetic.rsd", "well": "A01", "scan": scan,
            "col": col, "channel": f"Channel{col + 1}", "base": base,
            "kind": kind, "area": area, "height": 1.0, "color": "k",
            "left": scan - 3, "right": scan + 3, "gid": gid}


def picker(mains, std=None):
    pk = g.PeakPicker.__new__(g.PeakPicker)
    pk.records = list(mains)
    pk.std = std
    pk.length_bp = None
    pk._d2_cache = {}
    pk.doc = type("D", (), {"well": "A01", "n_scans": 5000, "acgt": None})()
    return pk


class TestDuplexTagging(unittest.TestCase):
    def test_four_peaks_tag_two_and_two(self):
        """The internal-standard pattern: 2 homoduplexes + 2 heteroduplexes."""
        mains = [mkrec(100 + 10 * i, 2, "main", 1000.0, gid=i) for i in range(4)]
        pk = picker(mains)
        msg = pk.mark_duplex()
        self.assertEqual([n for _, n in pk.std],
                         ["HOM1", "HOM2", "HET1", "HET2"])
        self.assertIn("Duplex species", msg)

    def test_separate_positions_do_not_merge(self):
        """Two het positions far apart stay two duplex sets: the tagged one is
        the last-picked main's own position only."""
        mains = [mkrec(100 + 10 * i, 2, "main", 1000.0, gid=i)
                 for i in range(2)]
        mains += [mkrec(500 + 10 * i, 2, "main", 900.0, gid=2 + i)
                  for i in range(2)]
        pk = picker(mains)
        pk.mark_duplex()
        self.assertEqual([s for s, _ in pk.std], [500, 510])
        self.assertEqual([n for _, n in pk.std], ["HOM1", "HOM2"])
        self.assertIsNone(pk.mass_action(mains[0]))
        self.assertIsNotNone(pk.mass_action(mains[2]))

    def test_three_peaks_are_one_homo_plus_two_hets(self):
        """A homozygote's re-annealed strands make 1+2, not a truncated 2+1
        standard -- naming it 2+1 would call a heteroduplex a homoduplex and
        halve a low MF."""
        mains = [mkrec(100 + 10 * i, 2, "main", 1000.0, gid=i) for i in range(3)]
        pk = picker(mains)
        pk.mark_duplex()
        self.assertEqual([n for _, n in pk.std], ["HOM1", "HET1", "HET2"])

    def test_two_peaks_tag_hom1_hom2_only(self):
        mains = [mkrec(100, 2, "main", 900.0, gid=0),
                 mkrec(110, 2, "main", 100.0, gid=1)]
        pk = picker(mains)
        pk.mark_duplex()
        self.assertEqual([n for _, n in pk.std], ["HOM1", "HOM2"])

    def test_untagged_channel_is_not_picked(self):
        """The duplex set comes from one channel only -- the last-picked one."""
        mains = [mkrec(100, 2, "main", 900.0, gid=0),
                 mkrec(110, 2, "main", 100.0, gid=1),
                 mkrec(500, 3, "main", 700.0, gid=2),
                 mkrec(510, 3, "main", 300.0, gid=3)]
        pk = picker(mains)
        pk.mark_duplex()
        self.assertEqual({m["col"] for m in pk.labelled_species()}, {3})

    def test_needs_two_mains(self):
        pk = picker([mkrec(100, 2, "main", 900.0, gid=0)])
        with self.assertRaises(ValueError):
            pk.mark_duplex()
        pk = picker([mkrec(100, 2, "plus_a", 5.0, gid=0)])
        with self.assertRaises(ValueError):
            pk.mark_duplex()
        # one main + one +A is not a duplex position
        pk = picker([mkrec(100, 2, "main", 900.0, gid=0),
                     mkrec(110, 2, "plus_a", 5.0, gid=1)])
        with self.assertRaises(ValueError):
            pk.mark_duplex()

    def test_extra_peaks_beyond_four_are_ignored(self):
        mains = [mkrec(100 + 10 * i, 2, "main", 1000.0, gid=i) for i in range(6)]
        pk = picker(mains)
        pk.mark_duplex()
        self.assertEqual(len(pk.std), 4)
        self.assertEqual([n for _, n in pk.std],
                         ["HOM1", "HOM2", "HET1", "HET2"])

    def test_mark_std_sets_length(self):
        mains = [mkrec(100, 2, "main", 900.0, gid=0),
                 mkrec(110, 2, "main", 100.0, gid=1)]
        pk = picker(mains)
        pk.mark_std(216.0)
        self.assertEqual(pk.length_bp, 216.0)
        self.assertIn("216", pk.mark_std(216.0))

    def test_clear_std(self):
        mains = [mkrec(100, 2, "main", 900.0, gid=0),
                 mkrec(110, 2, "main", 100.0, gid=1)]
        pk = picker(mains)
        pk.mark_std(216.0)
        pk.clear_std()
        self.assertIsNone(pk.std)
        self.assertIsNone(pk.length_bp)


class TestMassAction(unittest.TestCase):
    def mf_for(self, areas):
        """Build a duplex position and return mass_action.  `areas` is the
        homoduplex area(s) then the heteroduplex area(s), in migration order."""
        scans = [100, 112, 200, 212]
        mains, std = [], []
        for i, a in enumerate(areas):
            if a is None:
                continue
            mains.append(mkrec(scans[i], 2, "main", float(a), gid=i))
            std.append((scans[i], ["HOM1", "HOM2", "HET1", "HET2"][i]))
        pk = picker(mains, std=std)
        return pk.mass_action(mains[0])

    def test_homozygote_wt_reads_zero(self):
        """One homoduplex, no heteroduplex: homozygous wild type, MF 0."""
        ma = self.mf_for([1000, None, None, None])
        self.assertAlmostEqual(ma["mf"], 0.0)
        self.assertIsNone(ma["ai"])
        self.assertEqual(ma["n_homoduplex"], 1)

    def test_het_reads_half_not_quarter(self):
        """A true heterozygote: two equal homoduplexes plus re-annealed
        heteroduplex.  MF 0.5, not the 0.25 a plain area ratio of the two
        homoduplexes gives -- that factor of two is the whole point."""
        ma = self.mf_for([1000, 1000, 1000, 1000])
        self.assertAlmostEqual(ma["mf"], (1000 + 0.5 * 2000) / 4000)
        self.assertAlmostEqual(ma["mf"], 0.5)
        # without the 1/2-HET term the same well would read 0.25
        self.assertAlmostEqual(1000 / 4000, 0.25)

    def test_het_skewed_by_reannealing(self):
        """Unequal heteroduplex areas still read 0.5: the 1/2 term is applied
        to their sum, and a 50/50 het has A_HET = A_WT."""
        ma = self.mf_for([1000, 1000, 1500, 500])
        self.assertAlmostEqual(ma["mf"], (1000 + 0.5 * 2000) / 4000)
        self.assertAlmostEqual(ma["mf"], 0.5)

    def test_low_mf_no_mutant_homoduplex(self):
        """Below ~5 % MF essentially all mutant strands are in heteroduplex, so
        the well is 1 homoduplex + 2 heteroduplex and MF = 1/2 A_HET /
        (A_WT + A_HET)."""
        ma = self.mf_for([1000, None, 100, 100])
        self.assertAlmostEqual(ma["mf"], (0 + 0.5 * 200) / 1200)
        self.assertAlmostEqual(ma["mf"], 0.08333333333333333, places=9)
        self.assertEqual(ma["n_homoduplex"], 1)
        self.assertIsNone(ma["ai"])

    def test_low_mf_with_faint_mutant_homoduplex(self):
        ma = self.mf_for([1000, 20, 100, 100])
        self.assertAlmostEqual(ma["mf"], (20 + 0.5 * 200) / 1220)
        self.assertAlmostEqual(ma["mf"], 0.098360655737704915, places=9)
        self.assertEqual(ma["n_homoduplex"], 2)

    def test_larger_homoduplex_is_wild_type(self):
        ma = self.mf_for([1000, 100, None, None])
        self.assertEqual(ma["a_wt"], 1000.0)
        self.assertEqual(ma["a_mut"], 100.0)
        self.assertAlmostEqual(ma["mf"], 100 / 1100)

    def test_ai_uses_label_order_not_size(self):
        """AI = A_HOMO1/(A_HOMO1+A_HOMO2) -- unaffected by which is bigger."""
        ma = self.mf_for([100, 1000, None, None])
        self.assertAlmostEqual(ma["ai"], 100 / 1100)
        self.assertEqual(ma["a_wt"], 1000.0)
        self.assertEqual(ma["a_mut"], 100.0)

    def test_mf_needs_a_homoduplex(self):
        """Heteroduplexes alone carry no MF: the denominator would be all
        half-counted product."""
        mains = [mkrec(100, 2, "main", 1000.0, gid=0),
                 mkrec(200, 2, "main", 100.0, gid=1)]
        pk = picker(mains, std=[(100, "HET1"), (200, "HET2")])
        self.assertIsNone(pk.mass_action(mains[0]))

    def test_untagged_position_has_no_mf(self):
        mains = [mkrec(100, 2, "main", 1000.0, gid=0),
                 mkrec(112, 2, "main", 100.0, gid=1)]
        pk = picker(mains, std=None)
        self.assertIsNone(pk.mass_action(mains[0]))
        # the plain fraction still works
        self.assertAlmostEqual(pk.clust_frac(mains[0]), 100 / 1100, places=4)

    def test_zero_area_gives_none(self):
        mains = [mkrec(100, 2, "main", 0.0, gid=0),
                 mkrec(112, 2, "main", 0.0, gid=1)]
        pk = picker(mains, std=[(100, "HOM1"), (112, "HOM2")])
        self.assertIsNone(pk.mass_action(mains[0]))

    def test_plus_a_peak_has_no_mf(self):
        mains = [mkrec(100, 2, "main", 1000.0, gid=0),
                 mkrec(112, 2, "main", 100.0, gid=1),
                 mkrec(120, 2, "plus_a", 5.0, gid=2)]
        pk = picker(mains, std=[(100, "HOM1"), (112, "HOM2")])
        self.assertIsNone(pk.mass_action(mains[2]))


class TestExportAndTable(unittest.TestCase):
    def test_export_carries_mf(self):
        mains = [mkrec(100, 2, "main", 1000.0, gid=0),
                 mkrec(112, 2, "main", 100.0, gid=1),
                 mkrec(200, 2, "main", 400.0, gid=2),
                 mkrec(212, 2, "main", 400.0, gid=3)]
        pk = picker(mains, std=[(100, "HOM1"), (112, "HOM2"),
                                (200, "HET1"), (212, "HET2")])
        rows = pk.export_rows()
        self.assertEqual(len(rows), 4)
        for row in rows:
            self.assertIn("mf", row)
            self.assertIn("ai", row)
        self.assertAlmostEqual(rows[0]["mf"],
                               (100 + 0.5 * 800) / 1900, places=4)
        self.assertEqual(rows[0]["fraction"], round(100 / 1100, 4))

    def test_export_blank_mf_when_untagged(self):
        mains = [mkrec(100, 2, "main", 1000.0, gid=0),
                 mkrec(112, 2, "main", 100.0, gid=1)]
        pk = picker(mains, std=None)
        rows = pk.export_rows()
        self.assertEqual(rows[0]["mf"], "")
        self.assertEqual(rows[0]["ai"], "")

    def test_export_sorted_by_scan(self):
        mains = [mkrec(200, 2, "main", 400.0, gid=0),
                 mkrec(100, 2, "main", 1000.0, gid=1),
                 mkrec(112, 2, "main", 100.0, gid=2)]
        pk = picker(mains)
        self.assertEqual([r["scan"] for r in pk.export_rows()],
                         [100, 112, 200])


class TestRealA01(unittest.TestCase):
    """MF on the real OY rs1695 A01 well: two heterozygous allele pairs on the
    sample channel.  Tagged 2+0, MF must equal the plain area fraction there."""

    @classmethod
    def setUpClass(cls):
        if not A01.exists():
            raise unittest.SkipTest("A01.rsd not present")
        from analyzer_core import load_trace
        doc = load_trace(A01, base_order="ACTG")
        pk = g.PeakPicker(doc, A01, include_sh=False)
        cls.pk = pk
        cls.pairs = []
        for scan, vol in ((2123, 19326), (2143, 3125),
                          (2205, 20134), (2225, 3121)):
            rec = pk.pick(scan, vol=vol)
            assert rec is not None, f"nothing picked at {scan}"
            assert rec["col"] == 1, f"{scan}: col {rec['col']}, want 1 (Ch2)"
            cls.pairs.append(rec)

    def test_the_two_positions_stay_apart(self):
        """2123/2143 and 2205/2225 are two separate heterozygous positions ~100
        scans apart, so the duplex set covers one of them, not the whole
        channel."""
        self.pk.mark_duplex()
        self.assertEqual([n for _, n in self.pk.std], ["HOM1", "HOM2"])
        self.assertEqual([s for s, _ in self.pk.std], [2205, 2225])

    def test_mf_of_the_tagged_position(self):
        ma = self.pk.mass_action(self.pairs[3])
        self.assertIsNotNone(ma)
        self.assertAlmostEqual(ma["a_het"], 0.0, places=6)
        self.assertAlmostEqual(ma["mf"], ma["a_mut"] / (ma["a_wt"] + ma["a_mut"]))
        self.assertGreater(ma["mf"], 0.10)
        self.assertLess(ma["mf"], 0.16)

    def test_untagged_position_has_no_mf(self):
        self.pk.mark_duplex()
        self.assertIsNone(self.pk.mass_action(self.pairs[0]))

    def test_export_mf_matches(self):
        self.pk.mark_duplex()
        row = next(r for r in self.pk.export_rows() if r["scan"] == 2205)
        self.assertNotEqual(row["mf"], "")
        self.assertAlmostEqual(row["mf"], row["fraction"], places=3)
        self.assertAlmostEqual(row["ai"], 1 - row["fraction"], places=3)
        blank = next(r for r in self.pk.export_rows() if r["scan"] == 2123)
        self.assertEqual(blank["mf"], "")


if __name__ == "__main__":
    unittest.main(verbosity=2)
