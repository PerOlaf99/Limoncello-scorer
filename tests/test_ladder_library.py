"""The ladder library: user ladders on disk, validation, import, exclusions."""
import json
from pathlib import Path

import numpy as np
import pytest

import fragment_sizing as fs
import ladder_library as lib
import mark_library
from fragment_sizing import Ladder


@pytest.fixture(autouse=True)
def _private_library(monkeypatch, tmp_path):
    monkeypatch.setattr(mark_library, "LIBRARY_DIR", tmp_path / "library")


def _lad(name="My ladder", lengths=(50, 100, 150, 200, 250), **kw):
    return Ladder(name=name, lengths=tuple(float(b) for b in lengths), **kw)


# --------------------------------------------------------------------------- #
# storage
# --------------------------------------------------------------------------- #
def test_save_find_list_and_delete_a_user_ladder():
    key = lib.save_user_ladder(_lad(dye="FAM", color="#1f77b4", default_channel=2))
    assert key == "my_ladder"
    assert (lib.ladder_dir() / "my_ladder.json").is_file()
    got = lib.find("My Ladder")
    assert got.lengths == (50, 100, 150, 200, 250) and got.dye == "FAM"
    assert got.default_channel == 2 and got.color == "#1f77b4"
    assert "my_ladder" in fs.list_ladders()
    assert fs.load_ladder("my_ladder").name == "My ladder"
    kinds = {r["key"]: r["builtin"] for r in lib.all_ladders()}
    assert kinds["my_ladder"] is False and kinds["et400r"] is True
    lib.delete_user_ladder(key)
    assert lib.find("My ladder") is None


def test_saving_twice_needs_overwrite():
    lib.save_user_ladder(_lad())
    with pytest.raises(FileExistsError):
        lib.save_user_ladder(_lad())
    lib.save_user_ladder(_lad(lengths=(60, 90, 120, 150)), key="my_ladder",
                         overwrite=True)
    assert lib.find("my_ladder").lengths == (60, 90, 120, 150)


def test_bundled_ladders_cannot_be_replaced_or_deleted():
    with pytest.raises(ValueError):
        lib.save_user_ladder(_lad(name="ET400R"))
    with pytest.raises(ValueError):
        lib.delete_user_ladder("et400r")


def test_unique_name_avoids_clashes():
    lib.save_user_ladder(_lad())
    assert lib.unique_name("My ladder") == "My ladder (2)"
    assert lib.unique_name("Other") == "Other"


def test_a_broken_file_does_not_break_the_library():
    lib.ladder_dir().mkdir(parents=True)
    (lib.ladder_dir() / "bad.json").write_text("{not json", encoding="utf-8")
    lib.save_user_ladder(_lad())
    assert set(lib.user_ladders()) == {"my_ladder"}


# --------------------------------------------------------------------------- #
# validation
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("lad,expect", [
    (_lad(lengths=(50, 100, 150)), "at least 4"),
    (_lad(lengths=(50, 100, 100, 150, 200)), "duplicate"),
    (_lad(exclude=(75,)), "not in the ladder"),
    (_lad(lengths=(50, 100, 150, 200), exclude=(100,)), "fewer than 4"),
    (_lad(default_channel=5), "1-4"),
    (_lad(name=" "), "name"),
])
def test_validation_names_the_problem(lad, expect):
    problems = lib.validate(lad)
    assert any(expect in p for p in problems), problems
    with pytest.raises(ValueError):
        lib.save_user_ladder(lad)


# --------------------------------------------------------------------------- #
# import
# --------------------------------------------------------------------------- #
def test_import_json_one_and_many(tmp_path):
    one = tmp_path / "one.json"
    one.write_text(json.dumps({"name": "A", "dye": "LIZ", "lengths": [20, 40, 60, 80],
                               "exclude": [40]}), encoding="utf-8")
    [a] = lib.import_ladder_file(one)
    assert a.name == "A" and a.dye == "LIZ" and a.exclude == (40.0,)
    many = tmp_path / "many.json"
    many.write_text(json.dumps({"ladders": [{"name": "B", "lengths": [1, 2, 3, 4]},
                                            {"name": "C", "lengths": [5, 6, 7, 8]}]}),
                    encoding="utf-8")
    assert [x.name for x in lib.import_ladder_file(many)] == ["B", "C"]


@pytest.mark.parametrize("text,expect", [
    ("60,90,100,120,150", [60, 90, 100, 120, 150]),                 # one row
    ("60\n90\n100\n120\n", [60, 90, 100, 120]),                      # one column
    ("Size (bp),Height\n35,10\n50,12\n75,9\n100,11\n", [35, 50, 75, 100]),
    ("height;size\n10;35,5\n12;50\n9;75\n11;100\n", [35.5, 50, 75, 100]),  # EU csv
    ("1\nMapping_Set_II_V1_1\nET-ROX\tFAM\tTMR\tHEX\n"                # GP .tdb
     "60\t90\t100\t120\t150\t160\nA02\tA04\tA06\tA08\n", [60, 90, 100, 120, 150, 160]),
])
def test_import_text_formats(tmp_path, text, expect):
    path = tmp_path / "ladder.txt"
    path.write_text(text, encoding="utf-8")
    [lad] = lib.import_ladder_file(path)
    assert list(lad.lengths) == expect
    assert lad.name == "ladder" and "imported" in lad.source


def test_import_without_lengths_says_so(tmp_path):
    path = tmp_path / "empty.txt"
    path.write_text("no numbers here\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no fragment lengths"):
        lib.import_ladder_file(path)


def test_export_round_trips(tmp_path):
    lad = _lad(vendor="Acme", dye="FAM", exclude=(150,), default_channel=3,
               notes="lot 7")
    lib.export_ladder(lad, tmp_path / "x.json")
    [back] = lib.import_ladder_file(tmp_path / "x.json")
    assert back.lengths == lad.lengths and back.exclude == lad.exclude
    assert back.vendor == "Acme" and back.default_channel == 3 and back.notes == "lot 7"


# --------------------------------------------------------------------------- #
# excluded bands in sizing
# --------------------------------------------------------------------------- #
class _Doc:
    def __init__(self, acgt):
        self.acgt, self.well, self.path, self.base_order = acgt, "A01", "A01.rsd", "ACGT"


def test_excluded_bands_are_matched_but_do_not_anchor_the_curve():
    lengths = [60, 90, 100, 120, 150, 160, 170, 190, 200, 220]
    scan = lambda bp: 1500 + 9.0 * bp                       # noqa: E731
    x = np.arange(5000)
    acgt = 5.0 + np.random.default_rng(0).normal(0, 2.0, (5000, 4))
    for b in lengths:
        acgt[:, 0] += 900 * np.exp(-0.5 * ((x - scan(b)) / 4.0) ** 2)
    for b in (110, 175):
        acgt[:, 1] += 700 * np.exp(-0.5 * ((x - scan(b)) / 4.0) ** 2)
    lad = Ladder(name="t", lengths=tuple(map(float, lengths)), exclude=(150.0,))
    res = fs.size_trace(_Doc(acgt), lad, ladder_channel=1, sample_channel=2,
                        base_order="ACGT")
    assert 150.0 not in [b for _, b in res.anchors]
    assert len(res.anchors) == len(lengths) - 1
    got = sorted(r["length_bp"] for r in res.rows if r["height"] > 300)
    assert got == pytest.approx([110, 175], abs=0.3)


# --------------------------------------------------------------------------- #
# bundled tables, pinned to their sources (band count, first, last, excluded)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("key,n,first,last,exclude", [
    ("et400r", 20, 60, 400, []),            # GP 2.2 test_plate.tdb / gpcom.dll
    ("et550r", 22, 60, 550, []),            # GP 2.2 gpcom.dll
    ("et900r", 37, 60, 900, []),            # GP 2.2 gpcom.dll (loop count 0x25)
    ("genescan120_liz", 9, 15, 120, []),    # Thermo Pub. 4474504, ch. 2
    ("genescan350_rox", 12, 35, 350, [250]),
    ("genescan400hd_rox", 21, 50, 400, []),
    ("genescan500_rox", 16, 35, 500, [250]),
    ("genescan500_liz", 16, 35, 500, [250]),
    ("genescan600_liz", 36, 20, 600, []),
    ("genescan1200_liz", 68, 20, 1200, []),
    ("promega_ils600", 22, 60, 600, []),    # Promega TM255
    ("promega_ils500", 21, 60, 500, []),    # Promega TMD039 + WEN/CC5 letter
    ("beckman_ss400", 22, 60, 420, []),     # Beckman PN A32350
    ("beckman_ss600", 33, 60, 640, []),
    ("geneflo1000_rox", 25, 400, 1000, []),
])
def test_bundled_tables_match_their_sources(key, n, first, last, exclude):
    lad = fs.load_ladder(key)
    assert len(lad.lengths) == n
    assert (lad.lengths[0], lad.lengths[-1]) == (first, last)
    assert list(lad.exclude) == exclude
    assert lad.source and not lib.validate(lad)


def test_every_bundled_ladder_is_pinned():
    pinned = {"et400r", "et550r", "et900r", "genescan120_liz", "genescan350_rox",
              "genescan400hd_rox", "genescan500_rox", "genescan500_liz",
              "genescan600_liz", "genescan1200_liz", "promega_ils600",
              "promega_ils500", "beckman_ss400", "beckman_ss600",
              "geneflo1000_rox"}
    assert set(fs.BUILTIN_LADDERS) == pinned
