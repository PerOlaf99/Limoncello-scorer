"""The auto-genotyping GUI: the menu item, the channel-role check, the results
table and the export.  Skips where there is no display to put a window on."""
import csv
import os
from pathlib import Path

import pytest

pytest.importorskip("tkinter")
pytest.importorskip("matplotlib")

if not os.environ.get("DISPLAY"):
    # Module-level skip needs allow_module_level, otherwise recent pytest treats
    # it as a collection error and aborts the whole run instead of skipping.
    pytest.skip("no DISPLAY for the GUI test", allow_module_level=True)

import tkinter as tk  # noqa: E402
from tkinter import filedialog, messagebox, simpledialog, ttk  # noqa: E402

import analyzer_core  # noqa: E402
import genotyping  # noqa: E402
import scorer  # noqa: E402
from sequence_analyzer import LimoncelloAnalyzerApp  # noqa: E402

RUN = Path("/media/per/78B0C7DE1FA7081C/OY/OY_rs1695_T9_270910Run01")
WELLS = ["A01", "A02", "A03", "A04"]


@pytest.fixture()
def app(monkeypatch, tmp_path):
    paths = [RUN / f"{w}.rsd" for w in WELLS]
    if not all(p.exists() for p in paths):
        pytest.skip(f"{RUN} not present")
    for name in ("showinfo", "showwarning", "showerror"):
        monkeypatch.setattr(messagebox, name, lambda *a, **k: "ok")
    monkeypatch.setattr(filedialog, "asksaveasfilename", lambda *a, **k: "")
    # the durable library lives in the app folder: keep it out of the tests,
    # otherwise one test's marks are restored into the next one's pickers
    import mark_library
    monkeypatch.setattr(mark_library, "LIBRARY_DIR", Path(tmp_path) / "library")
    try:
        app = LimoncelloAnalyzerApp()
    except tk.TclError as e:
        pytest.skip(f"cannot open a window: {e}")
    app.update()
    for p in paths:
        app.docs[str(p.resolve())] = analyzer_core.load_trace(
            p, base_order=app.base_order_var.get())
        app.file_list.insert(tk.END, f"run/{p.name}")
    app.n_graphs.set(1)
    app.selected = paths
    app.files = paths          # the sample window, for the no-selection fallback
    yield app
    app.destroy()


def _entries(menu):
    """[(type, label)] of a menu, skipping separators."""
    out = []
    for i in range(menu.index("end") + 1):
        try:
            out.append((menu.type(i), menu.entrycget(i, "label")))
        except tk.TclError:
            pass
    return out


def _find_genotyping_menu(app):
    """The Genotyping menu: a Tk child menu that owns the manual-pick item."""
    for w in app.winfo_children():
        if not isinstance(w, tk.Menu):
            continue
        if any(lbl == "Manual peak picking…" for _, lbl in _entries(w)):
            return w
    raise AssertionError("no Genotyping menu on the window")


def test_the_menu_actually_offers_auto_genotyping(app):
    """The whole point of adopting V3: the feature has to be reachable from
    the Genotyping menu, not only importable from Python."""
    gen = _find_genotyping_menu(app)
    entries = _entries(gen)
    assert ("cascade", "Auto-genotyping") in entries, entries
    i = gen.index("end") + 1
    for k in range(gen.index("end") + 1):
        if gen.type(k) == "cascade" and gen.entrycget(k, "label") == "Auto-genotyping":
            i = k
    sub = gen.nametowidget(gen.entrycget(i, "menu"))
    sub_labels = [lbl for typ, lbl in _entries(sub)]
    assert "Auto-genotype sample window…" in sub_labels, sub_labels
    assert "Channel roles (standard / sample)…" in sub_labels, sub_labels
    assert "Save auto-genotype table…" in sub_labels, sub_labels


def test_mark_all_refuses_without_a_learned_shape(app, monkeypatch):
    """96 unmarked wells is not a licence to guess each one separately.

    The batch command used to fall back to bare equimolar detection, which is
    what put noise in a third of the ABCC2 wells. With nothing learned and
    nothing picked it now says so and points at the marks it needs.
    """
    shown = []
    monkeypatch.setattr(messagebox, "showinfo",
                        lambda title, msg="", **k: shown.append(msg))
    monkeypatch.setattr(simpledialog, "askinteger",
                        lambda *a, **k: genotyping.DEFAULT_IS_CHANNEL)
    app.genotyping_active = True
    app._auto_is_model = None
    app._gen_mark_std_batch()
    assert shown and "Mark a few wells by hand first" in shown[0], shown
    assert "Auto-genotype" in shown[0], shown
    for well in WELLS:
        assert not app._ensure_picker(RUN / f"{well}.rsd").std, well


def test_mark_all_uses_the_learned_shape_and_tags_its_marks(app,
                                                            monkeypatch):
    """With a model, the batch places the standard from it -- and says so."""
    app._auto_is_channel.set(genotyping.DEFAULT_IS_CHANNEL)
    app._auto_sample_channel.set(genotyping.DEFAULT_SAMPLE_CHANNEL)
    seed = {"A01": [2135, 2214, 2425, 2520], "A02": [2172, 2266, 2477, 2563]}
    for well, scans in seed.items():
        pk = app._ensure_picker(RUN / f"{well}.rsd")
        pk.std = [(x, "") for x in scans]
    app.semi_auto_is()
    assert app._auto_is_model is not None

    monkeypatch.setattr(simpledialog, "askinteger",
                        lambda *a, **k: genotyping.DEFAULT_IS_CHANNEL)
    monkeypatch.setattr(messagebox, "showinfo", lambda *a, **k: "ok")
    monkeypatch.setattr(messagebox, "showwarning", lambda *a, **k: "ok")
    app._gen_mark_std_batch()

    # Seed wells keep the operator's own marks; the rest are placed and, unlike
    # a hand mark, are tagged as machine-found.
    for well in WELLS:
        pk = app._ensure_picker(RUN / f"{well}.rsd")
        assert pk.std, well
        if well not in seed:
            assert pk.std_auto, well


def test_semi_auto_marks_are_drawn_like_hand_picks(app):
    """The placed standards become visible main picks, not hidden attributes.

    ``semi_auto_is`` wrote ``pk.std`` but no records, and the plot and pick
    table only render ``records`` -- so the standard propagated to the other
    wells stayed invisible ("only manual peak picking marks anything") until
    it mattered at genotype time.
    """
    app._auto_is_channel.set(genotyping.DEFAULT_IS_CHANNEL)
    app._auto_sample_channel.set(genotyping.DEFAULT_SAMPLE_CHANNEL)
    seed = {"A01": [2135, 2214, 2425, 2520], "A02": [2172, 2266, 2477, 2563]}
    for well, scans in seed.items():
        app._ensure_picker(RUN / f"{well}.rsd").std = [(x, "") for x in scans]
    app.semi_auto_is()

    for well in WELLS:
        pk = app._ensure_picker(RUN / f"{well}.rsd")
        mains = [r for r in pk.records if r["kind"] == "main"]
        if well in seed:
            # the operator set pk.std directly here, so no record was expected
            continue
        assert mains, well
        for r in mains:
            assert pk.duplex_of(r), (well, r["scan"])
    # idempotent: a second semi-auto pass must not stack duplicate marks
    first = {w: len(app._ensure_picker(RUN / f"{w}.rsd").records)
             for w in WELLS}
    app.semi_auto_is()
    for well in WELLS:
        assert len(app._ensure_picker(RUN / f"{well}.rsd").records) == first[well]


def test_same_channel_for_standard_and_sample_is_refused(app):
    """Swapping the two roles does not fail loudly -- it scores the sample's
    own peaks as the standard and returns confident nonsense. So the GUI has
    to refuse the pair before a single well is scored."""
    app._auto_is_channel.set(3)
    app._auto_sample_channel.set(3)
    errors = []
    app.update()
    app.auto_genotype_wells()
    # the guard returns before scoring, so nothing was written
    assert app._auto_rows == []


def test_autogenotype_produces_a_row_per_well(app):
    app._auto_is_channel.set(genotyping.DEFAULT_IS_CHANNEL)
    app._auto_sample_channel.set(genotyping.DEFAULT_SAMPLE_CHANNEL)
    app.auto_genotype_wells()
    app.update()
    assert len(app._auto_rows) == len(WELLS)
    for r in app._auto_rows:
        assert r["call"] in scorer.CALLS, r["call"]
        assert r["well"]
        if r["call"] == "no-call":
            assert r["reason"], r


def test_results_table_shows_the_call_and_the_reason(app):
    app._auto_is_channel.set(genotyping.DEFAULT_IS_CHANNEL)
    app._auto_sample_channel.set(genotyping.DEFAULT_SAMPLE_CHANNEL)
    app.auto_genotype_wells()
    app.update()
    cols = list(app._auto_tree["columns"])
    assert "call" in cols and "reason" in cols and "frac" in cols
    items = app._auto_tree.get_children()
    assert len(items) == len(WELLS)
    first = app._auto_tree.item(items[0])["values"]
    assert first[cols.index("call")]


def test_export_writes_the_auto_rows(app, tmp_path, monkeypatch):
    app._auto_is_channel.set(genotyping.DEFAULT_IS_CHANNEL)
    app._auto_sample_channel.set(genotyping.DEFAULT_SAMPLE_CHANNEL)
    app.auto_genotype_wells()
    out = tmp_path / "auto.csv"
    monkeypatch.setattr(genotyping.filedialog, "asksaveasfilename",
                        lambda *a, **k: str(out))
    app.auto_genotype_save()
    with open(out, newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == len(WELLS)
    for key in ("well", "call", "frac", "reason", "is_channel",
                "sample_channel"):
        assert key in rows[0], list(rows[0])


def test_without_learned_geometry_no_standard_is_invented(app):
    """With no marks and no model there is no standard to find.

    The old engine had T9's geometry baked in, so a bare call to
    auto-genotype placed a standard on any plate that happened to share it.
    That is the behaviour being removed: without something to learn from, the
    wells must come back no-call rather than resting on another assay's
    numbers.
    """
    app._auto_is_channel.set(genotyping.DEFAULT_IS_CHANNEL)
    app._auto_sample_channel.set(genotyping.DEFAULT_SAMPLE_CHANNEL)
    app.auto_genotype_wells()
    got = {r["well"]: r for r in app._auto_rows}
    assert all(r["std_scans"] == "" for r in got.values()), got
    assert any(r["call"] == "no-call" for r in got.values()), got


def test_semi_automatic_standard_then_calls_match_ground_truth(app):
    """The real workflow on the real T9 wells: two operator marks, learn,
    place the rest, call -- and agree with tests/data/rs1695_expected.csv,
    the manual ground truth the 95/96 validation was measured against.
    Call names differ between the fixture ('hom2') and the engine ('hom-2')."""
    expected = {}
    with open(Path(__file__).parent / "data" / "rs1695_expected.csv",
              newline="") as f:
        for row in csv.DictReader(f):
            expected[row["well"]] = row["expected"]

    app._auto_is_channel.set(genotyping.DEFAULT_IS_CHANNEL)
    app._auto_sample_channel.set(genotyping.DEFAULT_SAMPLE_CHANNEL)

    # Two wells marked by hand, taken from the operator's IS.csv for this run.
    seed = {"A01": [2135, 2214, 2425, 2520], "A02": [2172, 2266, 2477, 2563]}
    for well, scans in seed.items():
        pk = app._ensure_picker(RUN / f"{well}.rsd")
        pk.std = [(x, "") for x in scans]

    app.semi_auto_is()
    assert app._auto_is_model is not None
    app.auto_genotype_wells()

    got = {r["well"]: r["call"].replace("hom-", "hom") for r in app._auto_rows}
    want = {w: expected.get(w) for w in WELLS}
    assert got == want, (got, want)

    # The wells the operator did not mark were placed, and are editable marks.
    placed = {r["well"]: r for r in app._auto_rows}
    for well in ("A03", "A04"):
        assert placed[well]["std_scans"], placed[well]
    assert "strong" in app._auto_is_model_note or "learned" in \
        app._auto_is_model_note, app._auto_is_model_note


def _run_auto(app):
    app._auto_is_channel.set(genotyping.DEFAULT_IS_CHANNEL)
    app._auto_sample_channel.set(genotyping.DEFAULT_SAMPLE_CHANNEL)
    app.auto_genotype_wells()
    app.update()


def _mark_two(app):
    seed = {"A01": [2135, 2214, 2425, 2520], "A02": [2172, 2266, 2477, 2563]}
    for well, scans in seed.items():
        pk = app._ensure_picker(RUN / f"{well}.rsd")
        pk.std = [(x, "") for x in scans]
    app.semi_auto_is()
    _run_auto(app)


def _all_buttons(win):
    out = []
    for child in win.winfo_children():
        try:
            if isinstance(child, ttk.Button):
                out.append(child)
        except Exception:
            pass
        out.extend(_all_buttons(child))
    return out


def test_auto_table_has_its_own_save_button(app):
    _mark_two(app)
    labels = []
    for b in _all_buttons(app.center):
        try:
            labels.append(b.cget("text"))
        except Exception:
            pass
    assert any("Save table" in t for t in labels), labels


def test_table_shows_flags_and_is_peak_positions(app):
    _mark_two(app)
    idx = {c: i for i, c in enumerate(app._auto_tree["columns"])}
    assert "ispeaks" in idx and "flags" in idx, idx
    first = app._auto_tree.item(app._auto_tree.get_children()[0])["values"]
    # the IS peaks column carries the standard's scan positions, not a blank
    assert first[idx["ispeaks"]].replace("/", "").isdigit(), first
    # the two hand-marked wells are the tinted rows
    assert "manual" in app._auto_tree.item(
        app._auto_tree.get_children()[0], "tags")


def test_double_click_edits_a_call_and_the_export_writes_the_override(
        app, tmp_path, monkeypatch):
    _mark_two(app)
    tree = app._auto_tree
    well0 = tree.item(tree.get_children()[0], "values")[0]
    original = app._auto_rows[0]["call"]
    assert original != "no-call", (original, well0)

    # closing an open editor without applying leaves the call untouched
    app._auto_edit_call_at(0)
    assert app._auto_edit_win is not None
    data = dict(app._auto_rows[0])
    app._auto_edit_close()
    assert app._auto_edit_win is None
    assert app._auto_rows[0] == data

    app._auto_edit_call_at(0)
    win = app._auto_edit_win
    assert win is not None
    # the dialog lists the engine's vocabulary
    om = win.winfo_children()[1]
    assert isinstance(om, ttk.Combobox)
    win._var.set("no-call")
    win._apply()
    assert app._auto_edit_win is None
    assert app._auto_rows[0]["call"] == "no-call"
    assert tree.item(tree.get_children()[0], "values")[1] == "no-call"

    out = tmp_path / "auto.csv"
    monkeypatch.setattr(filedialog, "asksaveasfilename",
                        lambda *a, **k: str(out))
    app.auto_genotype_save()
    with open(out, newline="") as f:
        rows = list(csv.DictReader(f))
    by_well = {r["well"]: r["call"] for r in rows}
    assert by_well[well0] == "no-call", (by_well, well0, rows)

def test_marks_alone_are_enough_to_genotype_the_whole_plate(app):
    """Mark a few wells, press Auto, get the plate back -- no extra step.

    Semi automatic standard… was a separate requirement nobody enforced, so
    marking the IS and running Auto gave the marked wells and a wall of
    no-calls.  Auto now fits the same model from the same marks itself and
    writes the placements into the pick table, so the two hand marks carry
    the whole plate -- and it still has to agree with the manual ground
    truth in tests/data/rs1695_expected.csv.
    """
    app._auto_is_channel.set(genotyping.DEFAULT_IS_CHANNEL)
    app._auto_sample_channel.set(genotyping.DEFAULT_SAMPLE_CHANNEL)
    seed = {"A01": [2135, 2214, 2425, 2520], "A02": [2172, 2266, 2477, 2563]}
    for well, scans in seed.items():
        pk = app._ensure_picker(RUN / f"{well}.rsd")
        pk.std = [(x, "") for x in scans]

    assert app._auto_is_model is None          # Semi automatic was never run
    app.auto_genotype_wells()

    assert app._auto_is_model is not None, "Auto did not learn from the marks"
    assert app._auto_is_model_note
    expected = {}
    with open(Path(__file__).parent / "data" / "rs1695_expected.csv",
              newline="") as f:
        for row in csv.DictReader(f):
            expected[row["well"]] = row["expected"]
    got = {r["well"]: r["call"].replace("hom-", "hom") for r in app._auto_rows}
    assert got == {w: expected[w] for w in WELLS}, (got, expected)
    # the two unmarked wells were placed, and are visible as marks
    for r in app._auto_rows:
        assert r["std_scans"], r
        assert r["std_source"] != "none", r
    for well in ("A03", "A04"):
        pk = app._ensure_picker(RUN / f"{well}.rsd")
        assert pk.std and pk.std_auto, (well, pk.std, pk.std_auto)


def test_marks_without_enough_of_them_invent_nothing(app):
    """One marked well is not a shape to fit: nothing spreads from it."""
    app._auto_is_channel.set(genotyping.DEFAULT_IS_CHANNEL)
    app._auto_sample_channel.set(genotyping.DEFAULT_SAMPLE_CHANNEL)
    pk = app._ensure_picker(RUN / "A01.rsd")
    pk.std = [(x, "") for x in (2135, 2214, 2425, 2520)]
    app.auto_genotype_wells()
    assert app._auto_is_model is None
    got = {r["well"]: r for r in app._auto_rows}
    assert got["A01"]["std_scans"], got["A01"]   # the mark itself is reported
    for well in ("A02", "A03", "A04"):
        assert got[well]["std_scans"] == "", got[well]
    assert any(r["call"] == "no-call" for r in app._auto_rows), got


def test_the_marked_wells_report_counts_and_lists_the_marks(app, monkeypatch):
    shown = []
    monkeypatch.setattr(messagebox, "showinfo",
                        lambda title, msg="", **k: shown.append(msg))
    app.report_marked_wells()
    assert "No wells are marked" in shown[0], shown

    shown.clear()
    _mark_two(app)               # also places the standard on A03 and A04
    app.report_marked_wells()
    text = shown[-1]             # _mark_two's own semi-automatic dialog is in
    assert "4 wells marked" in text, text   # shown[0]; this is the report
    for well in ("A01", "A02", "A03", "A04"):
        assert well in text, text
    assert "hand-marked" in text, text
    assert "Learned:" in text, text


def test_the_marked_wells_is_a_menu_item(app):
    gen = _find_genotyping_menu(app)
    assert ("command", "Marked IS wells…") in _entries(gen), _entries(gen)


def test_auto_targets_falls_back_to_the_whole_sample_window(app):
    """No selection means the whole list, not an error.

    Selecting 96 files one at a time to genotype a plate was pure friction:
    the sample window already lists exactly the run being worked on.
    """
    app.selected = []
    assert app._auto_targets() == [Path(p) for p in app.files]
    app.selected = [RUN / "A02.rsd"]
    assert app._auto_targets() == [RUN / "A02.rsd"]


def test_marks_are_saved_and_restored_after_a_restart(app):
    """Marks belong on disk, not in the window: a restart used to throw away
    the work of marking a plate, and the next Auto-genotype had nothing to
    learn from."""
    import mark_library
    seed = {"A01": [2135, 2214, 2425, 2520], "A02": [2172, 2266, 2477, 2563]}
    for well, scans in seed.items():
        pk = app._ensure_picker(RUN / f"{well}.rsd")
        pk.std = [(x, "") for x in scans]
    app._library_store(*[(RUN / f"{w}.rsd", app._ensure_picker(RUN / f"{w}.rsd"))
                         for w in seed])
    saved = mark_library.load_run(RUN.name)
    assert set(seed) <= set(saved), saved
    for well, scans in seed.items():
        assert [p[0] for p in saved[well]["std"]] == scans, (well, saved[well])

    # restart: every in-memory picker, cache and model is gone
    app._gen_pickers.clear()
    app._library_cache.clear()
    app._auto_is_model = None
    app._auto_models.clear()

    for well, scans in seed.items():
        pk = app._ensure_picker(RUN / f"{well}.rsd")
        assert [int(x) for x, _ in pk.std] == scans, (well, pk.std)


def test_auto_genotype_reseeds_itself_from_the_library(app):
    """Mark a plate, restart, press Auto: the marks on disk are the seeds.

    This is the whole point of the library -- without it the restart lost the
    marks and the plate came back as a wall of no-calls."""
    import mark_library
    seed = {"A01": [2135, 2214, 2425, 2520], "A02": [2172, 2266, 2477, 2563]}
    wells = {}
    for well, scans in seed.items():
        pk = app._ensure_picker(RUN / f"{well}.rsd")
        pk.std = [(x, "") for x in scans]
        wells[well] = mark_library.entry_from_picker(pk)
    mark_library.save_run(RUN.name, wells)

    # restart: only the file survives
    app._gen_pickers.clear()
    app._library_cache.clear()
    app._auto_is_model = None
    app._auto_models.clear()

    _run_auto(app)

    assert app._auto_is_model is not None, "Auto did not learn from the file"
    expected = {}
    with open(Path(__file__).parent / "data" / "rs1695_expected.csv",
              newline="") as f:
        for row in csv.DictReader(f):
            expected[row["well"]] = row["expected"]
    got = {r["well"]: r["call"].replace("hom-", "hom") for r in app._auto_rows}
    assert got == {w: expected[w] for w in WELLS}, (got, expected)
    for r in app._auto_rows:
        assert r["std_scans"], r


def test_mark_peaks_as_standard_marks_every_well_i_picked_in(app,
                                                             monkeypatch):
    """Pick four peaks on four wells, press the command once: four marks.

    It used to mark only the last well clicked, so three wells stayed
    unmarked with nothing said -- and the auto-genotype had no seeds for
    them."""
    import numpy as np
    import mark_library
    app.enter_genotyping_picking()
    app.n_graphs.set(4)
    is_col = genotyping.acgt_index_for_channel(app.base_order_var.get(),
                                               genotyping.DEFAULT_IS_CHANNEL)
    for path in app.selected:
        doc = app._ensure_doc(path)
        found = genotyping.find_is_quartet(
            np.asarray(doc.acgt, dtype=float)[:, is_col])
        assert found, path.name
        peaks, _heights = found
        pk = app._ensure_picker(path)
        for scan in peaks:
            assert pk.pick(float(scan), only_col=is_col) is not None, \
                (path.name, scan, pk._reject)

    monkeypatch.setattr(simpledialog, "askstring", lambda *a, **k: "")
    app._gen_mark_std()

    for well in WELLS:
        assert app._ensure_picker(RUN / f"{well}.rsd").std, well
    assert "Standard set in 4 well(s)" in app.status_var.get(), \
        app.status_var.get()
    saved = mark_library.load_run(RUN.name)
    assert all(saved[w].get("std") for w in WELLS), sorted(saved)


# --------------------------------------------------------------------------- #
# correcting the table: the call, and the standard the call was read against
# --------------------------------------------------------------------------- #
def _row_index(app, well):
    for i, r in enumerate(app._geno_rows):
        if r.get("well") == well:
            return i
    raise AssertionError(f"{well} is not in the table: "
                         f"{[r.get('well') for r in app._geno_rows]}")


def test_the_table_addresses_its_columns_by_name(app):
    """The editors hit columns through the table's own column list.

    A Treeview numbers its columns by position (``"#2"`` is the second one
    here), so hard-coding an id in the click handler would silently edit the
    wrong cell as soon as a column is added."""
    _mark_two(app)
    assert app._geno_col("call") == "#2"
    assert app._geno_col("ispeaks") == "#12"
    assert app._geno_col("reason") == f"#{len(app._auto_cols)}"
    assert app._geno_col("no-such-column") == ""


def test_editing_a_library_only_row_takes_it_over(app):
    """A table rebuilt from the library is editable, not a dead copy.

    After a restart the table's rows come off disk and are not in the session
    list at all, so editing one used to change a dict nothing could reload --
    the correction vanished on the next save."""
    import mark_library
    _mark_two(app)                       # the library now holds all four wells
    app.selected = [RUN / "A01.rsd"]     # a fresh session over the same plate
    app.files = app.selected
    app._auto_rows = []
    app._build_auto_table()

    i = _row_index(app, "A03")
    assert all(r.get("well") != "A03" for r in app._auto_rows), "library row"
    app._auto_edit_call_at(i)
    win = app._auto_edit_win
    assert win is not None
    win._var.set("no-call")
    win._apply()

    assert any(r.get("well") == "A03" and r["call"] == "no-call"
               for r in app._auto_rows), app._auto_rows
    saved = mark_library.load_run(RUN.name)
    assert saved["A03"]["result"]["call"] == "no-call", saved["A03"]
    assert saved["A03"]["result"]["label_source"] == "edited"


def test_double_click_on_the_is_scans_corrects_the_standard(app):
    """The four scans in the table are the standard, and they are editable.

    Correcting them re-reads the well against exactly those peaks and keeps
    both the mark and the new answer on disk, so the row describes the
    standard now on screen rather than the one before."""
    import mark_library
    _mark_two(app)                        # also places the standard on A03

    i = _row_index(app, "A03")
    row = app._geno_row(i)
    before = app._std_scans_of(row)
    assert before, before                    # the row carries a standard
    assert row["std_source"] != "manual", row   # placed by the model

    app._auto_edit_is_at(i)
    win = app._auto_edit_win
    assert win is not None and len(win._is_vars) == 4
    assert [v.get() for v in win._is_vars][:len(before)] == \
        [str(s) for s in before]
    assert all(not v.get() for v in win._is_vars[len(before):])

    # leaving the dialog alone changes nothing
    app._auto_edit_close()
    assert app._auto_edit_win is None
    assert app._std_scans_of(app._geno_row(i)) == before

    # four numbers that are not four different, increasing scans are refused
    app._auto_edit_is_at(i)
    win = app._auto_edit_win
    for v, s in zip(win._is_vars, (2100, 2100, 2400, 2500)):
        v.set(str(s))
    win._apply()
    assert app._auto_edit_win is not None      # still open, nothing written
    assert app._std_scans_of(app._geno_row(i)) == before

    # four ascending scans re-read the well against them.  A well whose
    # bands partly merged comes through with fewer than four, so the
    # operator fills the rest in -- which is the point of the editor: it
    # says where the peaks are, not what the model managed to find.
    shifted = [s + 1 for s in before]
    while len(shifted) < 4:
        shifted.append(shifted[-1] + 100)
    for v, s in zip(win._is_vars, shifted):
        v.set(str(s))
    win._apply()
    assert app._auto_edit_win is None
    row = app._geno_row(i)
    assert app._std_scans_of(row) == shifted
    assert row["std_source"] == "manual", row

    saved = mark_library.load_run(RUN.name)
    assert [int(s) for s, _ in saved["A03"]["std"]] == shifted, saved["A03"]
    assert saved["A03"]["result"]["std_scans"] == "/".join(str(s)
                                                           for s in shifted)
    assert saved["A03"]["result"]["label_source"] == "edited"


def _cell_centre(bbox, k):
    """Centre of a cell along axis *k* (0=x, 1=y).

    identify_column()'s boundary sits a few pixels right of bbox()'s, so only
    the middle of a cell is unambiguously in that column -- clicking at the
    left edge lands in the previous one."""
    return bbox[k] + bbox[k + 2] // 2


def _double_click(widget, x, y):
    """Two real button presses, so Tk fires the <Double-1> itself."""
    for _ in range(2):
        widget.event_generate("<ButtonPress-1>", x=x, y=y)
        widget.event_generate("<ButtonRelease-1>", x=x, y=y)
    widget.update()


def test_a_real_double_click_on_a_table_cell_opens_its_editor(app):
    """The handler has to work from a real click, not only from a direct call.

    The click lands on a cell of a live, mapped Treeview, so this is the path
    the operator actually takes: identify the row and column under the pointer
    and open the matching editor -- or, if that mapping is wrong, silently do
    nothing, which a direct call to the editor cannot catch."""
    _mark_two(app)
    app.geometry("1600x1000")
    tree = app._auto_tree
    app.update()
    iid = tree.get_children()[0]

    tree.see(iid)
    app.update()
    bbox = tree.bbox(iid, app._geno_col("call"))
    assert bbox, "the call cell must be on screen to click it"
    _double_click(tree, _cell_centre(bbox, 0), _cell_centre(bbox, 1))
    assert app._auto_edit_win is not None, "a double-click on call did nothing"
    assert isinstance(app._auto_edit_win.winfo_children()[1], ttk.Combobox)
    app._auto_edit_close()

    tree.see(iid)
    app.update()
    bbox = tree.bbox(iid, app._geno_col("ispeaks"))
    assert bbox, "the IS-peaks cell must be on screen to click it"
    _double_click(tree, _cell_centre(bbox, 0), _cell_centre(bbox, 1))
    assert app._auto_edit_win is not None, "a double-click on IS did nothing"
    assert len(app._auto_edit_win._is_vars) == 4
    app._auto_edit_close()
