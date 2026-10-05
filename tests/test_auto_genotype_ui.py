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
from tkinter import filedialog, messagebox, simpledialog  # noqa: E402

import analyzer_core  # noqa: E402
import genotyping  # noqa: E402
import scorer  # noqa: E402
from sequence_analyzer import LimoncelloAnalyzerApp  # noqa: E402

RUN = Path("/media/per/78B0C7DE1FA7081C/OY/OY_rs1695_T9_270910Run01")
WELLS = ["A01", "A02", "A03", "A04"]


@pytest.fixture()
def app(monkeypatch):
    paths = [RUN / f"{w}.rsd" for w in WELLS]
    if not all(p.exists() for p in paths):
        pytest.skip(f"{RUN} not present")
    for name in ("showinfo", "showwarning", "showerror"):
        monkeypatch.setattr(messagebox, name, lambda *a, **k: "ok")
    monkeypatch.setattr(filedialog, "asksaveasfilename", lambda *a, **k: "")
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
    assert "Auto-genotype selected wells…" in sub_labels, sub_labels
    assert "Channel roles (standard / sample)…" in sub_labels, sub_labels
    assert "Save auto-genotype table…" in sub_labels, sub_labels


def test_mark_all_refuses_without_a_learned_shape(app, monkeypatch):
    """96 unmarked wells is not a licence to guess each one separately.

    The batch command used to fall back to bare equimolar detection, which is
    what put noise in a third of the ABCC2 wells. With nothing learned it now
    says so and points at the semi-automatic item.
    """
    shown = []
    monkeypatch.setattr(messagebox, "showinfo",
                        lambda title, msg="", **k: shown.append(msg))
    monkeypatch.setattr(simpledialog, "askinteger",
                        lambda *a, **k: genotyping.DEFAULT_IS_CHANNEL)
    app.genotyping_active = True
    app._auto_is_model = None
    app._gen_mark_std_batch()
    assert shown and "Semi automatic standard" in shown[0], shown
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