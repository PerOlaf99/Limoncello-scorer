"""The auto-genotyping GUI: the menu item, the channel-role check, the results
table and the export.  Skips where there is no display to put a window on."""
import csv
import os
from pathlib import Path

import pytest

pytest.importorskip("tkinter")
pytest.importorskip("matplotlib")

if not os.environ.get("DISPLAY"):
    pytest.skip("no DISPLAY for the GUI test")

import tkinter as tk  # noqa: E402
from tkinter import filedialog, messagebox  # noqa: E402

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


def test_the_real_t9_wells_match_the_recorded_ground_truth(app):
    """The A-row calls through the GUI agree with tests/data/rs1695_expected.csv,
    which is the manual ground truth the 95/96 validation was measured against.
    Call names differ between the fixture ('hom2') and the engine ('hom-2')."""
    expected = {}
    with open(Path(__file__).parent / "data" / "rs1695_expected.csv",
              newline="") as f:
        for row in csv.DictReader(f):
            expected[row["well"]] = row["expected"]
    app._auto_is_channel.set(genotyping.DEFAULT_IS_CHANNEL)
    app._auto_sample_channel.set(genotyping.DEFAULT_SAMPLE_CHANNEL)
    app.auto_genotype_wells()
    got = {r["well"]: r["call"].replace("hom-", "hom") for r in app._auto_rows}
    want = {w: expected.get(w) for w in WELLS}
    assert got == want, (got, want)