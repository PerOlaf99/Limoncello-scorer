"""The Tag-duplex menu item end to end: MF has to reach the status bar, the
MF column of the pick table and the saved CSV, through the real app, on a real
well.  Skips where there is no display to put a window on."""
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
from types import SimpleNamespace  # noqa: E402
from tkinter import filedialog, messagebox  # noqa: E402

RUN2 = Path("/media/per/78B0C7DE1FA7081C/OY/OY_rs1695_T9_270910Run01")

import analyzer_core  # noqa: E402
import genotyping  # noqa: E402
from sequence_analyzer import LimoncelloAnalyzerApp  # noqa: E402

RUN = Path("/media/per/78B0C7DE1FA7081C/OY/OY_rs1695_T9_270910Run01")
WELL = RUN / "A01.rsd"


@pytest.fixture()
def app(tmp_path, monkeypatch):
    if not WELL.exists():
        pytest.skip(f"{WELL} not present")
    for name in ("showinfo", "showwarning", "showerror"):
        monkeypatch.setattr(messagebox, name, lambda *a, **k: "ok")
    monkeypatch.setattr(filedialog, "asksaveasfilename", lambda *a, **k: "")
    try:
        app = LimoncelloAnalyzerApp()
    except tk.TclError as e:                       # no usable display
        pytest.skip(f"cannot open a window: {e}")
    app.update()
    app.docs[str(WELL.resolve())] = analyzer_core.load_trace(
        WELL, base_order="ACTG")
    app.file_list.insert(tk.END, "run/A01.rsd")
    app.n_graphs.set(1)
    app.selected = [WELL]
    app.enter_genotyping_picking()
    app.update()
    app._gen_sh.set(False)
    app._sync_gen_opts()
    yield app
    app.destroy()


def _pick(app, scan, vol):
    app._on_gen_pick(SimpleNamespace(button=1, xdata=float(scan),
                                     ydata=float(vol),
                                     inaxes=app._plot_axes[0]))
    app.update()


def _tag_one_position(app):
    _pick(app, 2205, 20134)
    _pick(app, 2225, 3121)
    app.status_var.set("")
    app._gen_mark_duplex()
    app.update()


def test_status_bar_shows_the_duplex_species_and_mf(app):
    _tag_one_position(app)
    status = app.status_var.get()
    assert "Duplex species" in status
    assert "MF" in status
    assert "HOM1@2205" in status and "HOM2@2225" in status


def test_mf_column_of_the_pick_table(app):
    _tag_one_position(app)
    cols = list(app.pick_tree["columns"])
    assert cols[-1] == "MF", cols
    rows = [app.pick_tree.item(i)["values"] for i in app.pick_tree.get_children()]
    cells = [r[cols.index("MF")] for r in rows if r[cols.index("scan")] == 2205]
    assert cells and cells[0] not in ("", None)


def test_mf_reaches_the_saved_csv(app, tmp_path, monkeypatch):
    _tag_one_position(app)
    out = tmp_path / "peaks.csv"
    monkeypatch.setattr(genotyping.filedialog, "asksaveasfilename",
                        lambda *a, **k: str(out))
    app._gen_save()
    with open(out, newline="") as f:
        rows = list(csv.DictReader(f))
    assert "mf" in rows[0] and "ai" in rows[0]
    r = next(x for x in rows if x["scan"] == "2205")
    assert r["mf"] not in ("", None)
    assert abs(float(r["mf"]) - float(r["fraction"])) < 1e-3
