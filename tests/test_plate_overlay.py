"""Plate overlay: every well aligned on its standard, coloured by call."""
import numpy as np
import pytest

import plate_overlay as po
import snp_genotyping as sg
from test_snp_genotyping import ORDER, _plate


@pytest.fixture(scope="module")
def plate():
    return _plate(seed=11)


def test_overlay_lines_up_the_sample(plate):
    truth, arrs = plate
    data = po.overlay_data(arrs, ORDER, "piecewise", "sample")
    assert len(data["traces"]) >= 90 and len(data["ref"]) == 4
    tops = [x[np.argmax(y)] for w, (x, y) in data["traces"].items() if truth[w] == "hom1"]
    assert np.std(tops) < 4
    assert all(data["calls"][w] == truth[w] for w in truth if truth[w] != "none")
    lo, hi = data["xlim"]
    assert lo < data["ref"][0] and hi > data["ref"][-1]


def test_not_aligned_shows_the_recorded_spread(plate):
    truth, arrs = plate
    rows = sg.genotype_plate(arrs, ORDER)
    data = po.overlay_data(arrs, ORDER, "none", "sample", rows=rows)
    tops = [x[np.argmax(y)] for w, (x, y) in data["traces"].items() if truth[w] == "hom1"]
    assert np.std(tops) > 50


def test_standard_channel_shows_the_four_bands(plate):
    truth, arrs = plate
    data = po.overlay_data(arrs, ORDER, "piecewise", "standard")
    x, y = next(iter(data["traces"].values()))
    for r in data["ref"]:
        near = (x > r - 4) & (x < r + 4)
        assert y[near].max() > 0.5 * y.max()


def test_plate_without_a_standard_draws_nothing():
    rng = np.random.default_rng(2)
    arrs = {f"A{c:02d}": rng.normal(100, 5, (4000, 4)) for c in range(1, 13)}
    data = po.overlay_data(arrs, ORDER)
    assert data["traces"] == {}
    assert "no internal standard" in data["note"]


def test_window_draws_and_a_click_shows_the_well(plate):
    tk = pytest.importorskip("tkinter")
    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("no display")
    root.withdraw()
    truth, arrs = plate
    shown = []
    try:
        win = po.PlateOverlay(
            root, runs=lambda: ["RUN1"],
            files_for=lambda run: {w: f"/data/{run}/{w}.rsd" for w in arrs},
            load=lambda p: arrs[p.rsplit("/", 1)[1][:-4]],
            order_for=lambda d: ORDER, show=shown.append, status=lambda t: None)
        lines = [ln for ln in win.ax.get_lines() if ln.get_gid()]
        assert len(lines) >= 90
        assert "standard on Ch3" in win.note.cget("text")

        class Ev:
            artist = lines[0]
        win._picked(Ev())
        assert shown == [[f"/data/RUN1/{lines[0].get_gid()}.rsd"]]
        win.method_var.set(po.METHODS[2][1])               # the 2005 sum normalisation
        win.refresh()
        assert len([ln for ln in win.ax.get_lines() if ln.get_gid()]) >= 90
    finally:
        root.destroy()
