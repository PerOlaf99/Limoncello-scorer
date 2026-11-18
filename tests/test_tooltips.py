"""Every button of the main window explains itself on hover, and the auto-tour
sits in the top bar where it can be found."""
import os

import pytest

from tooltips import label_key

tk = pytest.importorskip("tkinter")


def test_label_key_drops_symbols_and_spacing():
    assert label_key("  ⊕ Add…  ") == "add"
    assert label_key("▶ Tour") == "tour"
    assert label_key("⟲ Reset view") == "reset view"
    assert label_key("  Time (min)  ") == "time (min)"
    assert label_key("Align") == "align"
    assert label_key("−") == "−"


def test_help_keys_stay_distinct_once_normalised():
    import sequence_analyzer as sa
    keys = [label_key(k) for k in sa.BUTTON_HELP]
    assert len(keys) == len(set(keys))


@pytest.fixture()
def app(tmp_path, monkeypatch):
    if not os.environ.get("DISPLAY"):
        pytest.skip("no DISPLAY for the GUI test")
    import mark_library
    import sequence_analyzer as sa
    monkeypatch.setattr(mark_library, "LIBRARY_DIR", tmp_path / "library")
    try:
        a = sa.LimoncelloAnalyzerApp(initial_folders=[])
    except tk.TclError as e:
        pytest.skip(f"cannot open a window: {e}")
    a.update()
    yield a
    a.destroy()


def test_every_button_in_the_main_window_has_hover_help(app):
    from tooltips import missing_tooltips
    assert missing_tooltips(app) == []


def test_the_tour_button_is_in_the_top_bar_and_shows_its_state(app):
    btn = app._tour_btn
    assert btn.master is app._chan_bar
    assert btn.cget("text") == "▶ Tour"
    app.files = []                      # nothing to tour: it stops at once
    app.toggle_tour()
    assert btn.cget("text") == "▶ Tour" and not app.touring
