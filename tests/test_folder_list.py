"""Data folders vs the Samples list: added runs are only listed until opened;
double-click / Enter opens them, the broom empties the list; run comments
come from the instrument and from a note kept in the run folder."""
import os
import shutil
from pathlib import Path

import pytest

tk = pytest.importorskip("tkinter")
if not os.environ.get("DISPLAY"):
    pytest.skip("no DISPLAY for the GUI test", allow_module_level=True)

import sequence_analyzer as sa  # noqa: E402

M13 = Path(__file__).resolve().parents[1] / "example_data" / "M13"
WELLS = ("A01", "A02")
RUNS = ("PLATE_Run01", "PLATE_Run02", "PLATE_Run03")


@pytest.fixture()
def project(tmp_path):
    if not all((M13 / f"{w}.rsd").exists() for w in WELLS):
        pytest.skip("example M13 traces not present")
    root = tmp_path / "project"
    for run in RUNS:
        (root / run).mkdir(parents=True)
        for w in WELLS:
            shutil.copy(M13 / f"{w}.rsd", root / run / f"{w}.rsd")
    return root


@pytest.fixture()
def app(tmp_path, monkeypatch):
    import mark_library
    monkeypatch.setattr(mark_library, "LIBRARY_DIR", tmp_path / "library")
    try:
        a = sa.LimoncelloAnalyzerApp(initial_folders=[])
    except tk.TclError as e:
        pytest.skip(f"cannot open a window: {e}")
    a.update()
    yield a
    a.destroy()


def _runs_listed(app):
    return sorted({f.parent.name for f in app.files})


def test_added_runs_wait_in_data_folders(app, project):
    app._absorb_runs(project)
    app.refresh_file_list()
    assert [p.name for p in app._folder_paths] == list(RUNS)
    assert app.files == []


def test_open_two_runs_clear_and_open_the_next(app, project):
    app._absorb_runs(project)
    app.refresh_file_list()
    app.folder_list.selection_set(0, 1)              # Run01 + Run02, then Enter
    app.open_selected_runs()
    assert _runs_listed(app) == ["PLATE_Run01", "PLATE_Run02"]
    assert len(app.files) == 4
    assert app.folder_list.get(0).startswith("✓")
    assert not app.folder_list.get(2).startswith("✓")
    app.clear_file_list()                            # the broom
    assert app.files == [] and app.selected == []
    assert len(app._folder_paths) == 3               # runs stay available
    app.folder_list.selection_clear(0, tk.END)
    app.folder_list.selection_set(2)
    app.open_selected_runs()
    assert _runs_listed(app) == ["PLATE_Run03"]


def test_double_click_toggles_a_run(app, project):
    app._absorb_runs(project)
    app.refresh_file_list()
    app.update()

    class Ev:
        y = app.folder_list.bbox(1)[1] + 2           # row 2: Run02
    app._on_folder_activate(Ev())
    assert _runs_listed(app) == ["PLATE_Run02"]
    assert app.selected and app.selected[0].parent.name == "PLATE_Run02"
    app._on_folder_activate(Ev())                    # again: out of the list
    assert app.files == []


def test_command_line_folders_are_listed_at_once(tmp_path, project, monkeypatch):
    import mark_library
    monkeypatch.setattr(mark_library, "LIBRARY_DIR", tmp_path / "library")
    a = sa.LimoncelloAnalyzerApp(initial_folders=[project / "PLATE_Run01"])
    try:
        a.update()
        assert _runs_listed(a) == ["PLATE_Run01"]
    finally:
        a.destroy()


def test_run_note_and_instrument_comment(app, project):
    app._absorb_runs(project, activate=True)
    app.refresh_file_list()
    path = app.files[0]

    class Doc:
        footer_fields = {"COMMENT": "Hel plate 829\r\nGrad 50_65C"}
    assert app._instrument_comment(Doc()) == "Hel plate 829\nGrad 50_65C"
    assert app._read_comment(path) == ""
    app._run_note_file(path).write_text("re-run B05\n", encoding="utf-8")
    assert app._read_comment(path) == "re-run B05"
    # the same note for every well of the run
    other = [f for f in app.files if f.parent == path.parent and f != path][0]
    assert app._read_comment(other) == "re-run B05"
    summary = app._comment_summary(path, Doc())
    assert summary.splitlines() == ["Hel plate 829", "Grad 50_65C", "re-run B05"]
