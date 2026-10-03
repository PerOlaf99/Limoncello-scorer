"""GUI entry point for fragment-length sizing.

``LimoncelloAnalyzerApp.size_selected_fragments`` is thin: it asks for the
ladder and the two physical channels, sizes every selected well with the same
headless engine the CLI uses, and writes CSV/JSON/HTML.  The dialogs are
monkeypatched and the app is a stub, so the branching (cancel, write failure,
success) is covered without a display.
"""
from pathlib import Path

import numpy as np
import pytest

tkinter = pytest.importorskip("tkinter")
matplotlib = pytest.importorskip("matplotlib")
pytest.importorskip("tkinter.ttk")

import sequence_analyzer as sa  # noqa: E402
import fragment_sizing as fs  # noqa: E402
from analyzer_core import acgt_index_for_channel  # noqa: E402

N_SCANS = 6000
SCAN_OF = lambda bp: 1200.0 + 260.0 * np.log(bp)  # noqa: E731
TRUE_BP = [90.0, 180.0, 275.0, 420.0]


class _Var:
    def __init__(self, value=""):
        self.v = value

    def get(self):
        return self.v

    def set(self, value):
        self.v = value


class FakeDoc:
    def __init__(self, path):
        ladder = fs.load_ladder("genescan500_rox")
        x = np.arange(N_SCANS)
        acgt = np.zeros((N_SCANS, 4)) + 5.0
        for bp, seed in ((TRUE_BP, 0), (ladder.lengths, 1)):
            col = (acgt_index_for_channel("ACTG", 2) if bp is TRUE_BP
                   else acgt_index_for_channel("ACTG", 4))
            rng = np.random.default_rng(seed)
            y = acgt[:, col] + rng.normal(0.0, 2.0, N_SCANS)
            for b in bp:
                y += 850.0 * np.exp(-0.5 * ((x - SCAN_OF(b)) / 6.0) ** 2)
            acgt[:, col] = y
        self.acgt = acgt
        self.well = Path(path).stem
        self.path = path
        self.base_order = "ACTG"


class StubApp:
    def __init__(self, docs):
        self.selected = list(docs)
        self._docs = {str(Path(p)): FakeDoc(p) for p in docs}
        self.base_order_var = _Var("ACTG")
        self.status_var = _Var("")

    def _ensure_doc(self, path):
        return self._docs[str(Path(path))]


@pytest.fixture
def stub_app(tmp_path):
    return StubApp([tmp_path / "A01.rsd", tmp_path / "A02.rsd"])


def _patch_dialogs(monkeypatch, *, ladder="genescan500_rox",
                   channels=(4, 2), out="", askstring_cancel=False,
                   file_cancel=False):
    calls = {"info": [], "error": []}
    monkeypatch.setattr(sa.simpledialog, "askstring",
                        lambda *a, **k: None if askstring_cancel else ladder)
    it = iter(channels)
    monkeypatch.setattr(sa.simpledialog, "askinteger",
                        lambda *a, **k: next(it))
    monkeypatch.setattr(sa.filedialog, "asksaveasfilename",
                        lambda *a, **k: "" if file_cancel else out)
    monkeypatch.setattr(sa.messagebox, "showinfo",
                        lambda *a, **k: calls["info"].append(a))
    monkeypatch.setattr(sa.messagebox, "showerror",
                        lambda *a, **k: calls["error"].append(a))
    return calls


def test_sizes_selected_wells_to_csv(monkeypatch, stub_app, tmp_path):
    out = str(tmp_path / "sizes.csv")
    calls = _patch_dialogs(monkeypatch, out=out)
    sa.LimoncelloAnalyzerApp.size_selected_fragments(stub_app)

    assert Path(out).is_file()
    text = Path(out).read_text(encoding="utf-8-sig")
    assert "length_bp" in text
    assert "A01" in text and "A02" in text
    assert not calls["error"]
    assert "Sized" in stub_app.status_var.get()


def test_cancelling_the_ladder_prompt_writes_nothing(monkeypatch, stub_app, tmp_path):
    out = str(tmp_path / "sizes.csv")
    calls = _patch_dialogs(monkeypatch, out=out, askstring_cancel=True)
    sa.LimoncelloAnalyzerApp.size_selected_fragments(stub_app)
    assert not Path(out).exists()
    assert not calls["error"]


def test_cancelling_the_save_dialog_writes_nothing(monkeypatch, stub_app, tmp_path):
    out = str(tmp_path / "sizes.csv")
    _patch_dialogs(monkeypatch, out=out, file_cancel=True)
    sa.LimoncelloAnalyzerApp.size_selected_fragments(stub_app)
    assert not Path(out).exists()


def test_writes_an_html_report_when_asked(monkeypatch, stub_app, tmp_path):
    out = str(tmp_path / "sizes.html")
    _patch_dialogs(monkeypatch, out=out)
    sa.LimoncelloAnalyzerApp.size_selected_fragments(stub_app)
    html_text = Path(out).read_text(encoding="utf-8")
    assert "Leave-one-out RMS error" in html_text
    assert "<table" in html_text


def test_custom_lengths_ladder_is_accepted(monkeypatch, stub_app, tmp_path):
    out = str(tmp_path / "sizes.csv")
    _patch_dialogs(monkeypatch, ladder="50,100,150,200,250", out=out)
    sa.LimoncelloAnalyzerApp.size_selected_fragments(stub_app)
    assert Path(out).is_file()
