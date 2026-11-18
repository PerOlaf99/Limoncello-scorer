"""The user manual names menu items by their exact labels (``Menu ▸ Item``);
check every such path against the menus the app really builds, so the manual
cannot drift from the program again.  Skips where there is no display."""
import os
import re

import pytest

import user_manual

TOP = ("File", "View", "Base calling", "Manual genotyping", "Auto genotyping",
       "Fragment analysis",
       "Comments", "Help")


def test_manual_has_contents_and_plain_text():
    titles = [t for t, _ in user_manual.SECTIONS]
    assert len(titles) == len(set(titles))
    text = user_manual.plain_text()
    assert "Quick start" in text and "```" not in text


def test_code_blocks_are_closed():
    for title, body in user_manual.SECTIONS:
        assert body.count("```") % 2 == 0, title


def _menu_paths():
    pytest.importorskip("tkinter")
    if not os.environ.get("DISPLAY"):
        pytest.skip("no DISPLAY to build the menus")
    import tkinter as tk
    import sequence_analyzer as sa
    try:
        app = sa.LimoncelloAnalyzerApp()
    except tk.TclError as e:
        pytest.skip(f"cannot open a window: {e}")
    paths = []

    def walk(menu, path):
        end = menu.index("end")
        if end is None:
            return
        for i in range(end + 1):
            kind = menu.type(i)
            if kind in ("separator", "tearoff"):
                continue
            label = menu.entrycget(i, "label")
            if kind == "cascade":
                walk(app.nametowidget(menu.entrycget(i, "menu")), path + [label])
            else:
                paths.append(" ▸ ".join(path + [label]))
    for name, menu in app._menus.items():
        walk(menu, [name])
    app.destroy()
    return paths


def test_every_menu_path_in_the_manual_exists():
    paths = _menu_paths()
    assert set(TOP) == {p.split(" ▸ ")[0] for p in paths}
    text = user_manual.plain_text()
    # items, and submenus named on their own ("View ▸ Sort files")
    valid = set(paths)
    for p in paths:
        parts = p.split(" ▸ ")
        valid |= {" ▸ ".join(parts[:k]) for k in range(2, len(parts))}
    valid = sorted(valid, key=len, reverse=True)
    bad = []
    for m in re.finditer(r"(%s) ▸ " % "|".join(map(re.escape, TOP)), text):
        rest = text[m.start():m.start() + 120]
        if not any(rest.startswith(v) for v in valid):
            bad.append(rest.split("\n")[0][:80])
    assert not bad, bad
    # menu names that no longer exist
    assert "Genotyping menu" not in text
    assert not re.search(r"(?<!Manual )(?<!Auto )Genotyping ▸", text)
