"""Hover help for Tk widgets.

:class:`ToolTip` shows a small yellow note after the pointer rests on a
widget; :func:`install_tooltips` walks a window and gives every button,
check box, menu button, drop-down and spin box a note from a help table
keyed by the widget's own label, so a button's help is written once and
follows the button wherever it is packed.
"""
from __future__ import annotations

import re
import tkinter as tk
from typing import Callable, Dict, Iterable, Optional

DELAY_MS = 450


class ToolTip:
    """Small floating label that appears after hovering a widget.
    `getter` is called on hover and returns the text to show (or '')."""

    def __init__(self, widget, getter: Callable[[], str]):
        self.widget = widget
        self.getter = getter
        self.tip = None
        self._after = None
        widget._has_tooltip = True
        widget.bind("<Enter>", lambda e: self._schedule(), add="+")
        widget.bind("<Leave>", lambda e: self._cancel(), add="+")
        widget.bind("<ButtonPress>", lambda e: self._cancel(), add="+")
        widget.bind("<Motion>", lambda e: self._schedule(), add="+")

    def _schedule(self):
        # Never stack tips: drop any visible one and re-arm.
        self._cancel()
        self._after = self.widget.after(DELAY_MS, self._show)

    def _cancel(self):
        if self._after:
            try:
                self.widget.after_cancel(self._after)
            except tk.TclError:
                pass
            self._after = None
        if self.tip is not None:
            try:
                self.tip.destroy()
            except tk.TclError:
                pass
            self.tip = None

    def _show(self):
        self._after = None
        try:
            text = self.getter()
        except Exception:                                  # noqa: BLE001
            text = ""
        if not text:
            return
        if self.tip is not None:
            self.tip.destroy()
            self.tip = None
        self.tip = tk.Toplevel(self.widget, bg="#FFFFE1")
        self.tip.wm_overrideredirect(True)
        tk.Label(self.tip, text=text, bg="#FFFFE1", fg="#222222",
                 justify=tk.LEFT, font=("DejaVu Sans", 9), wraplength=420,
                 padx=4, pady=2).pack()
        # below-right of the pointer, but never off the screen's edge
        self.tip.update_idletasks()
        w, h = self.tip.winfo_reqwidth(), self.tip.winfo_reqheight()
        sw, sh = self.widget.winfo_screenwidth(), self.widget.winfo_screenheight()
        x = min(self.widget.winfo_pointerx() + 14, sw - w - 4)
        y = self.widget.winfo_pointery() + 14
        if y + h > sh - 4:
            y = self.widget.winfo_pointery() - h - 8
        self.tip.wm_geometry(f"+{max(0, x)}+{max(0, y)}")


def label_key(text: str) -> str:
    """A widget label reduced for lookup: symbols and spacing dropped,
    lower case ("  ⊕ Add…  " -> "add", "▶ Tour" -> "tour")."""
    t = re.sub(r"[^\w()%≥µ/+\-–−]+", " ", str(text or ""))
    t = re.sub(r"\s+", " ", t).strip()
    t = re.sub(r"^[–\-−]\s+(?=\w)", "", t)   # "– Hide file(s)" -> "Hide file(s)"
    return t.lower()


_KINDS = {"Button", "TButton", "Checkbutton", "TCheckbutton", "Radiobutton",
          "TRadiobutton", "Menubutton", "TMenubutton", "TCombobox", "Spinbox",
          "TSpinbox", "Entry", "TEntry"}


def _text_of(w) -> str:
    try:
        return str(w.cget("text"))
    except (tk.TclError, AttributeError):
        return ""


def install_tooltips(root, helps: Dict[str, str],
                     patterns: Iterable = (), by_widget: Optional[Dict] = None) -> int:
    """Give every button-like widget under *root* hover help.

    *helps* maps :func:`label_key` of a widget's label to its note;
    *patterns* are ``(regex, note)`` pairs tried on the label next (e.g.
    every "Ch3 T" channel switch); *by_widget* gives notes for widgets
    without a label (drop-downs, spin boxes).  Widgets that already have a
    tip are left alone, so this can be called again after a bar is built.
    Returns how many tips were added."""
    helps = {label_key(k): v for k, v in helps.items()}
    pats = [(re.compile(p), t) for p, t in patterns]
    by_widget = by_widget or {}
    added = 0
    stack = [root]
    while stack:
        w = stack.pop()
        try:
            stack.extend(w.winfo_children())
        except tk.TclError:
            continue
        if getattr(w, "_has_tooltip", False):
            continue
        note = by_widget.get(w)
        if note is None:
            try:
                kind = w.winfo_class()
            except tk.TclError:
                continue
            if kind not in _KINDS:
                continue
            key = label_key(_text_of(w))
            note = helps.get(key) if key else None
            if note is None and key:
                note = next((t for p, t in pats if p.search(key)), None)
        if note:
            ToolTip(w, (lambda n=note: n() if callable(n) else n))
            added += 1
    return added


def missing_tooltips(root) -> list:
    """Labels of button-like widgets under *root* that still have no tip
    (for tests and for spotting a new button that needs a note)."""
    out, stack = [], [root]
    while stack:
        w = stack.pop()
        try:
            stack.extend(w.winfo_children())
            kind = w.winfo_class()
        except tk.TclError:
            continue
        if kind in {"Button", "TButton", "Checkbutton", "TCheckbutton",
                    "Menubutton", "TMenubutton"} \
                and not getattr(w, "_has_tooltip", False):
            out.append(_text_of(w).strip() or kind)
    return out
