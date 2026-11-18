"""Tk windows for the ladder library: the manager, the editor, and the sizing
set-up dialog that picks a ladder and the two channels.

The storage and import logic is headless in :mod:`ladder_library`; this
module only presents it.
"""
from __future__ import annotations

import tkinter as tk
from tkinter import colorchooser, filedialog, messagebox, ttk
from typing import Optional

import ladder_library as lib
from fragment_sizing import Ladder

COMMON_DYES = ["ET-ROX", "ROX", "LIZ", "FAM", "HEX", "VIC", "JOE", "NED",
               "TAMRA", "TET", "Cy5", "WellRED D1", "custom"]


def _fmt(b: float) -> str:
    return f"{b:g}"


def describe(ladder: Ladder) -> str:
    lo, hi = ladder.lengths[0], ladder.lengths[-1]
    return (f"{ladder.name} — {ladder.dye}, {len(ladder.lengths)} bands, "
            f"{_fmt(lo)}–{_fmt(hi)} bp")


def _modal(win: tk.Toplevel, parent) -> None:
    win.transient(parent)
    win.update_idletasks()
    try:
        win.grab_set()
    except tk.TclError:                 # not viewable yet on some WMs
        pass
    win.focus_set()


# --------------------------------------------------------------------------- #
# editor
# --------------------------------------------------------------------------- #
class LadderEditor(tk.Toplevel):
    """Create or edit one ladder.  ``result`` is the Ladder, or None."""

    def __init__(self, parent, ladder: Optional[Ladder] = None,
                 title: str = "Ladder"):
        super().__init__(parent)
        self.title(title)
        self.result: Optional[Ladder] = None
        lad = ladder
        frm = ttk.Frame(self, padding=10)
        frm.pack(fill=tk.BOTH, expand=True)
        frm.columnconfigure(1, weight=1)

        def row(r, label, widget):
            ttk.Label(frm, text=label).grid(row=r, column=0, sticky="nw",
                                            padx=(0, 8), pady=3)
            widget.grid(row=r, column=1, sticky="ew", pady=3)

        self.name = tk.StringVar(value=lad.name if lad else "")
        self.vendor = tk.StringVar(value=lad.vendor if lad else "")
        self.dye = tk.StringVar(value=lad.dye if lad else "ROX")
        self.color = tk.StringVar(value=(lad.color if lad and lad.color
                                         else lib.default_color(self.dye.get())))
        ch = lad.default_channel if lad and lad.default_channel else ""
        self.channel = tk.StringVar(value=str(ch))
        self.exclude = tk.StringVar(value=", ".join(_fmt(b) for b in lad.exclude)
                                    if lad else "")
        self.source = tk.StringVar(value=lad.source if lad else "")
        row(0, "Name", ttk.Entry(frm, textvariable=self.name, width=40))
        row(1, "Vendor", ttk.Entry(frm, textvariable=self.vendor))
        dyebox = ttk.Combobox(frm, textvariable=self.dye, values=COMMON_DYES)
        dyebox.bind("<<ComboboxSelected>>", self._dye_changed)
        row(2, "Dye", dyebox)

        cf = ttk.Frame(frm)
        self.swatch = tk.Label(cf, width=4, relief="groove",
                               bg=self.color.get())
        self.swatch.pack(side=tk.LEFT)
        ttk.Button(cf, text="Choose…", command=self._pick_color).pack(
            side=tk.LEFT, padx=6)
        ttk.Label(cf, text="(display only)", foreground="#666").pack(side=tk.LEFT)
        row(3, "Colour", cf)
        row(4, "Default ladder channel",
            ttk.Combobox(frm, textvariable=self.channel, width=4,
                         values=["", "1", "2", "3", "4"], state="readonly"))

        lf = ttk.Frame(frm)
        self.lengths = tk.Text(lf, height=5, width=48, wrap="word")
        self.lengths.pack(fill=tk.BOTH, expand=True)
        if lad:
            self.lengths.insert("1.0", ", ".join(_fmt(b) for b in lad.lengths))
        ttk.Label(lf, text="Band lengths in bp, separated by commas, spaces "
                  "or new lines.", foreground="#666").pack(anchor="w")
        row(5, "Lengths (bp)", lf)
        row(6, "Exclude from fit", ttk.Entry(frm, textvariable=self.exclude))
        row(7, "Source", ttk.Entry(frm, textvariable=self.source))
        self.notes = tk.Text(frm, height=3, width=48, wrap="word")
        if lad and lad.notes:
            self.notes.insert("1.0", lad.notes)
        row(8, "Notes", self.notes)

        btn = ttk.Frame(frm)
        btn.grid(row=9, column=0, columnspan=2, sticky="e", pady=(8, 0))
        ttk.Button(btn, text="Save", command=self._ok).pack(side=tk.LEFT, padx=4)
        ttk.Button(btn, text="Cancel", command=self.destroy).pack(side=tk.LEFT)
        self.bind("<Escape>", lambda e: self.destroy())
        _modal(self, parent)

    def _dye_changed(self, _e=None):
        self.color.set(lib.default_color(self.dye.get()))
        self.swatch.configure(bg=self.color.get())

    def _pick_color(self):
        rgb, hexc = colorchooser.askcolor(self.color.get(), parent=self)
        if hexc:
            self.color.set(hexc)
            self.swatch.configure(bg=hexc)

    def build(self) -> Ladder:
        """The ladder the form describes (ValueError when it cannot be read)."""
        lengths = lib.parse_lengths_text(self.lengths.get("1.0", "end"))
        excl = lib.parse_lengths_text(self.exclude.get()) \
            if self.exclude.get().strip() else []
        if not lengths:
            raise ValueError("enter the band lengths")
        ch = self.channel.get().strip()
        return Ladder(name=self.name.get().strip(), lengths=tuple(lengths),
                      dye=self.dye.get().strip() or "custom",
                      source=self.source.get().strip(),
                      notes=self.notes.get("1.0", "end").strip(),
                      vendor=self.vendor.get().strip(),
                      color=self.color.get(), exclude=tuple(excl),
                      default_channel=int(ch) if ch else None)

    def _ok(self):
        try:
            lad = self.build()
            problems = lib.validate(lad)
        except ValueError as e:
            problems = [str(e)]
        if problems:
            messagebox.showerror("Ladder", "\n".join(problems), parent=self)
            return
        self.result = lad
        self.destroy()


# --------------------------------------------------------------------------- #
# manager
# --------------------------------------------------------------------------- #
class LadderManager(tk.Toplevel):
    """List the bundled and user ladders; new / import / edit / duplicate /
    delete / export."""

    COLS = ("name", "vendor", "dye", "bands", "range", "kind", "source")

    def __init__(self, parent):
        super().__init__(parent)
        self.title("Size-standard ladders")
        self.geometry("900x460")
        frm = ttk.Frame(self, padding=8)
        frm.pack(fill=tk.BOTH, expand=True)
        self.tree = ttk.Treeview(frm, columns=self.COLS, show="headings",
                                 height=10, selectmode="browse")
        widths = {"name": 190, "vendor": 140, "dye": 70, "bands": 50,
                  "range": 90, "kind": 70, "source": 260}
        for c in self.COLS:
            self.tree.heading(c, text=c.capitalize())
            self.tree.column(c, width=widths[c], anchor="w",
                             stretch=c in ("name", "source"))
        self.tree.pack(fill=tk.BOTH, expand=True)
        self.tree.bind("<<TreeviewSelect>>", self._show)
        self.tree.bind("<Double-1>", lambda e: self.edit())
        self.detail = tk.Text(frm, height=5, wrap="word", relief="flat",
                              background=self.cget("background"))
        self.detail.pack(fill=tk.X, pady=(6, 0))
        btn = ttk.Frame(frm)
        btn.pack(fill=tk.X, pady=(6, 0))
        for text, cmd in (("New…", self.new), ("Import…", self.import_file),
                          ("Edit…", self.edit), ("Duplicate", self.duplicate),
                          ("Delete", self.delete), ("Export…", self.export)):
            ttk.Button(btn, text=text, command=cmd).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn, text="Close", command=self.destroy).pack(side=tk.RIGHT)
        self.bind("<Escape>", lambda e: self.destroy())
        self._rows = {}
        self.refresh()
        _modal(self, parent)

    # ------------------------------------------------------------------ data
    def refresh(self, select: Optional[str] = None):
        self.tree.delete(*self.tree.get_children())
        self._rows = {}
        for row in lib.all_ladders():
            lad, key = row["ladder"], row["key"]
            iid = ("b:" if row["builtin"] else "u:") + key
            self._rows[iid] = row
            self.tree.insert("", tk.END, iid=iid, values=(
                lad.name, lad.vendor, lad.dye, len(lad.lengths),
                f"{_fmt(lad.lengths[0])}–{_fmt(lad.lengths[-1])}",
                "bundled" if row["builtin"] else "mine", lad.source))
        if select and select in self._rows:
            self.tree.selection_set(select)
            self.tree.see(select)
        self._show()

    def _current(self):
        sel = self.tree.selection()
        return self._rows.get(sel[0]) if sel else None

    def _show(self, _e=None):
        row = self._current()
        self.detail.configure(state="normal")
        self.detail.delete("1.0", tk.END)
        if row:
            lad = row["ladder"]
            txt = "Lengths: " + ", ".join(_fmt(b) for b in lad.lengths)
            if lad.exclude:
                txt += "\nExcluded from the fit: " + ", ".join(
                    _fmt(b) for b in lad.exclude)
            if lad.default_channel:
                txt += f"\nDefault ladder channel: {lad.default_channel}"
            if lad.notes:
                txt += "\n" + lad.notes
            self.detail.insert("1.0", txt)
        self.detail.configure(state="disabled")

    def _save(self, lad: Ladder, key=None, overwrite=False) -> Optional[str]:
        try:
            return lib.save_user_ladder(lad, key=key, overwrite=overwrite)
        except (ValueError, OSError) as e:
            messagebox.showerror("Ladder", str(e), parent=self)
            return None

    # --------------------------------------------------------------- actions
    def new(self, ladder: Optional[Ladder] = None, title="New ladder"):
        ed = LadderEditor(self, ladder, title=title)
        self.wait_window(ed)
        if ed.result is not None:
            key = self._save(ed.result)
            if key:
                self.refresh(select="u:" + key)

    def import_file(self):
        path = filedialog.askopenfilename(
            parent=self, title="Import ladder",
            filetypes=[("Ladder files", "*.json *.csv *.tsv *.txt *.tdb"),
                       ("All files", "*.*")])
        if not path:
            return
        try:
            ladders = lib.import_ladder_file(path)
        except (OSError, ValueError) as e:
            messagebox.showerror("Import ladder", str(e), parent=self)
            return
        for lad in ladders:
            # review every import: name it, pick the dye and colour
            from dataclasses import replace
            lad = replace(lad, name=lib.unique_name(lad.name))
            self.new(lad, title=f"Import ladder — check and save")

    def edit(self):
        row = self._current()
        if not row:
            return
        if row["builtin"]:
            if messagebox.askyesno(
                    "Ladder", "Bundled ladders cannot be changed. Make an "
                    "editable copy?", parent=self):
                self.duplicate()
            return
        ed = LadderEditor(self, row["ladder"], title="Edit ladder")
        self.wait_window(ed)
        if ed.result is None:
            return
        new_key = lib.key_of(ed.result.name)
        if new_key != row["key"]:            # renamed: write new, drop old
            key = self._save(ed.result)
            if key:
                lib.delete_user_ladder(row["key"])
        else:
            key = self._save(ed.result, key=row["key"], overwrite=True)
        if key:
            self.refresh(select="u:" + key)

    def duplicate(self):
        row = self._current()
        if not row:
            return
        from dataclasses import replace
        lad = replace(row["ladder"], name=lib.unique_name(row["ladder"].name + " copy"))
        self.new(lad, title="Copy of ladder")

    def delete(self):
        row = self._current()
        if not row:
            return
        if row["builtin"]:
            messagebox.showinfo("Ladder", "Bundled ladders cannot be deleted.",
                                parent=self)
            return
        if messagebox.askyesno("Delete ladder",
                               f"Delete {row['ladder'].name!r} from your library?",
                               parent=self):
            lib.delete_user_ladder(row["key"])
            self.refresh()

    def export(self):
        row = self._current()
        if not row:
            return
        path = filedialog.asksaveasfilename(
            parent=self, defaultextension=".json",
            initialfile=lib.key_of(row["ladder"].name) + ".json",
            filetypes=[("Ladder JSON", "*.json")])
        if path:
            lib.export_ladder(row["ladder"], path)


# --------------------------------------------------------------------------- #
# the sizing set-up
# --------------------------------------------------------------------------- #
class SizingSetup(tk.Toplevel):
    """Pick a ladder and the two physical channels.  ``result`` is
    ``{"key", "ladder", "ladder_channel", "sample_channel"}`` or None."""

    def __init__(self, parent, last: Optional[dict] = None):
        super().__init__(parent)
        self.title("Fragment-length sizing")
        self.result = None
        last = dict(last or {})
        frm = ttk.Frame(self, padding=10)
        frm.pack(fill=tk.BOTH, expand=True)
        ttk.Label(frm, text="Size standard (ladder):").grid(row=0, column=0,
                                                            sticky="w")
        self.choice = tk.StringVar()
        self.box = ttk.Combobox(frm, textvariable=self.choice, state="readonly",
                                width=52)
        self.box.grid(row=1, column=0, columnspan=3, sticky="ew", pady=(2, 6))
        self.box.bind("<<ComboboxSelected>>", self._ladder_changed)
        ttk.Button(frm, text="Manage ladders…", command=self._manage).grid(
            row=1, column=3, padx=(6, 0))
        self.lch = tk.IntVar(value=int(last.get("ladder_channel") or 4))
        self.sch = tk.IntVar(value=int(last.get("sample_channel") or 2))
        ttk.Label(frm, text="Ladder channel (1-4, physical):").grid(
            row=2, column=0, sticky="w")
        ttk.Spinbox(frm, from_=1, to=4, width=4, textvariable=self.lch,
                    state="readonly").grid(row=2, column=1, sticky="w")
        ttk.Label(frm, text="Sample channel (1-4, physical):").grid(
            row=3, column=0, sticky="w")
        ttk.Spinbox(frm, from_=1, to=4, width=4, textvariable=self.sch,
                    state="readonly").grid(row=3, column=1, sticky="w")
        self.hint = ttk.Label(frm, text="", foreground="#666")
        self.hint.grid(row=4, column=0, columnspan=4, sticky="w", pady=(6, 0))
        btn = ttk.Frame(frm)
        btn.grid(row=5, column=0, columnspan=4, sticky="e", pady=(10, 0))
        ttk.Button(btn, text="Size", command=self._ok).pack(side=tk.LEFT, padx=4)
        ttk.Button(btn, text="Cancel", command=self.destroy).pack(side=tk.LEFT)
        self.bind("<Return>", lambda e: self._ok())
        self.bind("<Escape>", lambda e: self.destroy())
        self._fill(select=str(last.get("ladder") or ""))
        _modal(self, parent)

    def _fill(self, select: str = ""):
        self._rows = lib.all_ladders()
        self.box["values"] = [describe(r["ladder"]) for r in self._rows]
        idx = 0
        sel = select.strip().lower()
        for i, r in enumerate(self._rows):
            if sel and (r["key"] == sel or r["ladder"].name.lower() == sel
                        or r["key"] == lib.key_of(sel)):
                idx = i
                break
        if self._rows:
            self.box.current(idx)
            self._update_hint()

    def _row(self):
        i = self.box.current()
        return self._rows[i] if 0 <= i < len(self._rows) else None

    def _update_hint(self):
        row = self._row()
        if row:
            lad = row["ladder"]
            src = f"Source: {lad.source}" if lad.source else ""
            self.hint.configure(text=src[:110])

    def _ladder_changed(self, _e=None):
        row = self._row()
        if row and row["ladder"].default_channel:
            self.lch.set(row["ladder"].default_channel)
        self._update_hint()

    def _manage(self):
        keep = self._row()["key"] if self._row() else ""
        m = LadderManager(self)
        self.wait_window(m)
        self._fill(select=keep)
        _modal(self, self.master)

    def _ok(self):
        row = self._row()
        if row is None:
            return
        if self.lch.get() == self.sch.get():
            messagebox.showerror("Fragment sizing", "The ladder and the sample "
                                 "must be on different channels.", parent=self)
            return
        self.result = {"key": row["key"], "ladder": row["ladder"],
                       "ladder_channel": int(self.lch.get()),
                       "sample_channel": int(self.sch.get())}
        self.destroy()


def ask_sizing_setup(parent, last: Optional[dict] = None) -> Optional[dict]:
    dlg = SizingSetup(parent, last)
    parent.wait_window(dlg)
    return dlg.result


def open_manager(parent) -> None:
    parent.wait_window(LadderManager(parent))
