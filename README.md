# Limoncello CE Analyzer (Python)

A desktop viewer and base caller for capillary-electrophoresis traces,
wired to the plate-validated **best basecaller** configs.

## Features

| Feature | Status |
|---------|--------|
| Multiple data folders | Yes — File → Add data folder |
| Wells list (.rsd / .scf / .ab1 / text) | Yes — multi-select |
| Number of graphs (1–8) | Yes — spinbox |
| Views: raw / processed / base-called / wrap | Yes — View menu |
| Basecaller versions | `mb1000_accuracy`, `mb1000_length`, `mb4000_accuracy`, `mb4000_length`, `pos_bonus07`, `pos_profile`, `hz_soften`, `raw_peaks` |
| Advanced params (knobs) | Analysis → ⚙ Basecall settings… (dialog) |
| Wrap rows, tour interval | View menu |
| Per-graph label (folder/file) | Yes — drawn in each plot's top-left |
| X axis | Tick numbers only on the bottom pane; **Time (min)** toggle (1.75 Hz) |
| Trace colors / dye sets | Combobox beside the channel toggles + View → Trace colors palette |
| Pan/zoom | Axis bars: drag to pan, wheel to zoom; right-click/Home resets |
| Random background pictures | Yes — `Background*.jpg` / `BG*.jpg` beside the app |
| Export FASTA / CSV / trace text | Yes |
| Save graph image (PNG/PDF/SVG) | Yes — File → Save graph image… |
| Save/load settings JSON | Yes |
| Run comments (per-file, Stored beside the data) | Yes — Comments menu |
| Run info panel | Yes — Comments menu |
| Genotyping | Planned — next to base calling, later release |

## Requirements

```bash
pip install numpy scipy matplotlib pillow   # pillow only for backgrounds
# tkinter: included on Windows/macOS Python installers
# Ubuntu/Debian: sudo apt install python3-tk
```

The basecaller package is **bundled** in the zip (see the folder tree below),
so `.rsd`/`.scf` files work on any machine out of the box. `.ab1` and text
traces don't need it at all.

```
Limoncello scorer.zip
├── Limoncello scorer/
│   ├── sequence_analyzer.py
│   ├── analyzer_core.py
│   ├── README.md
│   ├── requirements.txt
│   ├── Background.jpg            # empty-view pictures (bundled)
│   ├── BG2.jpg, BG3.jpg, BG4.jpg
│   └── example_data/M13/         # 8 M13 wells (A01–A08.rsd)
└── BEST_BASECALLER_RELEASE/
    └── cimarron_basecaller/      # base caller (imported by the app)
```

## Run

```bash
cd "Limoncello scorer"
python3 sequence_analyzer.py

# Or open with a folder already loaded
python3 sequence_analyzer.py --folder "/path/to/Limoncello scorer/example_data/M13"
```

To try it immediately: **File → Add data folder** →
`example_data/M13`. These are `.rsd` files, so the basecaller package is
required to view them (place `BEST_BASECALLER_RELEASE` beside this folder).

## Typical workflow

1. **File → Add data folder** → select a directory of `.rsd` files (add more folders if needed).
2. Select one or more wells in the list; set **Graphs to show**.
3. Choose **View** (raw / processed / called / wrap).
4. Pick a **basecaller version**; rarely-needed tuning lives in **Analysis → ⚙ Basecall settings…** → *Apply + redraw* or *Basecall selected*.
5. Inspect sequence pane; **File → Export sequence (FASTA)** for BLAST.

## Navigating the plot

- The bottom bar is the **X (scan)** axis, the right bar is the **Y (signal)** axis.
  Grab a bar's thumb and slide to pan that axis.
- X tick numbers appear only under the **bottom pane**, so stacked plots keep
  their height. Tick **Time (min)** to read the axis as minutes:seconds instead
  of scans (1.75 scans/s, ~0.57 s/scan).
- Roll the mouse wheel over a bar to zoom that axis; **double-click** a bar to
  reset it. Both bars drive every visible plot together.
- Lost in the zoom? **Right-click** the plot or a bar, press **Home**, or use the
  **⟲ Reset view** button.
- With no well selected, the empty plot area shows a randomly chosen picture
  from `Background*.jpg` / `BG*.jpg` in this folder (a new one each redraw).
- **File → Save graph image…** writes the current plot area to PNG/PDF/SVG.
  The old matplotlib toolbar was removed in favour of these controls.

## Keyboard shortcuts

| Key | Action |
|-----|--------|
| ↑ / ↓ / PgUp / PgDn | page through wells |
| Space | start / stop auto-tour |
| Delete | remove selected sample(s) |
| Home / right-click | reset the X & Y view |
| F1 | user manual |
| F11 | full screen / restore |
| Esc | cancel a long job, or leave full screen |

## Basecaller versions

`pos_bonus07` and `pos_profile` are now aliases of the unified
**instrument × mode** presets (`*_accuracy` / `*_length`). The list lives in
`analyzer_core.BASECALLER_VERSIONS`; the underlying parameters come from the
bundled `BEST_BASECALLER_RELEASE/configs.py`.

| Name | Role |
|------|------|
| **mb1000_accuracy** | MegaBACE 1000 · longest error-free run, high %ID (deep Q-trim). Default. |
| **mb1000_length** | MegaBACE 1000 · longest read at ID ≥95% (mid hard-zone) |
| **mb4000_accuracy** | MegaBACE 4000 · error-free mode + MB4000 spectral CHM |
| **mb4000_length** | MegaBACE 4000 · length mode + MB4000 spectral CHM |
| **pos_bonus07** | legacy alias of `mb1000_accuracy` (kept for saved settings) |
| **pos_profile** | legacy alias of `mb1000_length` (kept for saved settings) |
| **hz_soften** | accuracy base + mild mid hard-zone |
| **raw_peaks** | Envelope peaks only (debug) |

> Use `mb4000_*` only on MegaBACE 4000 traces — on 1000 data the 4000
> spectral matrix is the wrong chemistry and reads collapse.

Base order default **ACTG** (MegaBACE: Ch1=A, Ch2=C, Ch3=T, Ch4=G).
DYEnamic runs are **TGCA** — set it under ⚙ Basecall settings → Dye/channel.

## Trace colors / MegaBACE dye sets

The palette combo box (top of the plot area, near the **Channels:** toggles)
and **View → Trace colors** control the trace colors. Traces are drawn in the
run's **physical channel order** (Ch1..Ch4) — the first tick is always the
first channel — each labelled by the base it carries (that base order is the
chemistry under **⚙ Basecall settings → Dye/channel → Base order**, default
**ACTG**).

**Classic (default)** is the MegaBACE software look — colors **fixed per base
letter** (A green, C blue, T red, G black). With the default base order that
shows Ch1=A green, Ch2=C blue, Ch3=T red, Ch4=G black; other chemistries only
shift which channel shows which color (letters keep theirs):

| Option | Ch1 | Ch2 | Ch3 | Ch4 |
|--------|-----|-----|-----|-----|
| **Classic** *(default, base order ACTG)* | A green | C blue | T red | G black |

Sequencing dye sets are fixed per **base**, as the MegaBACE sequence analyser
shows them (A green, C blue, T red, G black) — the same traces as Classic, but
each kit shown in its own channel order:

| Option | Ch1 | Ch2 | Ch3 | Ch4 |
|--------|-----|-----|-----|-----|
| **Seq DYEnamic (T·G·C·A)** | T red | G black | C blue | A green |
| **Seq ET primer (A·C·T·G)** | A green | C blue | T red | G black |

Genotyping dye sets are fixed **per channel** (no base colors):

| Option | Ch1 | Ch2 | Ch3 | Ch4 |
|--------|-----|-----|-----|-----|
| **Genotyping (R·B·Blk·G)** — dye set 1 (ET-ROX, FAM, HEX/NED, TET) | Red | Blue | Black | Green |
| **Genotyping (G·B·R·Blk)** — Green/Blue/Red/Black order Fragment-readers are used to | Green | Blue | Red | Black |

Changing the scheme recolors the traces only; it does not change the base
labels or the basecalling dye order.

## Comments & run info

The **Comments** menu sits between *Analysis* and *Help*:

- **Run comments…** writes a note for the selected run (main-curve style). It is
  stored beside the data file as `<file>.comment.txt`, so it travels with the
  run without ever touching the binary `.rsd`/`.scf` header. A saved comment is
  shown in the sequence pane under the `>` header line.
- **Run info…** gives a read-only rundown of the selected well: source, scan
  count and run time, base order, basecaller preset, sequence statistics
  (length, Qmean/Qmin, N, peak spacing), per-channel signal maxima and the
  instrument current.

## Note on environment

This sandbox may not have `tkinter`; run the app on your **local** Windows/macOS/Linux
desktop.
