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
| Basecaller versions | `pos_bonus07`, `pos_profile`, `hz_soften`, `raw_peaks` |
| Advanced params (knobs) | Analysis → ⚙ Basecall settings… (dialog) |
| Wrap rows, tour interval | View menu |
| Per-graph label (folder/file) | Yes — drawn in each plot's top-left |
| X axis | Tick numbers only on the bottom pane; **Time (min)** toggle (1.75 Hz) |
| Pan/zoom | Axis bars: drag to pan, wheel to zoom; right-click/Home resets |
| Random background pictures | Yes — `Background*.jpg` / `BG*.jpg` beside the app |
| Export FASTA / CSV / trace text | Yes |
| Save graph image (PNG/PDF/SVG) | Yes — File → Save graph image… |
| Save/load settings JSON | Yes |
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

| Name | Role |
|------|------|
| **pos_bonus07** | Recommended dual-aware (default) |
| **pos_profile** | Max matched_bp |
| **hz_soften** | Mild mid-zone less deconv |
| **raw_peaks** | Envelope peaks only (debug) |

Dye order default **TGCA** (ET plate).

## Note on environment

This sandbox may not have `tkinter`; run the app on your **local** Windows/macOS/Linux
desktop.
