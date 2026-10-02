# Limoncello CE Analyzer (Python)

A desktop viewer and base caller for capillary-electrophoresis traces,
wired to the plate-validated **best basecaller** configs.

## Features

| Feature | Status |
|---------|--------|
| Multiple data folders | Yes — File → Add data folder |
| Wells list (.rsd / .scf / .ab1 / text) | Yes — multi-select |
| Number of graphs (1–8) | Yes — spinbox |
| Views: raw / processed / base-called / sequencing trace (ESD peaks) / wrap | Yes — View menu |
| Basecaller versions | `mb1000_accuracy`, `mb1000_length`, `mb4000_accuracy`, `mb4000_length`, `pos_bonus07`, `pos_profile`, `hz_soften`, `raw_peaks` |
| Advanced params (knobs) | Base calling → ⚙ Basecall settings… (dialog) |
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
| Undo a mistaken base call (clear overlay) | Yes — Base calling → Clear base calls (undo) |
| Sequencing trace (ESD peaks) view | Yes — Base calling menu + View menu |
| Genotyping — manual peak picking (click/area, +A tagging, CTC-CE duplex internal standard, CSV/Excel/JSON export for ML) | Yes — Genotyping menu |
| Fragments/alleles used for ML training library | Built by manual picks (with or without internal standard) |

## Requirements

```bash
pip install numpy scipy matplotlib pillow   # pillow only for backgrounds
pip install openpyxl                        # optional: Excel .xlsx peak export
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
3. Choose **View** (raw / processed / called / ESD peaks / wrap).
4. Pick a **basecaller version**; rarely-needed tuning lives in **Base calling → ⚙ Basecall settings…** → *Apply + redraw* or *Basecall selected*. Ran the call on fragment/genotyping data by mistake? **Base calling → Clear base calls (undo)**.
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

## ESD-style sequencing view & undoing a call

**Base calling** is the sequencing menu (once *Analysis*). Besides the call
actions it offers the two plot styles so a base call and its result live
together:

- **Sequencing trace (ESD peaks)** — the MegaBACE basecaller's processed
  picture: the four dye traces with one ring per peak, coloured like the
  trace that carries it. When the run is base-called the rings sit at the
  caller's peaks, so a heterozygous position shows as two rings at one scan
  and a mutation as a ring in an unexpected channel. No letters and no
  quality line clutter this view (rings fall back to simple per-channel peak
  detection until the run is called).
- **Plain channels (processed)** — the standard channel view.

By default the app no longer overlays the Q (0–100) quality line after a
call, so the electropherogram stays clean; re-enable it under
**View → Show quality profile (0-100)** if wanted.

If a sequence base call was run on genotyping/fragment data by mistake (the
letters, peak marks and quality curve look confusing), undo it with
**Base calling → Clear base calls (undo)**. Only the in-memory view is reset;
the `.rsd`/`.scf` data files are never modified.

The **Genotyping** menu is an independent top-level heading (between *Base
calling* and *Comments*). It holds **Manual peak picking…**, which needs no
second window and no pop-up: the main window stays exactly as it is — the same
stacked viewer, zoom bars, Reset view, µA overlay and **Channels** row — but
clicking a peak now records it for that well, and the *Called sequence* box
below the plot becomes the *Picked peaks* table. **Peak picking is built from
scratch:**

- Show several wells at once with the **Graphs** spinbox (4–6 at a time is
  the intended workflow), zoom in with the axis bars / mouse wheel, and click
  each peak on its own subplot. Then move the batch onward —
  **Genotyping → Next batch / Previous batch** pages 4–6 wells at a time — and
  pick the next wells. The table below accumulates **every well you visited**,
  so you save the whole run in one go.
- Click a peak (or just beside it) — the best available algorithm locates it,
  shades the peak area and logs scan, channel/base, height (V) and area
  (V·scan). Clicking nearer the peak's actual height separates two channels
  that share a scan. The recognition method is selectable: Best prominence,
  Simple local maxima, or Gaussian fit.
- **A peak's area is measured between its own two valleys**, and the area is
  integrated out to those valleys — a variant fraction is only as good as the
  two areas it is built from, so the peak keeps every tail scan that is really
  its own. The click radius stays generous so the peak is found reliably, but
  the baseline for height, area, onset/end and the variant fraction is only ever
  drawn between the nearest valley *on each side of that peak*. A tight pair of
  alleles therefore no longer borrows each other's area — the small allele of a
  het pair used to report a nearly equal area, and a fake ~50/50 fraction with
  it. The fractions in the table are always `small/(small+large)` of the picked
  areas, never of the raw window.
  - The valley search walks a slightly smoothed copy of the trace (the area
    itself is always integrated from the raw samples), so the broad dip that
    closes a peak is found while the one- or two-scan ripples on its shoulder
    are not mistaken for it, and a peak whose own valley lies far out keeps its
    full width instead of being cut at a fixed distance.
  - The one case the peak's own valley cannot settle is a **low minor allele
    sitting on a main peak's tail**, which has no dip of its own. There the
    area is still bounded at 1.5 × the run's peak spacing — the same reach as
    the heterozygote window — so it can never measure into where a neighbour
    could start, at the cost of reading that minor allele's fraction slightly
    low. The fractions in the table are always `small/(small+large)` of the
    picked areas, never of the raw window.
- **Hovering a subplot shows a crosshair and a live readout** in that
  subplot's corner: the scan under the pointer, its voltage, and the dominant
  peak within one base on any visible channel, marked `[picked]` once you have
  it — so you can line up on the exact hump before clicking. The reported peak
  is the tallest nearby one, not merely the closest, so a flat-channel wobble
  cannot win over the real peak.
- **The picked-peaks table is sorted by scan** (then channel) inside each well,
  so the order you clicked in never scrambles the listing; **Undo last pick**
  still removes the most recent click.
- **Channel identity follows the run's dye order, everywhere.** `Base order
  (instrument)` (ACTG, TGCA, GATC, CTAG) decides which dye sits on Ch1–Ch4;
  the channel checkboxes, both legends and the exported `Ch` column all use
  that same mapping, so hiding Ch1 hides the same trace in the viewer and while
  picking. Changing the dye order with picks already on screen keeps them and
  says so in the status bar.
- **Mark start/end from the 2nd derivative** (on by default): square ticks
  show where each picked peak lifts off its baseline. The 2nd derivative of
  the smoothed trace crosses the noise floor from flat to concave-up at the
  true start, and concave-up again on the return at the end. Uncheck for a
  clean look. The start/end scans are also written to the exported table.
- **Add +A** (on by default) tags the strongest satellite *trailing* the main
  peak — the Taq A-addition a few scans later. Nothing is tagged in front of
  the main peak: for a single-base-extension product that leading shoulder is
  another A-addition on the GC-clamp side, not stutter, so it is left unmarked
  rather than mislabelled. Turn the option off when clicking allele peaks so a
  second allele is not swallowed by the +A tag.
- **Internal standard** (CTC-CE duplex pattern): the four standard peaks are
  all **one** fragment (the same number of base pairs) — cycling-temperature
  capillary electrophoresis separates them by *sequence*: peaks 1–2 are the two
  **homoduplexes** (they differ by the single SNP base of the rs number), peaks
  3–4 the two **heteroduplexes** made in the PCR when Watson and Crick strands
  pair wrongly, leaving one mismatch base pair. Pick that well's four standard
  main peaks *last* and use **Genotyping → Mark peaks as standard…** (a shared
  fragment length in bp is optional): earlier ones become `HOM1`/`HOM2`, later
  ones `HET1`/`HET2`. Variant ratios come from the **relative areas** of these
  duplex peaks, so no bp ladder is involved.
- **Mutant/variant fraction** is computed for two main peaks of the *same
  channel* within one repeat (~one base): small/(small+large), shown in the
  table. The grouping window scales with the run's own peak spacing.
- **Mass-action MF** (CTCE) — **Genotyping → Tag duplex species for MF…** tags
  *one allelic position* (the last-picked main's own cluster on its channel, so
  two positions in one well stay apart) with the same `HOM1`/`HOM2`/`HET1`/`HET2`
  names the internal standard gets, and that position then reports the PCR
  mass-action mutant fraction in the **MF** column:

  ```
  MF = (A_MUT + ½ × A_HET) / (A_WT + A_MUT + A_HET)
  ```

  where `A_WT`/`A_MUT` are the homoduplex areas and `A_HET` the combined
  heteroduplex area. The ½ term is the point of the formula: a clean
  heterozygote reads **0.5**, not the 0.25 a plain area ratio of the two
  homoduplexes gives, and below ≈5 % MF — where the mutant strands have all
  re-annealed and *no mutant homoduplex is visible at all* — the whole low
  fraction is carried by `A_HET`. The **ai** column adds the allelic imbalance
  `A_HOMO1 / (A_HOMO1 + A_HOMO2)`, which needs no wild-type choice. Both are
  blank until a position is tagged; `fraction` keeps its plain meaning. The
  peak count decides the homoduplex/heteroduplex split: 4 = 2+2, 3 = 1+2, 2 = 2+0.
- **Undo last pick / Clear picks** manage the picks (`Ctrl+Z` undoes the last
  pick while peak picking is active); an already-picked **area
  cannot be picked again** (the click is refused with a status-bar message —
  undo it first to re-pick). Neighbouring peaks such as the two alleles of a
  heterozygote remain pickable; **Save peaks table…**
  writes CSV (Excel-ready, UTF-8 BOM), Excel `.xlsx` or JSON: file, well, scan,
  channel, base, kind (`main`/`+A`), start/end scan (2nd derivative),
  height, area, duplex label (`HOM1`/`HOM2`/`HET1`/`HET2`), length (bp),
  fraction, and the mass-action `mf` and `ai`. Rows are written **grouped by
  sample** (run folder, then well name)
  and **in scan order inside each sample**, never in click order, so one
  sample's peaks are never interleaved with another's. That table is a
  labelled **training library for ML** — picking
  works with or without an internal standard.
- **Exit peak picking** returns to the normal trace viewer (the same menu item
  toggles back and forth).

## Comments & run info

The **Comments** menu sits between *Genotyping* and *Help*:

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
