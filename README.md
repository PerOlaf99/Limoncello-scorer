# Limoncello CE Analyzer (Python)

A desktop viewer and base caller for capillary-electrophoresis traces,
wired to the plate-validated **best basecaller** configs.

The full user manual, with quick starts for each job, is in the app:
**Help → User manual…** (F1), with contents and search.

## Features

| Feature | Status |
|---------|--------|
| Multiple data folders | Yes — File → Add data folder |
| Wells list (.rsd / .scf / .ab1 / .fsa / text) | Yes — multi-select |
| Number of graphs (1–12) | Yes — Graphs box in the view bar, with a page indicator |
| Overlay runs of the same well | Yes — view bar / View menu (Run01 solid, Run02 dashed, …) |
| Views: raw / processed / base-called / sequencing trace (ESD peaks) / wrap | Yes — View menu |
| Basecaller versions | `mb1000_accuracy`, `mb1000_length`, `mb4000_accuracy`, `mb4000_length`, `pos_bonus07`, `pos_profile`, `hz_soften`, `raw_peaks` |
| Advanced params (knobs) | Base calling → ⚙ Basecall settings… (dialog) |
| Wrap rows, tour interval | View menu |
| Per-graph label | Well in bold + run name, in each plot's top-left |
| Plate map (96 / 384 wells, coloured by call; click to show wells) | Yes — ▦ Plate in the view bar / View → Plate map… |
| Dye order per run (declared by the file, or chosen per run) | Yes — Comments → Dye order of this run… |
| Picking-mode banner (what the mouse does, Exit / Esc) | Yes |
| X axis | Tick numbers only on the bottom pane; **Time (min)** toggle (1.75 Hz) |
| Trace colors / dye sets | Combobox beside the channel toggles + View → Trace colors palette |
| Pan/zoom | Wheel on the graphs zooms at the cursor (Ctrl: Y, Shift: pan), drag zooms a box, middle-drag pans; Zoom −/+ and a scan-range box; right-click/double-click/Home resets |
| Y scaling | Shared or each graph its own; Fit Y to view |
| Align graphs to each other | Shift+drag a graph, Alt+←/→ nudge, ⇄ Align (cross-correlation in view, or on the internal standard with stretch); display only — true scans kept |
| Random background pictures | Yes — `Background*.jpg` / `BG*.jpg` beside the app |
| Export FASTA / CSV / trace text | Yes |
| Save graph image (PNG/PDF/SVG) | Yes — File → Save graph image… |
| Save/load settings JSON | Yes |
| Run comments (per-file, Stored beside the data) | Yes — Comments menu |
| Run info panel | Yes — Comments menu |
| Undo a mistaken base call (clear overlay) | Yes — Base calling → Clear base calls (undo) |
| Sequencing trace (ESD peaks) view | Yes — Base calling menu + View menu |
| Genotyping — manual peak picking (click/area, +A tagging, CTC-CE duplex internal standard, CSV/Excel/JSON export for ML) | Yes — Manual genotyping menu |
| Fragments/alleles used for ML training library | Built by manual picks (with or without internal standard) |
| Mark library — durable internal-standard marks, engine calls and fitted standards (`library/`) | Yes — see [MARK_LIBRARY.md](MARK_LIBRARY.md) |
| Scoring & QC — genotype calls + confidence, per-well read QC, HTML report | Yes — `scorer.py` (headless CLI `limoncello-scorer`) |
| Batch auto-genotyping — find the internal-standard quartet and call every well | Yes — Auto genotyping menu |
| CTCE analysis — peaks of interest from how they move between run temperatures (or by shape at one temperature, or from picked peaks); Hom2/Hom3/heteroduplex blocks (main peak + stutter) found as a plate-learnt pattern, clipped peaks measured in another channel, minor-variant fraction by mass action | Yes — Auto genotyping → CTCE analysis… (`ctce.py`) |
| Fragment-length sizing — size PCR/ssDNA peaks against a ladder (bp), user-definable ladders | Yes — `scorer.py size` |

## Requirements

```bash
# Core runtime (base caller + scorer): numpy and scipy only
pip install numpy scipy

# Desktop viewer (the GUI) also needs matplotlib; tkinter ships with
# Windows/macOS Python. On Ubuntu/Debian: sudo apt install python3-tk
pip install matplotlib pillow      # pillow only for empty-view backgrounds

# Optional file formats
pip install openpyxl               # Excel .xlsx peak export
pip install biopython              # .scf reading
```

Or install the package with extras (from a clone):

```bash
pip install -e ".[gui]"    # viewer + base caller + scorer
pip install -e ".[all]"    # everything (gui + xlsx + scf)
```

The basecaller package is **vendored in this repository** (`cimarron_basecaller/`
plus `configs.py` at the repo root), so `.rsd`/`.scf` files work on any machine
out of the box. `.ab1` and text traces don't need it at all. If you keep the
caller elsewhere, point `BEST_BASECALLER_RELEASE` at the folder that contains
`cimarron_basecaller/`; `python scorer.py check` tells you what was found.

```
Limoncello-scorer/
├── sequence_analyzer.py       # Tk desktop viewer / basecaller GUI
├── analyzer_core.py           # loaders + basecaller wrapper (vendor discovery, --check)
├── genotyping.py              # headless PeakPicker engine + Tk picker dialog
├── fragment_sizing.py         # ladder registry + fragment-length sizing (headless)
├── ladder_library.py          # ladder library: bundled + user ladders, import (headless)
├── ladder_gui.py              # ladder manager, editor and sizing set-up windows (Tk)
├── scorer.py                  # genotype calls, read QC, fragment sizing, HTML report
├── basecall.py                # standalone FASTA/FASTQ basecaller CLI
├── configs.py                 # basecaller presets (GOLDEN / PRECISION / SOFT)
├── cimarron_basecaller/       # vendored base caller package (numpy/scipy)
├── pyproject.toml             # packaging + console scripts
├── LICENSE                    # MIT
├── tests/                     # pytest suite (loaders, scorer, genotyping)
├── Background.jpg             # empty-view pictures (bundled)
├── BG2.jpg, BG3.jpg, BG4.jpg
└── example_data/M13/          # 8 M13 wells (A01–A08.rsd)
```

## Run

```bash
python3 sequence_analyzer.py

# Or open with a folder already loaded
python3 sequence_analyzer.py --folder example_data/M13

# Check that the basecaller, tkinter and optional bits are present, then exit
python3 sequence_analyzer.py --check
```

To try it immediately: **File → Add data folder** →
`example_data/M13`. These are `.rsd` files, handled by the vendored basecaller.

## Typical workflow

1. **File → Add data folder** → select a directory of `.rsd` files (add more folders if needed).
2. Select one or more wells in the list; set **Graphs** in the view bar above the plots.
3. Choose **View** (raw / processed / called / ESD peaks / wrap).
4. Pick a **basecaller version**; rarely-needed tuning lives in **Base calling → ⚙ Basecall settings…** → *Apply + redraw* or *Basecall selected*. Ran the call on fragment/genotyping data by mistake? **Base calling → Clear base calls (undo)**.
5. Inspect sequence pane; **File → Export sequence (FASTA)** for BLAST.

## Navigating the plot

- **On the graphs:** the mouse wheel zooms X around the mouse pointer,
  **Ctrl+wheel** zooms Y and **Shift+wheel** pans. **Drag** to zoom to a range
  (a flat drag zooms X only, a box zooms X and Y); while picking peaks or
  measuring areas use **Ctrl+drag**, since the left button is the tool there.
  **Middle-drag** pans; **double-click**, **right-click** or **Home** resets.
- **View bar** (under the channel toggles): **Graphs** (1–12) with a page
  indicator, **Zoom −/+**, a **Scans from–to** box (type a range, Enter; it
  always shows the scans on screen), **Y scale** *Shared* / *Each graph*,
  **Fit Y to view** and **Overlay runs**, plus **⟲ Reset view**.
- **Y scale** *Shared* keeps one height scale so wells compare; *Each graph*
  scales every graph to its own peaks so a weak well is readable. **Fit Y to
  view** rescales the height to the scans on screen as you zoom.
- **Overlay runs** draws one graph per well with every run of that plate on top
  (Run01 solid, Run02 dashed, Run03 dotted); run folders of one plate are those
  whose names differ only by a trailing `RunNN`.
- The bars still work: bottom = X (scan), right = Y (signal); slide the thumb to
  pan, wheel over a bar to zoom that axis, double-click a bar to reset it. Every
  visible plot moves together.
- X tick numbers appear only under the **bottom pane**; tick **Time (min)** to
  read the axis as minutes:seconds (1.75 scans/s).
- The table under the graphs (genotyping calls, picked peaks or the called
  sequence) sits on a **draggable divider**; it is hidden when the wells on
  screen have nothing for it (**View → Always show table under graphs** keeps it).
- With no well selected, the empty plot area shows a randomly chosen picture
  from `Background*.jpg` / `BG*.jpg` in this folder (a new one each redraw).
- **File → Save graph image…** writes the current plot area to PNG/PDF/SVG.
  The old matplotlib toolbar was removed in favour of these controls.

### Aligning graphs

Slide one graph left or right so a feature lines up with the same feature in
the graphs above or below it — two runs of one plate, wells that migrated a
little faster, or a peak you want to compare across wells.

| How | What it does |
|-----|--------------|
| **Shift + drag** on a graph | slides that graph, live; the others stay put (works in every mode, also while picking) |
| **Alt + ← / →** | moves the graph under the mouse 1 scan (**Alt + Shift**: 10) |
| **Alt + Home** | resets the graph under the mouse |
| **⇄ Align → Align shown graphs to the top one** (view bar, or **View → Alignment**) | cross-correlates the scans **in view** and slides every graph onto the top one; with *Overlay runs* each overlaid run is also aligned onto its graph's first run |
| **⇄ Align → Align on internal standard** | maps each well's marked IS bands onto the top marked well's — an offset, plus a stretch when three or more bands agree on one |
| **⇄ Align → Reset this graph / Reset all alignments** | back to the unshifted view |

Things to know:

- **Alignment is display only.** Each graph keeps its true scan numbers;
  alignment changes *which* scans it shows (internally a per-file map
  `reference = a · scan + b`). Picked peaks, measured areas, IS marks, the
  genotyping table and every export stay in true scans, so you can align and
  then pick as usual.
- **Zoom onto the feature first** before *Align shown graphs to the top one*:
  it lines up whatever dominates the visible window (with the whole trace in
  view that is usually the large primer peak).
- A shifted graph says so in its label, e.g. **⇄ +23 scans**, and **×1.004**
  when it is also stretched.
- The tick numbers and the **Scans** boxes are the **bottom graph's own**
  scans.
- The shared zoom, scroll bars and mouse wheel move all graphs together,
  each keeping its own offset.
- Alignments are **saved in the mark library** per run and well (field
  `align`, see [MARK_LIBRARY.md](MARK_LIBRARY.md)) and come back when the
  plate is opened again; a reset removes them.

## Keyboard and mouse

| Key / mouse | Action |
|-----|--------|
| ↓ ↑ ← → PgDn PgUp | page through wells (whole wells with *Overlay runs*) |
| Space | start / stop the auto-tour (in the Samples list) |
| Delete | hide the selected files from the list |
| Home, double-click, right-click | reset the view |
| Ctrl+Z | undo the last pick or measurement |
| R + click | remove that picked peak |
| Enter in the *Scans* boxes | zoom to the typed scan range |
| wheel / Ctrl+wheel / Shift+wheel | zoom across at the cursor / zoom height / pan |
| drag (Ctrl+drag while picking) | zoom to a range |
| Shift+drag, Alt+← / Alt+→ (Alt+Shift ×10), Alt+Home | slide one graph to align it / reset it |
| middle-drag | pan |
| F1 | user manual (contents + search) |
| F11 | full screen / restore |
| Esc | cancel a long job, leave full screen, or leave peak picking |

## Basecaller versions

The settings dialog (Base calling → ⚙ Basecall settings…) offers four
callers by plain name; the internal preset names stay valid in saved settings
and on the command line (`sequence_analyzer.CALLER_CHOICES` / `CALLER_ALIASES`;
parameters in `configs.py`).

| Shown as | Preset | Role |
|---|---|---|
| **Accuracy (recommended)** | `mb1000_accuracy` | Fewest errors: longest error-free run, deep quality trim. Default. |
| **Length** | `mb1000_length` | Longest read at ≥ 95 % identity |
| **Length, gentler sharpening** | `hz_soften` | As Length, milder sharpening mid-read |
| **Quick look (no processing)** | `raw_peaks` | Tallest channel at each envelope peak (checking only) |

`mb4000_accuracy` / `mb4000_length` are the same callers (no MegaBACE 4000
colour matrix is bundled yet) and `pos_bonus07` / `pos_profile` are older
names; all four load as the matching caller. The per-base "quality" is the
called peak's height scaled 0–100, not an error probability.

> The `mb4000_*` names are kept for MegaBACE 4000 data, but until a 4000
> spectral (colour-separation) matrix is bundled they run the MB1000 settings.

**Dye order is per run.** A run's own files declare it when they can
(sequencing `.rsd` files name the base on each channel, e.g. TGCA for the
MB1000 M13 runs) and that is used automatically. Runs that declare nothing
(genotyping runs) use the general setting under ⚙ Basecall settings →
Dye / channel — default **ACTG** (MegaBACE: Ch1=A, Ch2=C, Ch3=T, Ch4=G);
DYEnamic runs are **TGCA**. **Comments → Dye order of this run…** overrides it
for one run, and **Comments → Run info…** says which order a well uses and why.

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
- **Processed traces (colour-separated)** — the standard channel view.

By default the app no longer overlays the Q (0–100) quality line after a
call, so the electropherogram stays clean; re-enable it under
**View → Show quality profile (0-100)** if wanted.

If a sequence base call was run on genotyping/fragment data by mistake (the
letters, peak marks and quality curve look confusing), undo it with
**Base calling → Clear base calls (undo)**. Only the in-memory view is reset;
the `.rsd`/`.scf` data files are never modified.

The **Manual genotyping** menu (between *Base calling* and *Auto
genotyping*) holds two modes, **Manual peak picking…** and
**Peak area by drag…**, neither of which needs a second window or a pop-up:
the main window stays exactly as it is — the same stacked viewer, zoom bars,
Reset view, µA overlay and **Channels** row — but clicking a peak now records
it for that well, and the *Called sequence* box below the plot becomes the
*Picked peaks* or *Measured areas* table. Only one is active at a time;
choosing the active mode again, or **Manual genotyping → Exit peak picking**,
returns to plain viewing. **Peak picking is built from scratch:**

- Show several wells at once with **Graphs** in the view bar (4–6 at a time
  works well), zoom in with the mouse wheel or a Ctrl+drag box, and click
  each peak on its own subplot. Then move the batch onward —
  **Manual genotyping → Next batch / Previous batch** (or ↓ / ↑) page by the
  Graphs count — and
  pick the next wells. The table below accumulates **every well you visited**,
  so you save the whole run in one go.
- Click a peak (or just beside it) — the best available algorithm locates it,
  shades the peak area and logs scan, channel/base, height (V) and area
  (V·scan). Clicking nearer the peak's actual height separates two channels
  that share a scan.
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
- **Mark start/end from the 2nd derivative** (off by default): square ticks
  show where each picked peak lifts off its baseline. The 2nd derivative of
  the smoothed trace crosses the noise floor from flat to concave-up at the
  true start, and concave-up again on the return at the end. Uncheck for a
  clean look. The start/end scans are also written to the exported table.
- **Add +A** (off by default) tags the strongest satellite *trailing* the main
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
  main peaks *last* and use **Manual genotyping → Mark peaks as standard…** (a shared
  fragment length in bp is optional): earlier ones become `HOM1`/`HOM2`, later
  ones `HET1`/`HET2`. Variant ratios come from the **relative areas** of these
  duplex peaks, so no bp ladder is involved.
- **Mutant/variant fraction** is computed for two main peaks of the *same
  channel* within one repeat (~one base): small/(small+large), shown in the
  table. The grouping window scales with the run's own peak spacing.
- **Mass-action MF** (CTCE) — **Manual genotyping → Tag duplex species for MF…** tags
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
  heterozygote remain pickable. **Peak area by drag…** instead measures the
  area above the straight line between two points you drag, in the same window
  and with the same mode plumbing — handy when a peak's detected bounds look
  wrong and you would rather bracket it by eye than argue with the algorithm.
  Undo / Clear / Save work the same way in either mode. **Save peaks table…**
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

## Data folders and the Samples list

*Add…* (File ▸ Add data folder…) lists every run found under the chosen folder
in **Data folders**, but puts no files in the **Samples** list yet. Double-click
a run — or Ctrl/Shift-click several and press Enter / **▼ List files** — to list
its files; opened runs are ticked and shaded, and double-clicking one again
takes its files out. The broom (**Clear list**) empties the Samples list and
keeps the runs, so a large project can be worked through a pair of runs at a
time. Folders given on the command line (`--folder`) are listed at once.

## Comments & run info

The **Comments** menu sits between *Genotyping* and *Help*:

- **Run comments…** shows the comment typed on the instrument when the run was
  started (the `.rsd` footer's COMMENT, line breaks kept) and lets you keep
  your own note for the run, stored in the run folder as
  `limoncello-run-note.txt` so it travels with the run without ever touching
  the binary files (older per-file `<file>.comment.txt` notes are still read).
  Both are shown in the sequence pane under the `>` header line and in Run info.
- **Run info…** gives a read-only rundown of the selected well: source, scan
  count and run time, base order, basecaller preset, sequence statistics
  (length, Qmean/Qmin, N, peak spacing), per-channel signal maxima and the
  instrument current. For `.rsd` it also decodes the run fields written in the
  instrument footer — run/plate ID, sample and well name, machine ID,
  application, the base caller used, method name/comment, laser mode,
  beamsplitter/filter part numbers and the dye per channel — so two
  similar-looking plates can be told apart from the file itself.

### Run conditions (voltage, time, temperature, PMT)

These settings **are stored in the `.rsd`**, in a fixed-layout binary block that
follows the text footer. Its first five `float32` values are, in order, the run
voltage, run time, injection time, injection voltage and temperature; a 4-byte
tag follows, then the two PMT voltages. No text export is needed.

```
Run Voltage       9        Run time          75 s
Injection time   20 s      Injection voltage 10 kV
Temperature      56 C      PMT Voltage1/2    750 / 750
```

The block is *located*, not read at a hard-coded offset, and a candidate is
accepted only when all seven values fall in a range an electrophoresis run can
actually take. That check matters: trace data can contain any bit pattern, so
an unvalidated fixed-offset read eventually reports nonsense on some well. A
scan of all 27,168 `.rsd` files in the sample runs (OY, SIDS, MB4000 demo)
yields a value for every file and a consistent set within each of the 283 runs.

If the instrument also left a raw-scan text export beside the trace
(`Text/<well>.txt`), a few fields the binary does not carry are merged in — the
plate name, instrument ID, `Number of lines` and the authoritative per-run
`Base order` (dye per channel):

```
Run Info for the file : D:\Data\plate01\A01.rsd
Run Voltage : 9
Comment : Hel plate colo 829 DNA
Grad 50_65C
Inject 10KV, 60 sec, run 9kv, CTCE (53-50)x20
```

Run info shows a **Run conditions** block with the headline settings first.
Exports are optional (most OY runs have none), so those extras simply do not
appear when the file is absent.

Parsing notes: the encoding is detected (the files are UTF-8 in practice, but
UTF-16 is handled), the key/value split uses the **first** colon only so that
Windows paths survive, and fields such as `Comment` may span several lines.

Most `OY/` runs on this system have no `Text/` export — only the `mt_nucl`
runs do. When none is found, the block says so explicitly rather than showing
empty values; re-export the run from the instrument software to populate it.

## CTCE analysis (`ctce.py`)

Cycling temperature CE separates the duplexes of a fragment by how they melt
(Ekstrøm, Bjørheim & Thilly, BMC Genetics 2007). **Auto genotyping → CTCE
analysis…** works on the runs of one plate (folders that differ only by
`RunNN`), in ONE channel, the sample channel:

1. **What melts** — between runs at different temperatures (read from each
   `.rsd`, so a hotter run done first still shifts the right way) the duplex
   peaks move later, each by its own amount; ssDNA and primers hold still.
   At one temperature the shape decides (a duplex has a sharp front and a
   tail), and picked peaks can serve as a template.
2. **Species blocks** — each duplex species is its main peak + shoulder +
   Taq stutter (~15 scans later), measured as one block (10 scans before the
   main peak to 30 after) on a near-flat baseline. The plate learns its own
   pattern as migration-time ratios: Hom2 the largest species, Hom3 beyond
   Hom2's cluster at the ratio the strong wells share, then up to two
   heteroduplexes. Peaks seen static between temperatures, or at a fixed
   ratio to the early peak in most wells without a stutter (ssDNA), are never
   a species. The measuring temperature is the one where most wells resolve
   all four.
3. **Clipped peaks** — where the sample channel hits the 16-bit limit, the
   clipped stretch is rebuilt from the strongest other channel carrying the
   same dye's crosstalk, scaled by the gain fitted on the unclipped part of
   the same well (flagged in the table).
4. **Minor-variant fraction** by mass action, p = (A_Hom3 + ½ A_het) / A_total,
   with ± SD over reruns at the measuring temperature.

Checked against hand-integrated areas (mtDNA #19_20, 96 wells, 55 °C): Hom2
and Hom3 areas within 1–3 % (median), minor fraction within 0.03 of the hand
value in 85 % of wells (mean |error| 0.020). On fragment #8 the agreement is
weaker (mean |error| ≈ 0.05), so check a new fragment against wells you have
integrated yourself.

```python
import ctce
for res in ctce.analyse_plate(["plate_Run01", "plate_Run02"]):
    print(res.summary())          # species, areas, fractions, minor_fraction, flags
```

## Scoring, QC and reports (`scorer.py`)

`scorer.py` is the headless half of the tool and needs neither tkinter nor
matplotlib. It does three jobs:

- **Read QC** — base-call a trace and report length, Qmean/Qmin, the fraction of
  weak calls, N count, peak spacing (mean and CV), per-channel signal maxima and
  a pass/fail flag. Quality values from the bundled caller are *signal heights*
  on a per-preset scale (not PHRED), so only scale-free measures gate pass/fail.
- **Genotype calls** — turn picked-peak rows (the CSV/XLSX/JSON table the
  genotyping picker exports, and the ML training library) into a per-position
  `call` (`het` / `hom-major` / `no-call`), the minor-allele `fraction` and a
  0–100 `confidence`, plus a one-line genotype per well.
- **Fragment-length sizing** — find a ladder's peaks, match them to the kit's
  known lengths, and read every sample peak back through the resulting
  `scan -> bp` curve. See the next section.

There is also a third job, still library-only:

- **Batch auto-genotyping** — **Auto genotyping ▸ Auto-genotype
  sample window…** calls every selected well with no clicking:
  `auto_genotype()` locates the internal-standard quartet with
  `find_is_quartet()`, measures the sample's four duplexes, and returns a call
  per well with a human-readable *reason* when it cannot call one. With the
  standard marked in three wells (A03, H01, D12) it reproduces 93 of the 96
  wells of the T9 plate against the manual ground truth: G09 is the
  documented floor case, and H04/H07 are known regressions since the plate
  model (the pre-model engine called them right). Results land in their own table below the plot — well, call, frac,
  the four duplex areas and their significances, flags and reason — which
  never overwrites the manual pick table, and **Save auto-genotype table…**
  writes them out as CSV/XLSX/JSON.
  Which channel carries the standard and which the sample is a property of the
  assay, not of the dye order, so **Channel roles (standard / sample)…** asks
  for it rather than assuming the T9 defaults. Swapping the pair does not fail
  loudly — it scores the sample's own peaks as the standard and returns
  confident nonsense — so the same channel cannot be set as both.

  Marked standards and calls are kept in the **mark library** (`library/`), so a
  restart never costs a marked plate, a standard fitted on one plate can seed a
  later one (**Auto genotyping ▸ Use a saved standard…**), and **Export library
  for ML…** writes every run's labelled rows out as training data. See
  [MARK_LIBRARY.md](MARK_LIBRARY.md).

```bash
# QC every well in a folder and write a self-contained HTML report
python scorer.py qc example_data/M13/*.rsd --report run.html

# Score an exported peaks table and write the scored rows back out
python scorer.py peaks picks.csv --out scored.csv --report genotypes.html

# Report which dependencies / presets are available
python scorer.py check
```

Console scripts are installed by `pip install -e .`:

```bash
limoncello            # the desktop viewer
limoncello-scorer     # the scorer/QC CLI above
limoncello-basecall   # standalone FASTA/FASTQ basecaller
```

Confidence is a documented **heuristic**, not a calibrated probability: 70% of
the score is how far the minor-allele fraction sits from the call boundary, 30%
is how well the two alleles stand clear of the noise, and an internal-standard
(HOM/HET duplex) match pins it to at least 80.

## Fragment-length sizing (`scorer.py size`)

The other half of genotyping: a PCR/ssDNA sample is co-run with a **size
standard** (ladder) labelled with a different fluorophore, and every sample
peak's length in bp is read off the ladder. This is *sizing*, not allele
calling — it needs no internal standard and no prior knowledge of the sample.

> **Do not size CTCE runs this way.** A CTCE channel's internal standard is
> four same-length conformers of one labelled fragment, not a size ladder, so
> any bp values read off it are meaningless. This cannot be reliably detected
> automatically: on this project's own OY/CTCE data such a channel aligns to
> 15–16 of the 16 GeneScan 500 bands with a leave-one-out error as low as 3.8%,
> *better* than a genuinely realistic GeneScan trace (5.0%). So sizing warns
> when the ladder is weak or few bands match, but a CTCE channel can look like
> a good fit — check the run's method and dyes. Use peak picking and MF for
> CTCE.

```bash
# what ladders are built in, and which MegaBACE kits are compatible
python scorer.py ladders

# size a lane against the bundled 400-1000 bp ROX ladder
python scorer.py size well.rsd --ladder geneflo1000_rox \
    --ladder-channel 4 --sample-channel 2 --out sizes.csv --report sizes.html

# a whole run: the ladder kit and channel are found automatically
python scorer.py size run_folder/*.rsd --out sizes.csv --report sizes.html

# a MegaBACE genotyping plate: ET400-R in ET-ROX on Ch1, FAM sample on Ch2
python scorer.py size plate/*.rsd --ladder et400r \
    --ladder-channel 1 --sample-channel 2 --out sizes.csv

# a ladder that is not bundled: inline lengths or a JSON file
python scorer.py size well.rsd --lengths 50,100,150,200,250,300
python scorer.py size well.rsd --ladder my_ladder.json
```

How it works, and what to watch:

- **Channels are physical.** Ladder and sample live on different channels of
  the same well. They are given as 1–4 in the plate's dye order (`--base-order`,
  default `ACTG`), the *same* convention as the basecaller and the picker.
- **The ladder and its channel are found per run** (`--ladder auto`, the
  default, and no `--ladder-channel`). Wells spread over the run are fitted in
  every channel against the bundled ROX kits (ET400-R, ET550-R, ET900-R,
  ILS 500/600, GeneScan 400HD/500, MapMarker 1000, Beckman SS400) until one
  kit/channel is clean in six of them (at most 32 wells; fewer than three clean
  wells names no kit). ROX bleeds into the
  neighbouring dye and fits there just as well at a fraction of the height,
  so the channel with the tallest bands wins. A clean fit needs most bands
  found, each predicted from the others to a small part of a band spacing,
  **and** few tall peaks left that the kit does not explain — the test that
  tells a ladder from a busy sample channel (a run of 50+ sample peaks can be
  fitted to any kit) and one kit from another that shares most of its
  lengths (ET550-R and ILS600 share 16 of 22). Measured on the lab's 167
  fragment runs (2026-11-25): a true ladder leaves 0–0.09 unexplained tall
  peaks per band, a wrong kit 0.19–0.41, sample channels 1.9–20. The cell ID
  plates are ET550-R, in Ch1 on most runs and Ch3 on some. When two kits fit a
  weak plate about equally the run's note says the kit is uncertain. Naming a
  kit and/or channel restricts the search.
- **Weak wells get help from the plate.** A well whose own fit fails (ladder
  under pull-up, a band missing) is fitted again from where the run's clean
  wells have their bands, allowing for its capillary running faster or slower.
  The rescued bands must stand taller than the other peaks between them
  (1.5×), so a channel of small equal peaks cannot be made into a ladder. A
  rescued well says so in its warnings.
  The check of detection and rescue against all 167 runs, with per-run
  findings, is in [`docs/FRAGMENT_SIZING_QC.md`](docs/FRAGMENT_SIZING_QC.md).
- **A band off its listed length is corrected per run.** Ladder lots differ:
  in the lab's ET550-R runs from September 2012 on, the 125 bp band migrates
  ~2 bp long in nearly every well. When a run has at least 8 clean wells, each
  band is measured against its neighbours (local cubic) across the run. A band
  that sits off by 0.75 to 4 bp, the same way in at least 75 % of the wells,
  gets the length the run measures; at most 3 bands per run, worst first, and
  only if it lowers the run's typical sizing error. The run's note lists the
  corrections; the band table has `band_bp` (listed) and `band_used_bp`.
  `band_correct=False` turns it off.
- **Matching is order-based.** Only peaks of ladder-band height are offered
  (dye bleed from the sample channels sits far lower), and they are aligned to
  the kit's lengths monotonically (a Needleman–Wunsch alignment), so a missing
  or extra band is skipped instead of shifting every length. Several starting
  pairings are tried and the alignment whose anchors lie best on a smooth
  migration curve wins; it is then grown band by band from the anchors it got
  clearly right, which steps over primer/dye peaks ahead of the ladder. The fit
  is a shape-preserving PCHIP curve through the matched anchors — the classic
  `scan -> bp` size-call curve.
- **Accuracy is reported honestly.** A residual on the fit is always zero, so
  the tool reports the **leave-one-out RMS error** in bp (each anchor predicted
  from the others). Below four anchors it reports `n/a` rather than a falsely
  precise zero. `size` warns when fewer than 60 % of the ladder's bands were
  found, and when the bands fit poorly (RMS above 15 % of the spacing to the
  nearest band; a correct ladder sits near 1 %) — usually a wrong ladder
  channel, a wrong length table, or a well with no size standard (e.g. a
  spectral-matrix well). Every sized row also carries its well's
  `ladder_rms_bp` and `ladder_ok`, so a table opened on its own still shows
  which lengths stand on a usable ladder. The GUI (Fragment analysis →
  Size selected wells…) remembers the last ladder and channels used.
- **Out-of-range peaks are flagged.** A peak whose scan falls outside the
  ladder's span is extrapolated linearly and written with `in_range = False`;
  those bp values are the least trustworthy.

### Ladders: the ladder library

Every ladder lives in one library: **bundled** commercial standards, and
**your own** (custom or imported) ladders, kept as one JSON file each in
`library/ladders/` next to the marks and preferences. In the GUI, **Fragment
analysis → Size-standard ladders…** (or *Manage ladders…* in the sizing
dialog) lists them all and can make a **New** ladder, **Import** one,
**Edit**, **Duplicate** (bundled ladders are read-only; duplicate to change
one), **Delete** and **Export** (send the JSON to a colleague). The sizing
dialog picks the ladder from a list, pre-fills its usual channel, and
remembers the last ladder and channels used.

Bundled standards, each with its band table from a primary source (cited in
the ladder's *Source*); bands the vendor says not to size with are excluded
from the curve:

| Ladder | Dye | Bands | Range | Source |
|---|---|---|---|---|
| MegaBACE ET400-R / ET550-R / ET900-R | ET-ROX | 20 / 22 / 37 | 60–400 / 550 / 900 | GE Genetic Profiler 2.2 (`test_plate.tdb`, `gpcom.dll`) |
| GeneScan 120 LIZ | LIZ | 9 | 15–120 | Thermo Fisher Pub. 4474504 |
| GeneScan 350 ROX (250 excluded) | ROX | 12 | 35–350 | " |
| GeneScan 400HD ROX | ROX | 21 | 50–400 | " |
| GeneScan 500 ROX / 500 LIZ (250 excluded) | ROX / LIZ | 16 | 35–500 | " |
| GeneScan 600 LIZ (v1 and v2.0) | LIZ | 36 | 20–600 | " |
| GeneScan 1200 LIZ | LIZ | 68 | 20–1200 | " |
| Promega ILS 600 | CXR | 22 | 60–600 | Promega TM255 |
| Promega WEN / CC5 ILS 500 | WEN | 21 | 60–500 | Promega TMD039 |
| Beckman CEQ Size Standard 400 / 600 | WellRED D1 | 22 / 33 | 60–420 / 640 | Beckman PN A32350 |
| Geneflo 1000 ROX | ROX | 25 | 400–1000 | ChimerX insert |
| MapMarker 400 / 1000 (TAMRA; also ROX, FAM, Cy5, DY-632) | TAMRA | 20 / 23 | 70–400 / 50–1000 | Eurogentec data sheets |
| MapMarker 500 (BV500) / 1500 | ROX | 16 / 28 | 35–500 / 50–1500 | Eurogentec / BioVentures |

A ladder file looks like this; only `lengths` matters to sizing, `exclude`
lists bands that are in the lane but must not anchor the curve, `color` is for
display, and `default_channel` pre-fills the ladder channel:

```json
{"name": "My FAM ladder", "vendor": "home-made", "dye": "FAM", "color": "#1f77b4",
 "lengths": [50, 100, 150, 200, 250], "exclude": [], "default_channel": 2,
 "source": "lab protocol 12", "notes": ""}
```

**Import** reads a ladder JSON (one ladder, a list, or `{"ladders": [...]}`)
or any text/CSV/TSV list of lengths — a row of numbers, one per line, or the
`size`/`bp`/`length` column of a table (a European `;`-separated file with
decimal commas works too) — and a Genetic Profiler `.tdb` table. Every import
opens in the editor first so it can be named and given its dye and colour.
From the command line:

```bash
python scorer.py ladders                                   # list all
python scorer.py ladders --import insert.csv --name "Lab ladder" \
    --dye FAM --default-channel 2
python scorer.py ladders --delete lab_ladder
```

### MegaBACE-compatible standards

The MegaBACE genotyping chemistry uses an **ET ROX** size standard (FAM donor /
ROX acceptor) that is excited by the blue 488 nm laser; the recommended
standards are **ET400-R, ET550-R and ET900-R**, and single-ROX standards also
work but need more label. The validated genotyping dye sets are

| Set | Channels (sample dyes) | Size standard |
|-----|------------------------|---------------|
| Genotyping filter set 1 | FAM · TET · HEX | ET ROX |
| Genotyping filter set 2 | FAM · HEX · NED/TAMRA (or JOE · TAMRA) | ET ROX |

(`MegaBACE Instrument Administrator's Guide`, Appendix A.) In a MegaBACE
genotyping `.rsd` the footer names the dye on each physical channel; on the
Genetic Profiler 2.2 demo plates that is ET-ROX, FAM, NED, HEX on Ch1–Ch4, so
the ladder is `--ladder-channel 1`.

**ET400-R is bundled** (`--ladder et400r`): its 20 band sizes (60, 90, 100,
120, 150, 160, 170, 190, 200, 220, 250, 270, 290, 300, 310, 330, 350, 360, 380,
400 bp) are taken from GE's own Genetic Profiler 2.2 test-plate definition
(`test_plate.tdb` on the GP 2.2 CD). On the two GP 2.2 demo plates every band
is matched in every well that has a ladder, at about 0.15 bp leave-one-out, and
the Genotyping Test Plate fragments (`CS_Ladder_info.tdb`: 13 fragments in
FAM, TAMRA and HEX) are found in every sample well with a ladder (81–82 per
dye) and come out within 0.15 bp of Genetic Profiler's sizes. `tests/test_sizing_gp22_demo.py` checks this when the demo
data is present (`LIMONCELLO_GP22_DEMO`).

### Five-dye runs (LIZ) and `.fsa` files

ABI fragment-analysis files (`.fsa`, and `.ab1`) are read with **every dye
channel** the file holds — five on a LIZ-capable instrument — and the dye
names it records (e.g. 6-FAM, VIC, NED, PET, LIZ). A `.fsa` uses the analyzed
(spectrally separated) data when present, as fragment-analysis software does.
Channels beyond the first four appear as extra switches in the viewer's
channel bar ("Ch5 LIZ"), and sizing takes any physical channel: the LIZ
ladders default to Ch5.

```bash
python scorer.py size run/*.fsa --ladder genescan500_liz --sample-channel 1
```

### Dye overlap, baseline, peak table (matching Genetic Profiler 2.2)

Every dye also shows up in its neighbours' channels; on the MegaBACE
genotyping plates HEX reads at 56–109 % in the NED channel and the ET-ROX
ladder at 8–14 % in FAM and NED. Uncorrected, every strong peak leaves
"pull-up" copies that get sized as fragments. Sizing therefore first removes
the overlap with a spectral matrix (`spectral.py`): the run's own Genetic
Profiler calibration (`<run>/xtm/xtalk.tdb`), else one measured from the
plate's single-dye matrix-standard wells, which are found automatically
(`--matrix auto|run|plate|none`). Measured from GP20's own matrix wells it
reproduces GP's "Demo" calibration to within 0.012 in every entry. Each
channel is then baseline-corrected (running minimum), and sample peaks are
found 50 units above the baseline (`--min-height`), like GP's peak list.

The table (`--out`, or the GUI's save) has one row per ladder band
(`kind = standard`: band, scan, height, width, the size the curve through the
other bands predicts for it, residual) and one per sized peak
(`kind = sample`: channel, dye, scan, size, height, width, area, S/N,
quality 0–100, flags such as *possible pull-up*, *broad*, *outside ladder*,
*flat top (saturated)*). In the GUI the bands and sizes are labelled on the
traces (**Fragment analysis** menu).

Checked against Genetic Profiler 2.2 itself (run under Wine on its demo
plate GP20, ET400-R, Local 3rd Order sizing): the same ladder band to within
one scan in 99.4 % of 1 800 bands; at GP's own peak positions our sizes
differ by a median −0.01 bp (95 % within 0.10 bp), with no dye-specific
offset (GP 2.2 applies no dye-mobility shift). On the ScoreCard test plate
97.7 % of the expected test fragments are found, 95 % within 0.17 bp.

**ET550-R and ET900-R are bundled too** (`et550r`, `et900r`): their band
tables are the ones built into Genetic Profiler 2.2 (`gpcom.dll`), and MRC
Holland's Coffalyser.Net lists the same sizes. They have not yet been checked
against real ET550/ET900 runs here.

## Development

```bash
pip install -e ".[all]" pytest
python -m pytest -q
```

The suite covers the trace loaders, the vendored basecaller, the scorer logic,
the auto-genotyping engine, drag-area measuring and the genotyping picker
engine; genotyping tests skip automatically where tkinter/matplotlib are
unavailable, and the GUI tests (`tests/test_duplex_ui.py`,
`tests/test_auto_genotype_ui.py`) additionally skip without a `DISPLAY`.
Current state: **174 passed, 7 skipped**. CI (`.github/workflows/ci.yml`) runs
the same command on Python 3.9, 3.11 and 3.12.

## Note on environment

This sandbox may not have `tkinter`; run the app on your **local** Windows/macOS/Linux
desktop. To see what is available without starting the GUI, run
`python sequence_analyzer.py --check` (or `python scorer.py check`).
