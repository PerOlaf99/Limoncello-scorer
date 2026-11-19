"""The in-app user manual (Help ▸ User manual…, F1).

Kept out of the GUI module so it can be read, edited and tested on its own.
``SECTIONS`` is a list of ``(title, body)``.  Body markup, rendered by
``sequence_analyzer.LimoncelloAnalyzerApp.show_help``:

* plain lines wrap; a line starting with ``• `` is a bullet;
* ``**text**`` is bold;
* lines between two ``````` lines are shown as-is in a fixed-width
  font (tables, formulas).

Menu paths are written ``Menu ▸ Item`` with the item's exact label, so
``tests/test_user_manual.py`` can check every one of them against the real
menus.  Dialog paths use ``→`` instead.
"""

TITLE = "Limoncello CE Analyzer — user manual"

SECTIONS = [
    ("Start here", """\
Limoncello shows capillary-electrophoresis traces (.rsd, .scf, ABI .ab1/.fsa — up to five dyes, e.g. LIZ on Ch5 — and text/CSV exports) and does three different jobs with them. They are kept apart on purpose — pick the one that matches your run:

• **Base calling** — read the DNA sequence from a sequencing run. → Quick start: base calling.
• **Genotyping (CTC-CE)** — call genotypes from duplex peaks against an internal standard, by clicking peaks or automatically for a whole plate. → Quick start: genotyping.
• **Fragment-length sizing** — give peaks a length in bp against a real size ladder. Never for CTC-CE runs. → Quick start: fragment sizing.

The window:
```
 ┌ menus: File  View  Base calling  Manual genotyping  Auto genotyping  Fragment analysis  Comments  Help
 ├ channel bar: Ch1..Ch4 (Ch5 on five-dye runs) · µA · Time (min) · ▶ Tour + interval · Align · Plate map · Colors
 ├ view bar: Graphs · wells shown · Zoom · Scans · Y scale · Fit Y · Overlay · Align · Plate · Reset
 ├ (while picking) a coloured banner: what the mouse does now · Exit
 │ Data folders │
 │ Samples list │     graphs (one per well, stacked)
 │ Sort, Jump   │
 │              ├──── drag this divider ────
 │              │     table: genotypes / picked peaks / called sequence
 └ status line
```
Fluorescence is shown in volts (V); the instrument current in µA.

New here? Open a folder (File ▸ Add data folder…), click a well in the Samples list, and roll the mouse wheel over the graph to zoom."""),

    ("Quick start: viewing traces", """\
1. File ▸ Add data folder… and choose a run folder — or a parent folder holding a whole project. Its runs appear under Data folders; double-click a run (or select several and press Enter / ▼ List files) to list its files under Samples / files.
2. Click a well. Ctrl-click or Shift-click to select several, then set **Graphs** in the view bar to how many to stack (1–12).
3. Zoom with the mouse wheel on the graph; drag to zoom to a range; double-click to reset. (All gestures: Moving around the graphs.)
4. Page through the plate with ↓ / ↑ (or ← / →, PgDn / PgUp). The view bar shows which wells are on screen, e.g. “A01–A04 · 4 of 96 files”.
5. To compare two runs of the same plate: tick **Overlay runs**. Each graph then shows one well with every run of it on top — Run01 solid, Run02 dashed, Run03 dotted.
6. **▦ Plate** (or View ▸ Plate map…) shows the run as a 96-well plate coloured by genotype call; click a well, a row letter or a column number to show those wells."""),

    ("Quick start: base calling", """\
1. Load the sequencing run and select the wells.
2. Base calling ▸ Basecall selected (or Basecall all in list). The default caller, mb1000_accuracy, suits MegaBACE 1000 runs.
3. The letters appear on the peaks (View ▸ Base-called (with peaks)) and the sequence in the table under the graphs — use **Show sequence** there if the genotyping table is showing.
4. File ▸ Export sequence (FASTA)… writes the reads for BLAST.
• Dye order: sequencing runs declare it in their own files, and that is used automatically (Comments ▸ Run info… shows it). For a run that declares none, the general setting applies (⚙ Basecall settings… → Dye / channel → Base order (instrument); DYEnamic runs are TGCA, ET-primer runs ACTG); Comments ▸ Dye order of this run… overrides it for one run.
• Called a fragment run by mistake? Base calling ▸ Clear base calls (undo) removes the letters; the data files are never changed."""),

    ("Quick start: genotyping (CTC-CE)", """\
Every CTC-CE well carries an internal standard: four peaks of ONE fragment — two homoduplexes and two heteroduplexes (see Manual genotyping: the internal standard). Genotypes are read from the sample's duplex peaks relative to it.

Automatic, for a whole plate:
1. Check Auto genotyping ▸ Channel roles (standard / sample)… — which channel carries the standard and which the sample (the T9 kit: Ch3 standard, Ch2 sample).
2. Mark the standard by hand in **3–5 wells spread across the plate** (different columns; the standard's spacing drifts across a plate): Manual genotyping ▸ Manual peak picking…, click the four standard peaks in each, then Manual genotyping ▸ Mark peaks as standard….
3. Auto genotyping ▸ Semi automatic standard… learns the plate's standard from your marks and places it in the other wells.
4. Auto genotyping ▸ Auto-genotype sample window… calls every well. Results go to the table under the graphs, with a reason for every no-call.
5. Double-click a call to correct it. Auto genotyping ▸ Save auto-genotype table… saves the table; calls and corrections are also kept in the library/ folder.

By hand: Manual genotyping ▸ Manual peak picking…, click each peak, then Manual genotyping ▸ Save peaks table…. Manual genotyping ▸ Peak area by drag… lets you draw the baseline yourself."""),

    ("Quick start: fragment sizing", """\
Only for runs with a real multi-fragment size ladder in its own channel.
1. Select the wells.
2. Fragment analysis ▸ Size selected wells….
3. Pick the ladder from the list (Manage ladders… adds or imports one), its channel, the sample channels (all but the ladder's, or one) and the dye-overlap correction, and press Size. The dialog remembers your last choice.
4. The table (CSV/JSON) lists every ladder band (scan, height, predicted bp) and every sized peak (channel, dye, scan, size, height, width, area, S/N, quality 0–100, flags). The sizes are also labelled on the traces.
5. Read the quality before trusting a length: rows from a well whose ladder did not fit carry ladder_ok = False, and a channel without a ladder gives a warning and no lengths.
Do NOT size CTC-CE runs this way: their internal standard is four conformers of one fragment, not a ladder — see Fragment-length sizing."""),

    ("Data folders and the samples list", """\
Data folders (top left)
• Every run found under an added folder is listed here, but its files are not in the Samples list until you open the run: double-click it, or Ctrl/Shift-click several and press Enter or ▼ List files. Opened runs are ticked (✓) and shaded; double-click one again to take its files out.
• Work on a large project a few runs at a time: open Run01 and Run02, do the work, then the broom (Clear list) and open the next pair.
• – Remove run takes the selected runs out of Data folders; File ▸ Clear folders empties both lists.

Samples / files
• Each line is one trace, shown as run/file. Click to plot; Ctrl/Shift-click to plot several; Ctrl+A selects them all.
• The broom (Clear list) empties the list — the runs stay in Data folders. Hide selected (Delete key) hides single wells.
• **Sort** (the box above the list, or View ▸ Sort files):
```
  Well, row order       A01, A02 … A12, B01 …   runs of one well next to each other
  Well, column order    A01, B01 … H01, A02 …   column by column
  File name             alphabetical
  Run, then file        run by run
```
• **Jump/well** filters the list by well, file or run name as you type.
• – Hide file(s) or the Delete key hides the selected files from the list.
• ▶ Tour (top bar, next to Align) starts the auto-tour, which pages through the list every so many milliseconds (the box beside it); the button turns into ■ Stop while it runs, and Space starts and stops it too. View ▸ Tour interval sets common speeds.
• Rest the pointer on any button for a moment to see what it does."""),

    ("Moving around the graphs", """\
On the graphs
```
  wheel                 zoom in/out across, around the mouse pointer
  Ctrl + wheel          zoom the height
  Shift + wheel         pan left / right
  drag                  zoom to the dragged range (flat drag: across only;
                        a box: across and height)
  Ctrl + drag           the same while picking peaks or measuring areas
                        (there the plain left button is the tool)
  middle-button drag    pan
  double-click, right-click, Home      reset the view
```
The view bar
• **Graphs** — how many graphs to stack (1–12). Beside it: which wells are on screen; in red when more are selected than fit.
• **Zoom − / +** and **Scans from – to**: type a scan range and press Enter. The boxes always show the scans on screen.
• **Y scale**: **Shared** — one height scale for all graphs, to compare wells. **Each graph** — every graph scaled to its own peaks, so a weak well is not flattened by a strong one.
• **Fit Y to view** — the height follows the scans on screen, so zooming into a region also scales up to its peaks.
• **Overlay runs** — one graph per well with every run of that plate on top (Run01 solid, Run02 dashed, Run03 dotted, Run04 dash-dot). Runs of one plate are folders whose names differ only by a trailing RunNN. ↓ / ↑ then page whole wells.
• **▦ Plate** — the plate map (see View menu).
• **⟲ Reset view** — back to the whole trace.

The scroll bars still work: the bottom bar pans across, the right bar up and down; the wheel over a bar zooms that direction; double-click a bar to reset it. All graphs always move together.
X numbers appear under the bottom graph only; tick **Time (min)** to read minutes:seconds instead of scans (1.75 scans per second).

The table under the graphs (genotypes, picked peaks or the called sequence) sits on a divider you can drag. It disappears when the wells on screen have nothing for it; View ▸ Always show table under graphs keeps it.

Aligning graphs
Slide one graph left or right so a feature lines up with the same feature in the others — e.g. two runs of a plate, or wells that ran slightly faster. The data never moves: each graph keeps its true scan numbers (picks, tables and marks are unaffected); only which scans it shows is shifted.
• **Shift + drag** on a graph slides that graph, live.
• **Alt + ← / →** moves the graph under the mouse one scan; **Alt + Shift + ← / →** ten. **Alt + Home** resets it.
• **⇄ Align** (view bar) or View ▸ Alignment:
```
  Align shown graphs to the top one (scans in view)
        lines up what is in view by cross-correlation — zoom onto the
        feature you want aligned first (e.g. the standard's bands);
        with Overlay runs it also aligns each run onto its graph's first
  Align on internal standard
        maps each well's marked IS bands onto the top marked well's:
        offset, plus stretch when 3 or more bands agree
  Reset this graph (Alt+Home) / Reset all alignments
```
• A shifted graph says so in its label, e.g. “⇄ +23 scans” (and “×1.004” when stretched). The tick numbers and the Scans boxes are the bottom graph's own scans. Alignments are saved in the library/ folder per run and well and come back when the plate is opened again; a reset removes them."""),

    ("View menu", """\
Plot styles
```
  Raw traces                     the detector signal as recorded
  Processed (ACGT)               colour-separated channels
  Base-called (with peaks)       letters on the called peaks
  Sequencing trace (ESD peaks)   the base caller's picture: a ring on every
                                 called peak, in its base's colour; a het is
                                 two rings at one scan
  Wrap (one well, N rows)        one long trace over N rows (View ▸ Wrap rows)
```
On/off
• View ▸ Show base letters on peaks, View ▸ Show quality numbers (zoomed), View ▸ Show quality profile (0-100) (a green quality line, off by default), View ▸ Show instrument current (µA) (also the µA tick in the channel bar).
• View ▸ Y scale, View ▸ Fit Y to view, View ▸ Overlay runs of the same well, View ▸ Always show table under graphs — as in the view bar.

Plate map
• View ▸ Plate map… (or ▦ Plate in the view bar) draws one run as a plate — 96 wells, or 384 when the run has wells beyond H or 12 — each coloured by its genotype call: hom-1 blue “1”, hom-2 orange “2”, het green “H”, no-call grey “–”, white = not scored yet, light grey = no file. Wells on screen have a red outline.
• Click a well to show it; Ctrl-click adds or removes one; click a row letter or a column number to show that whole row or column (Graphs is raised to fit). With Overlay runs on, a click shows that well in every run of the plate.
• It stays open beside the main window and follows what is on screen and every new call.

Channels and colours
• The channel bar ticks Ch1–Ch4 in the run's physical order, each labelled with the base it carries in the run on screen (runs can differ — see Base calling menu). Hiding a channel hides it everywhere, including from picking and area measurement.
• View ▸ Trace colors (or the box at the right of the channel bar) chooses a colour scheme:
```
  Classic                    MegaBACE look, colour per BASE:
                             A green · C blue · T red · G black
  Seq DYEnamic (T·G·C·A)     the same colours for a DYEnamic run
  Seq ET primer (A·C·T·G)    the same colours for an ET-primer run
  Genotyping (R·B·Blk·G)     colour per CHANNEL: Ch1 red · Ch2 blue ·
                             Ch3 black · Ch4 green
  Genotyping (G·B·R·Blk)     green · blue · red · black on Ch1–4
  Chromas, High-contrast, Monochrome
```
Paging
• View ▸ Page forward / View ▸ Page backward step through the list (↓ / ↑); View ▸ Start/stop auto-tour plays it. View ▸ Reset view (X & Y) and View ▸ Redraw do what they say."""),

    ("Base calling menu", """\
• Base calling ▸ Basecall selected — calls the wells on screen.
• Base calling ▸ Basecall all in list — calls every listed file (asks first above 96).
• Base calling ▸ ⚙ Basecall settings… — the caller and its tuning; **Apply + redraw** re-runs the shown wells, **Basecall selected with these** calls with the new settings.

Dye order — which base each channel carries — is decided per run:
1. an order chosen for that run (Comments ▸ Dye order of this run…);
2. otherwise the order the run's files declare (sequencing .rsd files do; genotyping runs leave it out);
3. otherwise the general setting: ⚙ Basecall settings… → Dye / channel → Base order (instrument). MegaBACE default ACTG; DYEnamic runs TGCA.
Comments ▸ Run info… says which of the three applies to the well on screen.
• Base calling ▸ Sequencing trace (ESD peaks) / Base calling ▸ Plain channels (processed) — the two plot styles, so a call and its picture sit together.
• Base calling ▸ Clear base calls (undo) — removes letters, peak marks and the quality line (e.g. after calling fragment data by mistake). Data files are never changed.

Callers
```
  mb1000_accuracy   MegaBACE 1000 · longest error-free run, highest identity (default)
  mb1000_length     MegaBACE 1000 · longest read at ≥95 % identity
  mb4000_accuracy   currently the same as mb1000_accuracy (no MegaBACE 4000
  mb4000_length     colour matrix is bundled yet) — named for 4000 data
  pos_bonus07       old name of mb1000_accuracy (kept for saved settings)
  pos_profile       old name of mb1000_length
  hz_soften         mb1000_length with milder sharpening
  raw_peaks         minimal processing, envelope peaks (for checking only)
```"""),

    ("Manual genotyping menu", """\
Two click modes share the main window — the same graphs, zoom, channels and paging; only the left mouse button changes meaning and the table under the graphs becomes the picks table. A coloured banner above the graphs says which mode is on (amber: peak picking, blue: area by drag) and what the mouse does. Leave with **Esc**, the banner's **Exit** button, Manual genotyping ▸ Exit peak picking, or by choosing the mode again.

• Manual genotyping ▸ Manual peak picking… — click a peak (or beside it): it is found, shaded and logged with scan, channel, height and area. Hovering shows a crosshair and a readout of the nearest peak.
• Manual genotyping ▸ Peak area by drag… — hold the left button and drag between two points: the straight line between them is the baseline and the area above it is measured, exactly as drawn. For crowded or odd-shaped peaks. The channel measured is the one with the largest area in the span among the channels switched on — untick the standard's channel to measure only the sample.
• Manual genotyping ▸ Mark peaks as standard… — marks the internal standard in every well you picked peaks in: the earlier two picked peaks become HOM1/HOM2, the later two HET1/HET2. An optional fragment length (bp) can be entered. Marks are saved to library/ and seed auto-genotyping.
• Manual genotyping ▸ Mark IS in all shown wells (also the button in the picking banner) — places the standard in every well on screen from the shape learned by Semi automatic standard…; wells with picks but no mark are marked from their picks. Automatic marks are labelled as such in the table.
• Manual genotyping ▸ Marked IS wells… — lists the wells that carry a mark, i.e. what the plate will be fitted to.
• Manual genotyping ▸ Undo last pick / measurement (Ctrl+Z) and Manual genotyping ▸ Clear picks / measurements. To remove one particular peak: hold R and click it (if it was a standard band it is dropped from the standard too). Right-click does NOT remove — it resets the zoom.
• Manual genotyping ▸ Tag duplex species for MF… — tags the last-picked position's peaks HOM1/HOM2/HET1/HET2 so the MF column can be filled in (below).
• Manual genotyping ▸ Save peaks table… — CSV (Excel-ready), .xlsx or JSON; rows grouped by sample and ordered by scan.
• Manual genotyping ▸ Next batch → / Manual genotyping ▸ ← Previous batch — page by the Graphs count (same as ↓ / ↑). The table keeps every well you visited, so a whole run is saved in one go.
• Manual genotyping ▸ Picking options… ▸ Add +A (A-addition) peak (off by default) tags the Taq A-addition trailing a main peak; leave it off when two alleles sit close together. Manual genotyping ▸ Picking options… ▸ Mark start/end from the 2nd derivative (off by default) marks where each peak lifts off and returns. Manual genotyping ▸ Picking options… ▸ Pick on any channel (auto), or Lock picks to Ch1–Ch4 to pick only on one channel.

The internal standard
The four standard peaks are ONE fragment (one length in bp). Cycling-temperature CE separates them by sequence: two homoduplexes (peaks 1–2, differing by the SNP base) and two heteroduplexes (peaks 3–4, a mismatched Watson/Crick pair). Allele fractions come from the relative areas of these duplex peaks — no size ladder is involved.

Areas and fractions
• A picked peak's area runs between its own two valleys (at most 0.95 × the run's peak spacing), so two close alleles do not share area.
• Two main peaks within max(8 scans, 1.6 × peak spacing) count as a heterozygote pair; the fraction is small / (small + large).
• MF (mass-action fraction) for a tagged position:
```
  MF = (A_MUT + ½·A_HET) / (A_WT + A_MUT + A_HET)
  ai = A_HOMO1 / (A_HOMO1 + A_HOMO2)
```
  A_WT and A_MUT are the homoduplex areas, A_HET both heteroduplexes together. A clean heterozygote reads MF 0.5; below about 5 % the mutant appears only in A_HET. Peak count decides the split: 4 = 2 homo + 2 hetero, 3 = 1 + 2, 2 = 2 + 0."""),

    ("Auto genotyping menu", """\
• Auto genotyping ▸ Semi automatic standard… — learns this plate's internal standard from your marked wells (at least 2; 3–5 spread across the columns is better) and places it in the others. Each placement is reported with how it was found, so position-only finds can be checked.
• Auto genotyping ▸ Use a saved standard… — gives this plate a standard learned on an earlier plate of the same assay, so it needs no marks of its own. A plate that has marks still fits those first.
• Auto genotyping ▸ Forget learned standard — drops the learned or borrowed standard.
• Auto genotyping ▸ Auto-genotype sample window… — calls every well in the list (only your selection, if you made one). The table under the graphs shows well, call, fraction, the four duplex areas and their significance, flags and, for a no-call, the reason. Double-click a call or its IS peaks to correct them; corrections are saved.
• Auto genotyping ▸ Channel roles (standard / sample)… — which channel holds the standard and which the sample. A wrong pair does not fail loudly — it scores the sample as the standard — so both on one channel is refused.
• Auto genotyping ▸ Minimum minor allele for a het… — the smallest minor-allele fraction counted as a heterozygote (0 = the engine's own rule). A smaller second band is reported but not called het. Remembered with the plate.
• Auto genotyping ▸ Save auto-genotype table… — CSV / .xlsx / JSON.
• Auto genotyping ▸ Export library for ML… — every stored call and correction as training rows (see MARK_LIBRARY.md).
• Auto genotyping ▸ CTCE analysis… — cycling temperature CE (see the section CTCE analysis).
Marks, calls and learned standards are kept in the library/ folder beside the program, so a restart does not lose a marked plate."""),

    ("CTCE analysis", """\
Cycling temperature capillary electrophoresis separates the duplexes of a fragment by how they melt (Ekstrøm, Bjørheim & Thilly 2007). Raising the temperature moves the partly melted duplex peaks later, while ssDNA, primers and the labelled oligo keep their relative position. Too cold: all duplexes still run as one peak; in the window: they resolve (four peaks for two variants — two homoduplexes, two heteroduplexes — up to n² for n variants); too hot: they stop moving apart and merge again.

Auto genotyping ▸ CTCE analysis… finds the peaks of interest of a plate and measures them:
• **Plate** — the runs of one plate (folders that differ only by RunNN). Each run's temperature is read from its .rsd, and the runs are ordered by it, so a hotter run done first still shifts the right way.
• **Sample channel** — the analysis uses ONE channel, the strongest by default. **Internal standard** — choose its channel when the plate carries one; its bands then line the runs up (offset and stretch). Without it the runs are lined up on the early primer peak.
• How the peaks of interest are found:
```
  two or more temperatures   peaks that move LATER as it warms, each by its
                             own amount, are the melting duplexes; peaks that
                             hold still (ssDNA drifts slightly earlier) are
                             static
  one temperature            by shape: a duplex has a sharp front and a tail,
  (optimal run, or reruns)   ssDNA is symmetric; reruns count as replicates
                             and the fractions get a ± SD
  your picked peaks          tick "Use my picked peaks as template": click the
                             peaks of interest in one or a few wells with
                             Manual genotyping ▸ Manual peak picking…, and the
                             same pattern is found in every well (it may sit
                             a little earlier or later in each well)
```
• **Measure Homo1 / Homo2 / heteroduplex blocks** (on by default) — the quantification step for the whole plate:
```
  a species       = its main peak + shoulder + Taq stutter (~15 scans
                    later), measured together as one block: 10 scans
                    before the main peak to 30 after, on a near-flat
                    baseline; touching blocks share the valley after
                    the first one's stutter
  merging hets    = Hetero1 and Hetero2 whose valley stays above 20 % of
                    the smaller one are measured as one area under one
                    line (split at the valley only to show each): each
                    counts half towards the homoduplexes, so the split
                    does not matter
  the pattern     = learnt by the plate itself, as migration-time ratios
                    (stable from capillary to capillary): Homo1 the largest
                    species, Homo2 beyond Homo1's cluster at the ratio the
                    strong wells share, then up to two heteroduplexes
  never a species = peaks seen not to move between temperatures, and
                    peaks at the same ratio to the early peak in most
                    wells without a stutter (ssDNA) -- except a peak far
                    taller than anything else there: that is Homo1 sitting
                    where the other run had a peak
  temperature     = the one where most wells resolve all four species
                    (ties: the warmer); reruns there are replicates
```
  The minor-variant fraction is p = (A_Homo2 + ½ A_het) / A_total. Checked against hand-integrated areas of 96 wells (#19_20, 55 °C): Homo1 and Homo2 within 1–3 %, minor fraction within 0.03 in 85 % of wells. On another fragment (#8) the agreement is weaker (within 0.05 in about 70 %): check a new fragment against wells you have integrated yourself. Picked peaks (template) take precedence over the pattern.
• **Out-of-scale peaks** — when the sample channel clips (flat top at the 16-bit limit), the clipped stretch is rebuilt from the strongest other channel that carries the same dye's crosstalk, scaled by the gain fitted on the unclipped part of the same well. The well is flagged "clipped peak measured in channel …".
• **Peaks moving over one another** — as the temperature changes the species move and other peaks do not, so a peak resolved at one temperature can hide under a species block at another and inflate it; a big one changes the minor fraction. Every well is therefore also measured at the plate's other temperatures (where they resolve the same species), each compared with its plate's median there, and a well whose fraction moves differently from the plate is flagged "temperatures disagree …", naming the species whose share moved most — look there, at both temperatures.
• **Add stutter to its parent peak** — (without the block measurement) a smaller peak a few scans behind a duplex peak (Taq stutter / +A) is reported and, ticked, counted with its parent. Peaks you picked are never merged.
• The table: per well the temperatures, method, peaks of interest per run, the run measured, the species (Homo1/Homo2/Hetero1/Hetero2), the fractions (± SD over replicates), the minor-variant fraction and flags (one temperature, merged again at the top temperature, clipped peak, species not found, temperatures disagree, below detection, no product). **Save table…** writes it as CSV.
• **Double-click a well** (or Show well): its runs are overlaid, lined up as the analysis did, with the peaks marked — red ▼ species / peaks of interest (Homo1 … or P1, P2 … when picked), orange ▽ the rest of the block (stutter), grey ○ static.
"""),

    ("Fragment-length sizing", """\
Fragment analysis ▸ Size selected wells… sizes the sample's peaks in bp against a ladder run in its own channel, as Genetic Profiler does. It is a different analysis from peak picking and CTC-CE genotyping.
• Dye overlap: every dye also shows up in its neighbours' channels (on MegaBACE plates HEX reads at 56–109 % in the NED channel), which leaves pull-up copies of strong peaks. Sizing removes it first with a spectral overlap matrix — the run's own Genetic Profiler calibration (xtm/xtalk.tdb) or one measured from the plate's single-dye matrix-standard wells, found automatically. Choose with Fragment analysis ▸ Dye-overlap correction; Fragment analysis ▸ Dye-overlap matrix of this run… shows the matrix, where it came from and the dye on each channel.
• Every channel is baseline-corrected before peaks are found; a peak must stand 50 units above the baseline. Peaks below half the smallest ladder band (the primer/dye blob) are left out.
• Each sized peak gets a quality score 0–100 from the ladder fit, signal/noise, peak shape, dye pull-up, position inside the ladder and clipping; the flags column says which lowered it.
• Fragment analysis ▸ Show fragment sizes on traces labels the ladder bands (bp) and the sized peaks (green: quality ≥ 80, orange: lower); Fragment analysis ▸ Label peaks with quality ≥ sets the cut-off and Fragment analysis ▸ Dye-corrected traces shows the channels with the overlap removed.
• Ladder: picked from the ladder library — bundled commercial standards (MegaBACE ET400-R/ET550-R/ET900-R, GeneScan 120/350/400HD/500/600/1200, Promega ILS 600 and WEN/CC5 ILS 500, Beckman CEQ 400/600, MapMarker 400/500/1000/1500 (Eurogentec/BioVentures), Geneflo 1000) and your own.
• Fragment analysis ▸ Size-standard ladders… manages the library (Fragment analysis ▸ Import ladder… and Fragment analysis ▸ New custom ladder… go straight there): New, Import (a ladder JSON, a list/CSV of lengths, or a Genetic Profiler .tdb), Edit, Duplicate, Delete, Export. Bundled ladders are read-only; duplicate one to change it. Your ladders are kept in library/ladders/.
• A ladder can list bands to exclude from the fit (the vendor's “do not use for sizing” bands are already excluded), a display colour, and its usual channel.
• Channels are physical (1–4, or 5 on a LIZ instrument). On a MegaBACE genotyping plate the ET-ROX ladder is usually Ch1 (Run info shows the dye per channel).
• The fit is reported as leave-one-out error (rms / max, in bp) — an honest number, since a fitted curve always passes its own points.
• No ladder found, or only one fragment matched: you get a warning and NO lengths. Treat that as “not measured”.
• **Never size CTC-CE runs this way.** Their internal standard is a few same-length conformers, not a ladder, and it can line up with a GeneScan table convincingly enough to give believable but meaningless numbers. Use peak picking and MF instead.
• Results: a CSV/JSON table or an HTML report; Fragment analysis ▸ Save the last sizing table… saves it again.
• Checked against Genetic Profiler 2.2 on its own demo plates: the same ladder bands in 99.4 % of cases, sizes within 0.1 bp for 95 % of peaks, no dye-specific offset."""),

    ("Comments and run info", """\
• Comments ▸ Dye order of this run… — fix the dye order for the run of the selected well (Automatic = what its files declare, else the general setting). Its traces are reloaded at once.
• Comments ▸ Run comments… — the comment typed on the instrument when the run was started (read from the run's files), and your own note for the run, kept in the run folder as limoncello-run-note.txt (the data files are not touched). Both are shown under the graphs and in Run info.
• Comments ▸ Run info… — what the file says about the run: well, scan count, run time, dye order, caller, per-channel signal, current; for .rsd also plate ID, machine, method, laser, filters, dye per channel, and the run conditions (voltage, injection, temperature, PMT). A “Dye order note” appears when the file's dye order differs from the one in use — the first thing to check when traces look channel-swapped."""),

    ("Files and export", """\
• File ▸ Export sequence (FASTA)… — the called bases per well.
• File ▸ Export peak table (CSV)… — called bases with scan position and quality.
• File ▸ Export trace text (V + µA)… — the raw values.
• File ▸ Save graph image… — the graphs as PNG, PDF or SVG.
• File ▸ Save settings JSON… / File ▸ Load settings JSON… — your caller setup.
• library/ (beside the program) holds IS marks, calls, corrections, learned standards and your view preferences. Comment files sit beside the data files."""),

    ("Keyboard and mouse", """\
```
  ↓  ↑  ←  →  PgDn  PgUp   page through wells (whole wells with Overlay runs)
  Ctrl+A                   select every sample in the list (not in text boxes)
  Space                    start / stop the auto-tour (in the Samples list)
  Delete                   hide the selected files from the list
  Home                     reset the view
  Alt+← / Alt+→            slide the graph under the mouse one scan
                           (Alt+Shift: ten); Alt+Home resets it
  Ctrl+Z                   undo the last pick or measurement
  R + click                remove that picked peak
  Enter (Scans boxes)      zoom to the typed scan range
  F1                       this manual          F11   full screen / restore
  Esc                      cancel a long job, leave full screen, or leave peak picking

  wheel / Ctrl+wheel / Shift+wheel     zoom across / zoom height / pan
  drag (Ctrl+drag while picking)       zoom to a range
  Shift + drag                         slide that graph (alignment)
  middle-drag                          pan
  double-click / right-click           reset the view
```
In this manual: Ctrl+F jumps to the search box; Enter finds the next match."""),

    ("When something looks wrong", """\
• **Only one graph although several wells are selected** — raise Graphs in the view bar; the text beside it says how many are selected.
• **Channels look swapped** — the dye order is wrong for that run: Comments ▸ Run info… shows which order is used and where it came from; Comments ▸ Dye order of this run… fixes it for that run.
• **A weak well looks flat** — Y scale: Each graph, or zoom in with Fit Y to view on.
• **The table under the graphs is gone** — the wells on screen have nothing for it; View ▸ Always show table under graphs keeps it.
• **Auto-genotyping finds no standard in some wells** — mark a few more wells, spread over different columns, and run Auto genotyping ▸ Semi automatic standard… again. Check Auto genotyping ▸ Channel roles (standard / sample)… too.
• **Letters on a fragment run** — Base calling ▸ Clear base calls (undo).
• **Clicking does not pick a peak** — picking needs Manual genotyping ▸ Manual peak picking… (the amber banner shows it is on); with Ctrl held the click zooms instead.
• **Where is well H07?** — ▦ Plate shows the plate; click the well."""),
]


def plain_text():
    """The manual as plain text (for tests and for copying out)."""
    out = [TITLE, "=" * len(TITLE), ""]
    for i, (title, body) in enumerate(SECTIONS, 1):
        out += [f"{i}. {title}", body.replace("```\n", "").replace("\n```", ""), ""]
    return "\n".join(out)
