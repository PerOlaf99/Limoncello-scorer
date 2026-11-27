#!/usr/bin/env python3
"""Regular fragment-length sizing for capillary-electrophoresis traces.

This is the *sizing* half of genotyping, deliberately separate from the
SNP/allele caller in :mod:`scorer`.  A sample is amplified (PCR) so every
fragment is flanked by the same primer pair; one primer carries a fluorophore.
A size standard ("ladder") -- a set of fragments of known length, labelled with
a *different* fluorophore -- is run in the same lane.  Sizing means:

1. find the ladder's peaks in its channel,
2. match them, in order, to the kit's known fragment lengths,
3. fit a migration curve ``scan -> bp`` through those anchors, and
4. read every sample peak back through the curve to get its length in bp.

Nothing here is MegaBACE-specific except the channel convention: ``acgt`` is
always A,C,G,T while the physical channels follow the plate's dye order
("ACTG" on a MegaBACE), so a channel is always resolved through
:func:`analyzer_core.acgt_index_for_channel`.  The module is headless
(numpy/scipy only) so :mod:`scorer` can use it without Tk or matplotlib.

Ladders
-------
A ladder is just a name, a fluorophore and an ascending list of lengths.  Two
well-known ROX standards are bundled (see :data:`BUILTIN_LADDERS`); any other
kit is accepted either as a JSON file or as a comma-separated list on the
command line, so a lab can use whatever ladder it actually has -- including a
home-made one with a fluorophore of its choice.  MegaBACE-compatible standards
whose fragment tables are *not* bundled are listed in
:data:`KNOWN_MEGABACE_STANDARDS` for reference.

CLI::

    python scorer.py ladders
    python scorer.py size well.rsd --ladder geneflo1000_rox \
        --ladder-channel 4 --sample-channel 2 --out sizes.csv
    python scorer.py size well.rsd --lengths 35,50,75,100,139,150,200,250
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from analyzer_core import acgt_index_for_channel, channel_dye, physical_channel

# numpy dropped ``trapz`` in 2.0 (renamed ``trapezoid``); accept both.
_trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz")

# Physical channel order on a MegaBACE plate.  A channel number is meaningless
# without it: Ch2 is C here, but would be another base on a different dye set.
DEFAULT_BASE_ORDER = "ACTG"

# Sizing tolerances.  A sample peak whose length lands outside the ladder's
# span is reported but flagged: the curve extrapolates there and its accuracy
# degrades quickly past the outermost anchor.
DEFAULT_LADDER_PROMINENCE_FRAC = 0.015
DEFAULT_LADDER_HEIGHT_FRAC = 0.02
DEFAULT_SAMPLE_PROMINENCE_FRAC = 0.02
DEFAULT_SAMPLE_HEIGHT_FRAC = 0.03
# A ladder needs at least this many matched anchors before the sizing curve
# (and its leave-one-out accuracy estimate) mean anything.
MIN_LADDER_ANCHORS = 3

# A size standard's bands have to stand clear of the baseline by a wide margin
# before any of them may be believed.  Measured on this machine's M13 traces the
# injection-region noise is 6-21 units against bands of 800-2100, so a real
# ladder clears 100x comfortably; a channel holding nothing but noise does not.
MIN_LADDER_SNR = 20.0

# A correctly matched ladder predicts each band from the others to a small
# fraction of the spacing to its nearest band (ET400-R on the GP 2.2 demo
# plates: ~0.1-0.2 bp of 10 bp, about 1 %); a band matched to the wrong peak,
# or a well with no size standard, is off by a good part of a spacing.  A
# sparse ladder on a curved migration fits less tightly (synthetic GeneScan 500
# with two merged bands: 13 %).  RMS over the bands, as a fraction of that
# spacing:
MAX_LADDER_RMS_SPACING = 0.15


@dataclass(frozen=True)
class Ladder:
    """A size standard: a name, the dye it is labelled with, and its lengths."""

    name: str
    lengths: Tuple[float, ...]
    dye: str = "ROX"
    source: str = ""
    notes: str = ""
    vendor: str = ""
    color: str = ""                     # display colour, e.g. "#d62728"
    # bands that are in the lane but must not anchor the curve (known
    # anomalous migrators, e.g. GeneScan 500's 250 and 340 bp)
    exclude: Tuple[float, ...] = ()
    default_channel: Optional[int] = None   # usual physical ladder channel

    def __post_init__(self):
        if not self.lengths:
            raise ValueError(f"Ladder {self.name!r} has no fragment lengths")
        if any(b <= 0 for b in self.lengths):
            raise ValueError(f"Ladder {self.name!r} lengths must be positive")
        if list(self.lengths) != sorted(self.lengths):
            object.__setattr__(self, "lengths", tuple(sorted(self.lengths)))
        object.__setattr__(self, "exclude", tuple(float(b) for b in self.exclude))

    @property
    def span(self) -> Tuple[float, float]:
        return self.lengths[0], self.lengths[-1]


@dataclass(frozen=True)
class Peak:
    scan: int
    height: float
    prominence: float


@dataclass
class SizingResult:
    ladder: Ladder
    ladder_channel: int
    sample_channel: int
    base_order: str
    anchors: List[Tuple[int, float]]
    quality: Dict[str, object]
    rows: List[dict]
    well: str = ""
    file: str = ""
    warnings: List[str] = field(default_factory=list)
    standards: List[dict] = field(default_factory=list)   # the ladder bands found
    matrix: str = ""                                       # dye-overlap correction used
    detection: str = ""                                    # how the ladder/channel was chosen
    sample_channels: List[int] = field(default_factory=list)
    band_corrections: Dict[float, float] = field(default_factory=dict)  # listed bp -> bp used
    path: str = ""                                         # the file sized (set by size_files)

    def summary(self) -> dict:
        return {
            "well": self.well,
            "file": self.file,
            "ladder": self.ladder.name,
            "ladder_dye": self.ladder.dye,
            "ladder_channel": self.ladder_channel,
            "sample_channel": self.sample_channel,
            "base_order": self.base_order,
            "n_anchors": len(self.anchors),
            "n_sample_peaks": len(self.rows),
            "in_range_peaks": sum(1 for r in self.rows if r["in_range"]),
            **{f"ladder_{k}": v for k, v in self.quality.items()},
        }


# --------------------------------------------------------------------------- #
# built-in ladders
# --------------------------------------------------------------------------- #
# Only standards whose fragment table can be sourced are bundled.  Everything
# else is supplied by the user (``--lengths`` or a JSON file); see
# KNOWN_MEGABACE_STANDARDS for the MegaBACE-compatible kits and why their
# tables are not hard-coded here.
BUILTIN_LADDERS: Dict[str, dict] = {
    "genescan500_rox": {
        "name": "GeneScan 500 ROX",
        "dye": "ROX",
        "lengths": [35, 50, 75, 100, 139, 150, 160, 200, 250, 300,
                    340, 350, 400, 450, 490, 500],
        "source": "Applied Biosystems GeneScan 500 ROX Size Standard; "
                  "Thermo Fisher DNA Fragment Analysis by CE User Guide "
                  "(Pub. 4474504 Rev B), ch. 2 'Size standards'",
        "vendor": "Applied Biosystems",
        "color": "#d62728",
        "exclude": [250],
        "notes": "16 fragments; 100-mers (100,200,300,400,500) are spiked. "
                 "The vendor says not to use 250 for sizing.",
    },
    "geneflo1000_rox": {
        "name": "Geneflo 1000 ROX",
        "dye": "ROX",
        "lengths": list(range(400, 1001, 25)),
        "source": "ChimerX Geneflo 1000 Size Standard (cat. 3127), insert",
        "vendor": "ChimerX",
        "color": "#d62728",
        "notes": "25 equally spaced bands, 400-1000 bp; 100-mers spiked.",
    },
    "et550r": {
        "name": "MegaBACE ET550-R",
        "dye": "ET-ROX",
        "lengths": [60, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 310,
                    325, 350, 375, 400, 425, 450, 475, 500, 525, 550],
        "source": "GE MegaBACE Genetic Profiler 2.2 (gpcom.dll, built-in "
                  "'ET550-R Size Standard'); same list in MRC Holland's "
                  "Coffalyser.Net size-marker table",
        "vendor": "GE Healthcare / Amersham",
        "color": "#d62728",
        "default_channel": 1,
        "notes": "22 fragments, 60-550 bp. Coffalyser.Net leaves 60, 75 and "
                 "550 out of its sizing curve; GE flags none.",
    },
    "et900r": {
        "name": "MegaBACE ET900-R",
        "dye": "ET-ROX",
        "lengths": [60, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 310,
                    325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600,
                    610, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875,
                    900],
        "source": "GE MegaBACE Genetic Profiler 2.2 (gpcom.dll, built-in "
                  "'ET900-R Size Standard'); same list in MRC Holland's "
                  "Coffalyser.Net size-marker table",
        "vendor": "GE Healthcare / Amersham",
        "color": "#d62728",
        "default_channel": 1,
        "notes": "37 fragments, 60-900 bp, cat. 25-6900-01. Coffalyser.Net "
                 "leaves 60, 75 and 900 out of its sizing curve; GE flags none.",
    },
    "genescan120_liz": {
        "name": "GeneScan 120 LIZ",
        "dye": "LIZ",
        "lengths": [15, 20, 25, 35, 50, 62, 80, 110, 120],
        "source": "Thermo Fisher DNA Fragment Analysis by CE User Guide (Pub. 4474504 Rev B), ch. 2 \'Size standards\'",
        "vendor": "Applied Biosystems",
        "color": "#ff7f0e",
        "default_channel": 5,
        "notes": "9 fragments, 15-120 nt (denatured); for SNaPshot/SNP work.",
    },
    "genescan350_rox": {
        "name": "GeneScan 350 ROX",
        "dye": "ROX",
        "lengths": [35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350],
        "exclude": [250],
        "source": "Thermo Fisher DNA Fragment Analysis by CE User Guide (Pub. 4474504 Rev B), ch. 2 \'Size standards\'",
        "vendor": "Applied Biosystems",
        "color": "#d62728",
        "notes": "12 fragments, 35-350 nt. The vendor says not to use 250 for sizing.",
    },
    "genescan400hd_rox": {
        "name": "GeneScan 400HD ROX",
        "dye": "ROX",
        "lengths": [50, 60, 90, 100, 120, 150, 160, 180, 190, 200, 220, 240,
                    260, 280, 290, 300, 320, 340, 360, 380, 400],
        "source": "Thermo Fisher DNA Fragment Analysis by CE User Guide (Pub. 4474504 Rev B), ch. 2 \'Size standards\'",
        "vendor": "Applied Biosystems",
        "color": "#d62728",
        "notes": "21 fragments, 50-400 nt; high density, for microsatellites.",
    },
    "genescan500_liz": {
        "name": "GeneScan 500 LIZ",
        "dye": "LIZ",
        "lengths": [35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350,
                    400, 450, 490, 500],
        "exclude": [250],
        "source": "Thermo Fisher DNA Fragment Analysis by CE User Guide (Pub. 4474504 Rev B), ch. 2 \'Size standards\'",
        "vendor": "Applied Biosystems",
        "color": "#ff7f0e",
        "default_channel": 5,
        "notes": "16 fragments, 35-500 nt. The vendor says not to use 250 for sizing.",
    },
    "genescan600_liz": {
        "name": "GeneScan 600 LIZ (v1 / v2.0)",
        "dye": "LIZ",
        "lengths": [20, 40, 60, 80, 100, 114, 120, 140, 160, 180, 200, 214,
                    220, 240, 250, 260, 280, 300, 314, 320, 340, 360, 380, 400,
                    414, 420, 440, 460, 480, 500, 514, 520, 540, 560, 580, 600],
        "source": "Thermo Fisher DNA Fragment Analysis by CE User Guide (Pub. 4474504 Rev B), ch. 2 'Size standards'; GeneScan 600 LIZ v2.0 product sheet (Pub. 4425787 Rev C)",
        "vendor": "Applied Biosystems",
        "color": "#ff7f0e",
        "default_channel": 5,
        "notes": "36 fragments, 20-600 nt; v1 and v2.0 have the same peaks. A "
                 "primer peak can mask the smallest band.",
    },
    "genescan1200_liz": {
        "name": "GeneScan 1200 LIZ",
        "dye": "LIZ",
        "lengths": [20, 30, 40, 60, 80, 100, 114, 120, 140, 160, 180, 200, 214,
                    220, 240, 250, 260, 280, 300, 314, 320, 340, 360, 380, 400,
                    414, 420, 440, 460, 480, 500, 514, 520, 540, 560, 580, 600,
                    614, 620, 640, 660, 680, 700, 714, 720, 740, 760, 780, 800,
                    820, 840, 850, 860, 880, 900, 920, 940, 960, 980, 1000, 1020,
                    1040, 1060, 1080, 1100, 1120, 1160, 1200],
        "source": "Thermo Fisher DNA Fragment Analysis by CE User Guide (Pub. 4474504 Rev B), ch. 2 \'Size standards\'",
        "vendor": "Applied Biosystems",
        "color": "#ff7f0e",
        "default_channel": 5,
        "notes": "68 fragments, 20-1200 nt (no 1140).",
    },
    "promega_ils600": {
        "name": "Promega ILS 600",
        "dye": "CXR",
        "lengths": [60, 80, 100, 120, 140, 160, 180, 200, 225, 250, 275, 300,
                    325, 350, 375, 400, 425, 450, 475, 500, 550, 600],
        "source": "Promega MSI Analysis System v1.2 Technical Manual TM255 "
                  "(Rev 2/23), section 1.B",
        "vendor": "Promega",
        "color": "#d62728",
        "notes": "22 fragments, 60-600 bases; multiples of 100 are about twice as bright.",
    },
    "promega_ils500": {
        "name": "Promega WEN / CC5 ILS 500",
        "dye": "WEN",
        "lengths": [60, 65, 80, 100, 120, 140, 160, 180, 200, 225, 250, 275,
                    300, 325, 350, 375, 400, 425, 450, 475, 500],
        "source": "Promega PowerPlex Fusion Technical Manual TMD039 (WEN ILS "
                  "500); Promega customer letter: WEN ILS has the same fragment "
                  "sizes as CC5 ILS 500",
        "vendor": "Promega",
        "color": "#ff7f0e",
        "notes": "21 fragments, 60-500 bases. WEN and CC5 fragments migrate "
                 "slightly differently (matters for bins, not for sizes).",
    },
    "beckman_ss400": {
        "name": "Beckman CEQ Size Standard 400",
        "dye": "WellRED D1",
        "lengths": [60, 70, 80, 90, 100, 120, 140, 160, 180, 190, 200, 220,
                    240, 260, 280, 300, 320, 340, 360, 380, 400, 420],
        "source": "Beckman Coulter GenomeLab GeXP Pre-Installation Guide "
                  "(PN A32350 AB), consumables: DNA Size Standard Kit-400, P/N 608098",
        "vendor": "Beckman Coulter",
        "color": "#9467bd",
        "notes": "22 fragments, 60-420 nt.",
    },
    "beckman_ss600": {
        "name": "Beckman CEQ Size Standard 600",
        "dye": "WellRED D1",
        "lengths": [60, 70, 80, 90, 100, 120, 140, 160, 180, 190, 200, 220,
                    240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460,
                    480, 500, 520, 540, 560, 580, 600, 620, 640],
        "source": "Beckman Coulter GenomeLab GeXP Pre-Installation Guide "
                  "(PN A32350 AB), consumables: DNA Size Standard Kit-600, P/N 608095",
        "vendor": "Beckman Coulter",
        "color": "#9467bd",
        "notes": "33 fragments, 60-640 nt.",
    },
    "mapmarker400": {
        "name": "MapMarker 400 (TAMRA)",
        "dye": "TAMRA",
        "lengths": [70, 80, 90, 100, 120, 140, 160, 180, 190, 200, 220, 240,
                    260, 280, 300, 320, 340, 360, 380, 400],
        "source": "Eurogentec technical data sheet MapMarker 400 "
                  "(GRE-MW-0190-80-0719-V2); BioVentures MM-TMR product page",
        "vendor": "BioVentures / Eurogentec",
        "color": "#7f7f7f",
        "notes": "20 fragments, 70-400 nt. Sold labelled with TAMRA (MW-0190-80TMR, "
                 "MM-TMR), ROX, FAM or DY-632 -- duplicate and change the dye. "
                 "Coffalyser.Net leaves 70-90 out of its sizing curve; the vendor "
                 "flags none.",
    },
    "mapmarker1000": {
        "name": "MapMarker 1000 (TAMRA)",
        "dye": "TAMRA",
        "lengths": [50, 75, 100, 125, 150, 200, 250, 300, 350, 400, 450, 475,
                    500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000],
        "source": "Eurogentec technical data sheet MapMarker 1000 "
                  "(GRE-MW-0195-80-0719-V2); BioVentures MM-TMR 1000 product page",
        "vendor": "BioVentures / Eurogentec",
        "color": "#7f7f7f",
        "notes": "23 fragments, 50-1000 nt. Sold labelled with TAMRA "
                 "(MW-0195-80TMR, MM-TMR 1000), ROX, FAM, Cy5 or DY-632. "
                 "Coffalyser.Net leaves 50 out of its sizing curve; the vendor "
                 "flags none.",
    },
    "mapmarker500": {
        "name": "MapMarker 500 / BV500 (35-500)",
        "dye": "ROX",
        "lengths": [35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350,
                    400, 450, 490, 500],
        "source": "Eurogentec technical data sheet MapMarker 500 "
                  "(GRE-MW-0191-80-0719-V2); BioVentures MM-BV500 product page",
        "vendor": "BioVentures / Eurogentec",
        "color": "#d62728",
        "notes": "16 fragments, 35-500 nt, made to match GeneScan 500's sizes. "
                 "Eurogentec lists 139 (as GeneScan 500); BioVentures' page says "
                 "130, most likely a typo. Sold with ROX or DY-632; no TAMRA version.",
    },
    "mapmarker1500": {
        "name": "MapMarker 1500 (ROX)",
        "dye": "ROX",
        "lengths": [50, 75, 100, 125, 150, 200, 250, 300, 350, 400, 450, 475,
                    500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000,
                    1100, 1200, 1300, 1400, 1500],
        "source": "BioVentures MM-1500-ROX / MM-1500-632 product pages; "
                  "Eurogentec MapMarker product page",
        "vendor": "BioVentures / Eurogentec",
        "color": "#d62728",
        "notes": "28 fragments, 50-1500 nt. Sold with ROX or DY-632; no TAMRA version.",
    },
    "et400r": {
        "name": "MegaBACE ET400-R",
        "dye": "ET-ROX",
        "lengths": [60, 90, 100, 120, 150, 160, 170, 190, 200, 220, 250, 270,
                    290, 300, 310, 330, 350, 360, 380, 400],
        "source": "GE MegaBACE Genetic Profiler 2.2, test_plate.tdb "
                  "(Genotyping Test Plate definition)",
        "vendor": "GE Healthcare / Amersham",
        "color": "#d62728",
        "default_channel": 1,
        "notes": "20 fragments, 60-400 bp.  Checked on the GP 2.2 demo plates: "
                 "every band matched, ~0.15 bp leave-one-out, test-plate "
                 "fragments within 0.15 bp of GP's sizes.",
    },
}

# MegaBACE 1000-compatible size standards.  These are the kits a MegaBACE user
# is most likely to have; their exact fragment tables are printed on the kit
# inserts and are not reproduced from a primary source here (a half-remembered
# ladder is worse than none -- it silently mis-sizes every peak).  Supply the
# lengths from your own insert with ``--lengths`` or a ladder JSON file.
# Resolution figures come from the MegaBACE Long Read Matrix specification.
KNOWN_MEGABACE_STANDARDS = [
    {"key": "et400r", "name": "MegaBACE ET400-R (bundled)", "dye": "ET-ROX",
     "range_bp": "60-400", "resolution_bp": 0.3,
     "catalog": "GE/Cytiva MegaBACE ET400-R Size Standard"},
    {"key": "et550r", "name": "MegaBACE ET550-R (bundled)", "dye": "ET-ROX",
     "range_bp": "~20-550", "resolution_bp": 0.5,
     "catalog": "GE/Cytiva MegaBACE ET550-R Size Standard"},
    {"key": "et900r", "name": "MegaBACE ET900-R (bundled)", "dye": "ET-ROX",
     "range_bp": "~20-900", "resolution_bp": 2.0,
     "catalog": "GE/Cytiva MegaBACE ET900-R, cat. 25-6900-01"},
    {"key": "geneflo1000_rox", "name": "Geneflo 1000 ROX (bundled)",
     "dye": "ROX", "range_bp": "400-1000", "resolution_bp": None,
     "catalog": "ChimerX 3127 (usable as a ROX-channel standard)"},
]


def list_ladders() -> List[str]:
    """Names usable with :func:`load_ladder`: the bundled ladders, then the
    user's own library (:mod:`ladder_library`)."""
    names = sorted(BUILTIN_LADDERS)
    try:
        import ladder_library
        names += [k for k in ladder_library.user_ladders() if k not in names]
    except Exception:                                   # noqa: BLE001
        pass
    return names


def _coerce_lengths(lengths: Sequence) -> Tuple[float, ...]:
    out = []
    for value in lengths:
        if isinstance(value, str):
            value = value.strip()
            if not value:
                continue
        out.append(float(value))
    return tuple(out)


def ladder_from_lengths(lengths: Sequence, name: str = "custom",
                        dye: str = "ROX", source: str = "",
                        notes: str = "") -> Ladder:
    return Ladder(name=name, lengths=_coerce_lengths(lengths), dye=dye,
                  source=source, notes=notes)


def parse_lengths(spec: str) -> Ladder:
    """A comma/space separated list of fragment lengths -> a custom ladder."""
    if not spec or not spec.strip():
        raise ValueError("Empty ladder length list")
    parts = [p for p in spec.replace(";", ",").replace(" ", ",").split(",") if p]
    return ladder_from_lengths(parts, name="custom", dye="custom")


def ladder_from_dict(data: dict, default_name: str = "custom",
                     default_source: str = "") -> Ladder:
    """A :class:`Ladder` from its JSON form (see :mod:`ladder_library`)."""
    if isinstance(data, list):
        data = {"lengths": data}
    ch = data.get("default_channel")
    return Ladder(
        name=str(data.get("name") or default_name),
        lengths=_coerce_lengths(data.get("lengths", [])),
        dye=str(data.get("dye") or "custom"),
        source=str(data.get("source") or default_source),
        notes=str(data.get("notes") or ""),
        vendor=str(data.get("vendor") or ""),
        color=str(data.get("color") or ""),
        exclude=_coerce_lengths(data.get("exclude") or []),
        default_channel=int(ch) if ch not in (None, "") else None,
    )


def load_ladder(spec) -> Ladder:
    """Resolve *spec* to a :class:`Ladder`.

    *spec* may be a bundled or library ladder name (case-insensitive, see
    :mod:`ladder_library`), a comma-separated list of lengths, or a path to a
    JSON file::

        {"name": "My ROX ladder", "dye": "ROX",
         "lengths": [50, 100, 150, 200], "notes": "..."}
    """
    if isinstance(spec, Ladder):
        return spec
    if not isinstance(spec, str) or not spec.strip():
        raise ValueError("No ladder given")

    text = spec.strip()
    path = Path(text)
    if path.is_file():
        data = json.loads(path.read_text(encoding="utf-8"))
        return ladder_from_dict(data, default_name=path.stem,
                                default_source=str(path))

    key = text.lower()
    if key in BUILTIN_LADDERS:
        return ladder_from_dict(BUILTIN_LADDERS[key])
    try:
        import ladder_library
        found = ladder_library.find(text)
    except Exception:                                   # noqa: BLE001
        found = None
    if found is not None:
        return found

    if any(ch.isdigit() for ch in text):
        return parse_lengths(text)

    raise ValueError(
        f"Unknown ladder {spec!r}. Known: {', '.join(list_ladders())} "
        "(or pass lengths / a JSON path)")


# --------------------------------------------------------------------------- #
# peak detection (numpy/scipy only -- mirrors genotyping's picker, not imports)
# --------------------------------------------------------------------------- #
def _baseline_noise(trace) -> float:
    """Robust per-point noise SD of the flat injection region.

    The first few hundred scans are before anything has entered the capillary,
    so this is baseline.  Measured after the same median/Savitzky-Golay
    smoothing ``detect_peaks`` applies, because that is the signal the
    thresholds actually see.
    """
    from scipy.signal import medfilt, savgol_filter

    y = np.asarray(trace, dtype=float).ravel()
    if y.size < 32:
        return 0.0
    head = y[:min(400, y.size)]
    if head.size >= 9:
        head = savgol_filter(medfilt(head, 5), 9, 2, mode="interp")
    # MAD -> SD, so a stray spike in the injection region does not set the scale
    mad = float(np.median(np.abs(head - np.median(head))))
    return mad * 1.4826


def detect_peaks(trace, *, min_prominence_frac: float = 0.02,
                 min_height_frac: float = 0.03, cut: int = 0,
                 smooth: bool = True) -> List[Peak]:
    """Prominent local maxima of one channel, above a robust local baseline.

    The baseline is the median of the first ~400 scans (the injection region),
    which is where the trace is flattest; prominence/height thresholds are
    fractions of the largest baseline-corrected sample so they scale with the
    run's gain.
    """
    from scipy.signal import find_peaks

    y = np.asarray(trace, dtype=float).ravel()
    n = y.size
    if n == 0:
        return []
    if smooth and n >= 9:
        from scipy.signal import medfilt, savgol_filter
        y = savgol_filter(medfilt(y, 5), 9, 2, mode="interp")
    base = float(np.median(y[:min(400, n)]))
    z = y - base
    mx = float(z.max())
    if not np.isfinite(mx) or mx <= 0:
        return []
    idx, props = find_peaks(z, prominence=mx * min_prominence_frac,
                            height=mx * min_height_frac)
    peaks = []
    for k, i in enumerate(idx):
        if i < cut:
            continue
        peaks.append(Peak(int(i), float(z[i]), float(props["prominences"][k])))
    return peaks


def _peak_areas(trace, peaks: Sequence[Peak]) -> List[float]:
    """Trapezoid area of each peak down to the local valley between neighbours.

    The integration window runs from the midpoint toward each neighbour to the
    next, and the baseline is the floor of that window -- so two partially
    overlapping peaks split their shared valley instead of double-counting it.
    """
    y = np.asarray(trace, dtype=float).ravel()
    scans = [p.scan for p in peaks]
    areas = []
    for i, scan in enumerate(scans):
        left = int((scans[i - 1] + scan) / 2) if i > 0 else max(0, scan - 20)
        right = (int((scan + scans[i + 1]) / 2) if i + 1 < len(scans)
                 else min(len(y), scan + 20))
        seg = y[left:max(right, left + 1)]
        if seg.size < 2:
            areas.append(0.0)
            continue
        floor = float(np.min(seg))
        areas.append(max(0.0, float(_trapz(seg - floor, dx=1.0))))
    return areas


# --------------------------------------------------------------------------- #
# ladder <-> length matching and the sizing curve
# --------------------------------------------------------------------------- #
def _log_linear_fit(u: np.ndarray, v: np.ndarray) -> Tuple[float, float]:
    """Least-squares ``v = a + b*u`` with a guaranteed positive slope.

    In capillary electrophoresis scan number grows roughly with log(bp), so
    this is a better initial guess than a raw linear fit; ``b`` is forced
    positive because longer fragments must migrate later, never earlier.
    """
    if u.size < 2 or float(np.ptp(u)) == 0.0:
        return 0.0, 1.0
    b, a = np.polyfit(u, v, 1)
    if not np.isfinite(b) or b <= 0:
        b = max((v[-1] - v[0]) / (u[-1] - u[0]), 1e-9)
        a = float(v[0] - b * u[0])
    return float(a), float(b)


def _align_dp(scans: Sequence[int], lengths: Sequence[float],
              pred: np.ndarray) -> List[Tuple[int, int]]:
    """Monotonic alignment of detected peaks to expected lengths.

    A Needleman-Wunsch alignment: match, skip a detected peak, or skip an
    expected length.  The match cost is the scan distance to that length's
    predicted position, and every gap costs one average peak spacing, so a
    real peak is always cheaper to match than to skip.
    """
    n, k = len(scans), len(lengths)
    if n == 0 or k == 0:
        return []
    span = max(float(scans[-1]) - float(scans[0]), 1.0)
    gap = span / max(n, k)
    # plain Python lists, one row at a time: this runs tens of thousands of
    # times per plate, and element access on numpy arrays is ~7x slower
    pred = [float(p) for p in pred]
    sc = [float(s) for s in scans]
    prev = [0.0] * (k + 1)
    for j in range(1, k + 1):
        prev[j] = prev[j - 1] + gap
    back = [[2] * (k + 1)]          # 0 match, 1 skip scan, 2 skip length
    for i in range(1, n + 1):
        s = sc[i - 1]
        row = [prev[0] + gap] + [0.0] * k
        b = [1] + [0] * k
        for j in range(1, k + 1):
            best, c = prev[j - 1] + abs(pred[j - 1] - s), 0
            t = prev[j] + gap
            if t < best:
                best, c = t, 1
            t = row[j - 1] + gap
            if t < best:
                best, c = t, 2
            row[j] = best
            b[j] = c
        back.append(b)
        prev = row
    pairs = []
    i, j = n, k
    while i or j:
        choice = back[i][j]
        if choice == 0 and i and j:
            pairs.append((i - 1, j - 1))
            i, j = i - 1, j - 1
        elif choice == 1:
            i -= 1
        else:
            j -= 1
    pairs.reverse()
    return pairs


def _align_from(scans: List[int], lengths: List[float], coef: np.ndarray,
                iterations: int, tf=np.log) -> List[Tuple[int, int]]:
    """Iterate DP alignment against ``scan = polyval(coef, tf(bp))`` and a
    refit of that model on the matches: linear to start, quadratic once five
    or more anchors can carry the curvature (a straight line through 60-400
    bp misplaces the end bands by more than a band spacing)."""
    best: List[Tuple[int, int]] = []
    tl = tf(np.asarray(lengths, float))
    for _ in range(max(1, iterations)):
        pairs = _align_dp(scans, lengths, np.polyval(coef, tl))
        if len(pairs) < 2:
            best = pairs or best
            break
        u = tl[[j for _, j in pairs]]
        v = np.asarray([scans[i] for i, _ in pairs], float)
        if len(pairs) >= 5:
            coef = np.polyfit(u, v, 2)
        else:
            a, b = _log_linear_fit(u, v)
            coef = np.array([b, a])
        if pairs == best:
            break
        best = pairs
    return best


def _trusted(scans: List[int], lengths: List[float],
             pairs: List[Tuple[int, int]], tol_bp: float = 1.0) -> List[Tuple[int, int]]:
    """The anchors of *pairs* that the curve through all the others predicts
    to within *tol_bp*."""
    xs = [scans[i] for i, _ in pairs]
    ys = [lengths[j] for _, j in pairs]
    keep = []
    for h, pr in enumerate(pairs):
        try:
            c = SizeCurve(xs[:h] + xs[h + 1:], ys[:h] + ys[h + 1:])
        except ValueError:
            continue
        if abs(float(c.bp([xs[h]])[0]) - ys[h]) <= tol_bp:
            keep.append(pr)
    return keep


def _extend(scans: List[int], lengths: List[float],
            pairs: List[Tuple[int, int]], tol: float = 0.35) -> List[Tuple[int, int]]:
    """Grow an alignment band by band from anchors already trusted.

    Each unmatched band next to a matched one is predicted from its nearest
    four anchors (a local quadratic in bp) and takes the closest unused peak
    lying between its matched neighbours, if that peak is within *tol* of a
    band spacing of the prediction.  Local fits follow the migration curve
    where one global curve extrapolated over a hundred bases cannot."""
    match = {j: i for i, j in pairs}
    k = len(lengths)
    bp = np.asarray(lengths, float)
    changed = True
    while changed:
        changed = False
        for j in range(k):
            if j in match or not ((j - 1) in match or (j + 1) in match):
                continue
            near = sorted(match, key=lambda m: abs(bp[m] - bp[j]))[:4]
            if len(near) < 2:
                continue
            u = bp[near]
            v = np.asarray([scans[match[m]] for m in near], float)
            coef = np.polyfit(u, v, 2 if len(near) >= 4 else 1)
            pred = float(np.polyval(coef, bp[j]))
            slope = float(np.polyval(np.polyder(coef), bp[j]))
            if slope <= 0:
                continue
            nb = [bp[m] for m in (j - 1, j + 1) if 0 <= m < k]
            spacing = slope * min(abs(bp[j] - x) for x in nb)
            lo = max((scans[match[m]] for m in match if m < j), default=-np.inf)
            hi = min((scans[match[m]] for m in match if m > j), default=np.inf)
            used = set(match.values())
            cand = [i for i, sc in enumerate(scans) if lo < sc < hi and i not in used]
            if not cand:
                continue
            i = min(cand, key=lambda i: abs(scans[i] - pred))
            if abs(scans[i] - pred) <= tol * spacing:
                match[j] = i
                changed = True
    return sorted(((i, j) for j, i in match.items()), key=lambda p: p[1])


def _alignment_score(scans: List[int], lengths: List[float],
                     pairs: List[Tuple[int, int]], n_scans: Optional[int] = None) -> float:
    """Lower is better: the scatter (bp) of the anchors about a smooth
    migration curve -- a quadratic in bp or in log(bp), whichever fits --
    plus 1 bp per ladder band left unmatched.

    A global smooth curve rather than the PCHIP leave-one-out: on a sparse
    ladder PCHIP through a handful of points cannot follow the true shape,
    and its end anchors can only be extrapolated, so it would favour
    alignments that drop bands."""
    if len(pairs) < 3:
        return float("inf")
    v = np.asarray([scans[i] for i, _ in pairs], float)
    bp = np.asarray([lengths[j] for _, j in pairs], float)
    scans_per_bp = max((v[-1] - v[0]) / max(bp[-1] - bp[0], 1e-9), 1e-9)
    if scans_per_bp <= 1e-9 or np.any(np.diff(v) <= 0):
        return float("inf")
    best, best_fit = float("inf"), None
    for tf in (np.log, np.asarray):
        u = tf(bp)
        deg = 2 if len(pairs) >= 4 else 1
        coef = np.polyfit(u, v, deg)
        resid = v - np.polyval(coef, u)
        err = float(np.sqrt(np.mean(resid ** 2))) / scans_per_bp
        if err < best:
            best, best_fit = err, (tf, coef)
    matched = {j for _, j in pairs}
    missing = [b for j, b in enumerate(lengths) if j not in matched]
    if n_scans is not None and missing:
        # bands the curve puts past the end of the trace are not missing:
        # the run stopped before them (a long kit on a short run).  Without
        # this, an evenly spaced kit slides by a few bands to "explain" its
        # unreached long end with peaks that are really its short end.
        tf, coef = best_fit
        last = max(bp)
        missing = [b for b in missing
                   if b < last or float(np.polyval(coef, tf(np.asarray([b], float)))[0]) <= n_scans]
    return best + len(missing)


def align_peaks(scans: Sequence[int], lengths: Sequence[float],
                iterations: int = 5, n_scans: Optional[int] = None) -> List[Tuple[int, float]]:
    """Match ladder peak scans to the kit's known lengths.

    Bootstraps a migration model (scan linear in log(bp), or in bp), aligns
    with DP, refits the model on the matches, and repeats -- the iteration is what absorbs the curvature
    a single global fit cannot.  One stray peak at either end (a primer peak,
    dye bleed from a sample channel) would wreck a bootstrap taken from the
    outermost peaks, so several end pairings are tried -- each of the first and
    last few peaks against each of the first and last few lengths -- and the
    alignment whose anchors sit best on a smooth curve wins; a repair pass then
    re-aligns from the anchors that alignment got clearly right.
    *n_scans* (the trace length) lets bands past the end of the run go
    unmatched without penalty.  Returns ``(scan, length)`` anchors.
    """
    scans = [int(s) for s in scans]
    lengths = [float(b) for b in lengths]
    if not scans or not lengths:
        return []
    n, k = len(scans), len(lengths)
    # Migration is close to log(bp) for short fragments and close to linear
    # in bp over a few hundred bases; seed with both and let the fit decide.
    starts = []
    for tf in (np.log, np.asarray):
        tl = tf(np.asarray(lengths, float))
        if n == k:
            starts.append((tf, _log_linear_fit(tl, np.asarray(scans, float))))
        ends = 3
        for i0 in range(min(ends, n)):
            for i1 in range(max(i0 + 1, n - ends), n):
                for j0 in range(min(2, k)):
                    for j1 in range(max(j0 + 1, k - 2), k):
                        starts.append((tf, _log_linear_fit(
                            np.asarray([tl[j0], tl[j1]]),
                            np.asarray([scans[i0], scans[i1]], float))))
    best, best_score, seen = [], float("inf"), set()

    def consider(pairs):
        nonlocal best, best_score
        key = tuple(pairs)
        if key in seen:
            return
        seen.add(key)
        score = _alignment_score(scans, lengths, pairs, n_scans)
        if score < best_score or not best:
            best, best_score = pairs, score

    for tf, (a, b) in starts:
        consider(_align_from(scans, lengths, np.array([b, a]), iterations, tf))
    # Walks: assume four consecutive peaks at one end are the four end bands
    # (allowing for a few stray peaks or a missing band there) and grow the
    # alignment band by band from them.  A local prediction steps over dye
    # bleed sitting between bands, which a global curve fit does not.
    seed = 4
    if n >= seed and k >= seed:
        for s in range(min(3, n - seed + 1)):
            for b in range(min(2, k - seed + 1)):
                consider(_extend(scans, lengths,
                                 [(n - 1 - s - t, k - 1 - b - t) for t in range(seed)][::-1]))
        for s in range(min(max(3, n - k + 3), n - seed + 1)):
            for b in range(min(2, k - seed + 1)):
                consider(_extend(scans, lengths, [(s + t, b + t) for t in range(seed)]))
    # Repair: junk ahead of the ladder (primer and dye peaks) can push every
    # seed to pair the first bands one peak too early.  The anchors the rest
    # of the curve predicts well are right; keep those and grow the alignment
    # from them band by band (see _extend).
    for _ in range(3):
        good = _trusted(scans, lengths, best)
        if len(good) < 5 or len(good) == len(best):
            break
        before = best_score
        consider(_extend(scans, lengths, good))
        if best_score >= before:
            break
    return [(scans[i], lengths[j]) for i, j in best]


# Ladder bands of one kit come out at broadly similar heights; dye bleed from
# the sample channels, primer remnants and noise sit well below them.  Peaks
# under this fraction of a typical band (the median of the band-count tallest
# peaks) are not offered to the aligner.
LADDER_BAND_REL_HEIGHT = 0.3


def ladder_candidates(peaks: Sequence[Peak], n_bands: int,
                      rel: float = LADDER_BAND_REL_HEIGHT) -> List[Peak]:
    """The detected peaks tall enough to be bands of an *n_bands* ladder."""
    if not peaks or n_bands <= 0:
        return list(peaks)
    tall = sorted((p.height for p in peaks), reverse=True)[:n_bands]
    ref = float(np.median(tall))
    keep = [p for p in peaks if p.height >= rel * ref]
    # a ladder fading badly at one end: loosen rather than drop real bands
    if len(keep) < 0.6 * min(n_bands, len(peaks)):
        keep = [p for p in peaks if p.height >= 0.5 * rel * ref]
    return keep


class SizeCurve:
    """Piecewise-cubic ``scan -> bp`` curve through the ladder anchors.

    PCHIP is used because it is shape-preserving: the curve stays monotone
    between anchors and never overshoots into a non-physical length.  Outside
    the anchor span it extrapolates linearly from the outermost segment, and
    callers are expected to flag those points (they are the least trustworthy).
    """

    def __init__(self, scans: Sequence[int], lengths: Sequence[float]):
        from scipy.interpolate import PchipInterpolator

        pairs = sorted(zip((int(s) for s in scans), (float(b) for b in lengths)))
        xs, ys = [], []
        for x, y in pairs:                       # dedupe scans, keep first
            if xs and x <= xs[-1]:
                continue
            xs.append(x)
            ys.append(y)
        if len(xs) < 2:
            raise ValueError("A sizing curve needs at least two ladder anchors")
        self.scan_min, self.scan_max = float(xs[0]), float(xs[-1])
        self._pchip = PchipInterpolator(xs, ys, extrapolate=False)
        self._xs = np.asarray(xs, float)
        self._ys = np.asarray(ys, float)

    def bp(self, scan) -> np.ndarray:
        scan = np.atleast_1d(np.asarray(scan, dtype=float))
        out = np.empty_like(scan)
        inside = (scan >= self.scan_min) & (scan <= self.scan_max)
        if inside.any():
            out[inside] = self._pchip(scan[inside])
        if (~inside).any():
            out[~inside] = self._extrapolate(scan[~inside])
        return out

    def _extrapolate(self, scan: np.ndarray) -> np.ndarray:
        lo = scan < self.scan_min
        hi = ~lo
        out = np.empty_like(scan)
        if lo.any():
            slope = ((self._ys[1] - self._ys[0]) /
                     (self._xs[1] - self._xs[0]))
            out[lo] = self._ys[0] + slope * (scan[lo] - self._xs[0])
        if hi.any():
            slope = ((self._ys[-1] - self._ys[-2]) /
                     (self._xs[-1] - self._xs[-2]))
            out[hi] = self._ys[-1] + slope * (scan[hi] - self._xs[-1])
        return out

    @property
    def in_range(self):
        def _check(scan) -> np.ndarray:
            s = np.atleast_1d(np.asarray(scan, dtype=float))
            return (s >= self.scan_min) & (s <= self.scan_max)
        return _check


def ladder_fit_quality(scans: Sequence[int],
                       lengths: Sequence[float]) -> Dict[str, object]:
    """Leave-one-out accuracy of the sizing curve, in bp.

    A residual on the fit itself is identically zero (PCHIP passes through its
    anchors), so the honest measure is how well each anchor is predicted from
    the others.  Needs at least four anchors to say anything; below that the
    errors are reported as ``None`` rather than as a falsely precise zero.
    """
    scans = [int(s) for s in scans]
    lengths = [float(b) for b in lengths]
    n = len(scans)
    quality: Dict[str, object] = {
        "anchors": n,
        "rms_error_bp": None,
        "max_error_bp": None,
        "rms_error_rel": None,
        "rms_error_spacing": None,
    }
    if n < 4:
        return quality
    errors = []
    for hold in range(n):
        xs = [scans[i] for i in range(n) if i != hold]
        ys = [lengths[i] for i in range(n) if i != hold]
        try:
            curve = SizeCurve(xs, ys)
        except ValueError:
            continue
        errors.append(float(curve.bp([scans[hold]])[0] - lengths[hold]))
    if errors:
        err = np.abs(errors)
        denom = np.abs(lengths)
        quality["rms_error_bp"] = float(np.sqrt(np.mean(np.square(errors))))
        quality["max_error_bp"] = float(err.max())
        quality["rms_error_rel"] = float(np.sqrt(np.mean((err / denom) ** 2)))
        # each error as a fraction of the distance to the nearest other band:
        # a band matched to the wrong peak is off by a sizeable part of a
        # spacing, whatever the ladder's density
        ls = np.asarray(lengths, float)
        gaps = np.abs(ls[:, None] - ls[None, :])
        np.fill_diagonal(gaps, np.inf)
        quality["rms_error_spacing"] = float(np.sqrt(np.mean((err / gaps.min(axis=1)) ** 2)))
    return quality


# --------------------------------------------------------------------------- #
# the whole plate: which ladder, in which channel, and where its bands sit
# --------------------------------------------------------------------------- #
# A share of the kit's bands this large, fitted this tightly, is a ladder.
MIN_BAND_SHARE = 0.8

# A ladder channel holds the ladder and little else: tall peaks inside the
# ladder's span that the kit does not account for, per band of the kit.
# Measured (2026-11-25, cell ID and non-ladder runs): ET550-R in its own
# channel 0-0.09; a wrong kit on that ladder (ILS600/ILS500 on ET550-R)
# 0.19-0.41; sample and noise channels of runs without a ladder 1.9-20.
MAX_UNEXPLAINED = 0.25


def unexplained_tall(peaks: Sequence[Peak], anchors: Sequence[Tuple[int, float]],
                     n_bands: int) -> Optional[float]:
    """Tall peaks (at least half the median matched band) that are not
    anchors, from the first anchor to the end of the trace, per band of the
    kit.  Past the last anchor counts too: a shorter kit nested in the real
    one (ET550-R in ET900-R, ILS500 in ILS600) fits the bands it shares and
    would otherwise leave the rest of the ladder unseen.  Ahead of the first
    anchor is not counted -- primer and dye peaks belong there."""
    if not anchors or not n_bands:
        return None
    h = {p.scan: p.height for p in peaks}
    ms = {int(s) for s, _ in anchors}
    band = float(np.median([h.get(s, 0.0) for s in ms]))
    lo = min(ms)
    n = sum(1 for sc, ht in h.items() if lo <= sc and sc not in ms and ht >= 0.5 * band)
    return n / float(n_bands)


# The ROX-channel standards a MegaBACE lab is likely to have used (a
# narrower search, e.g. for a MegaBACE-only archive).
MEGABACE_KITS = ("et400r", "et550r", "et900r", "promega_ils600", "promega_ils500",
                 "genescan500_rox", "genescan400hd_rox", "mapmarker1000",
                 "beckman_ss400")
DETECT_KITS = MEGABACE_KITS        # kept for older callers


def auto_kits() -> List[str]:
    """Kits tried when the ladder is "auto": every bundled ladder and every
    ladder in the user's own library (:func:`list_ladders`)."""
    return list_ladders()


def _dye_fits_channel(lad: "Ladder", channel: int, n_channels: int) -> bool:
    """A LIZ kit sits in the fifth channel of a five-dye instrument; every
    other kit does not.  Decides between kits whose tables are identical
    (GeneScan 500 ROX / LIZ, MapMarker 500)."""
    liz = (lad.dye or "").upper() == "LIZ"
    return liz == (n_channels >= 5 and channel == 5)


# a run's kit/channel needs a clean fit in at least this many test wells
MIN_DETECT_WELLS = 3


# a run that stopped early may show only part of a kit -- but at least this share
MIN_KIT_REACHED = 0.6


def bands_reached(anchors: Sequence[Tuple[int, float]], lengths: Sequence[float],
                  n_scans: int) -> int:
    """How many of the kit's bands the run reached: those not longer than
    what the curve through *anchors* puts at the trace's last scan.  A run
    stopped before a long kit's last bands cannot show them, and they must
    not count as missing."""
    if len(anchors) < 2 or n_scans <= 0:
        return len(lengths)
    try:
        end_bp = float(SizeCurve([s for s, _ in anchors], [b for _, b in anchors])
                       .bp([n_scans - 1])[0])
    except ValueError:
        return len(lengths)
    return sum(1 for b in lengths if b <= end_bp + 1e-6)


def ladder_fit_ok(quality: Dict[str, object], n_bands: int) -> bool:
    """Most of the kit's bands found, each predicted from the others to a
    small part of a band spacing, and (when measured) few tall peaks in the
    channel left unexplained by the kit (:data:`MAX_UNEXPLAINED`).  When the
    run stopped before the kit's last bands (``bands_reached``), only the
    bands it reached count -- provided that is most of the kit."""
    reached = quality.get("bands_reached")
    if reached is not None and reached < n_bands:
        if reached < MIN_KIT_REACHED * n_bands:
            return False
        n_bands = int(reached)
    n = int(quality.get("anchors") or 0)
    rel = quality.get("rms_error_spacing")
    unexpl = quality.get("unexplained_tall")
    return (n >= max(MIN_LADDER_ANCHORS + 1, MIN_BAND_SHARE * n_bands)
            and rel is not None and rel <= MAX_LADDER_RMS_SPACING
            and (unexpl is None or unexpl <= MAX_UNEXPLAINED))


def fit_ladder_channel(trace, ladder, peaks: Optional[List[Peak]] = None) -> Dict[str, object]:
    """Ladder fit of one channel alone (no sample sizing): the anchors, their
    leave-one-out quality, and the median height of the matched bands.
    *peaks*: the channel's :func:`detect_peaks` at ladder settings, if known."""
    if isinstance(ladder, str):
        ladder = load_ladder(ladder)
    if peaks is None:
        peaks = detect_peaks(trace, min_prominence_frac=DEFAULT_LADDER_PROMINENCE_FRAC,
                             min_height_frac=DEFAULT_LADDER_HEIGHT_FRAC)
    noise = _baseline_noise(trace)
    top = max((p.height for p in peaks), default=0.0)
    out = {"anchors": [], "quality": ladder_fit_quality([], []), "band_height": 0.0,
           "ok": False}
    if not peaks or noise <= 0 or top < MIN_LADDER_SNR * noise:
        return out
    cand = ladder_candidates(peaks, len(ladder.lengths))
    n = len(ladder.lengths)
    if len(cand) < MIN_BAND_SHARE * MIN_KIT_REACHED * (n - len(ladder.exclude)) - 0.5:
        # too few band-height peaks for this kit: no alignment can find the
        # share of its bands a clean fit needs
        out["sparse"] = True
        return out
    if len(cand) > (1.0 + MAX_UNEXPLAINED) * n + 12:
        # far more band-height peaks than the kit has bands: whatever fits,
        # too many tall peaks would be left unexplained -- skip the alignment
        out["crowded"] = True
        return out
    matched = align_peaks([p.scan for p in cand], ladder.lengths, n_scans=len(trace))
    anchors = [(s, b) for s, b in matched
               if not any(abs(b - x) < 1e-6 for x in ladder.exclude)]
    q = ladder_fit_quality([s for s, _ in anchors], [b for _, b in anchors])
    q["unexplained_tall"] = unexplained_tall(peaks, anchors, len(ladder.lengths))
    q["bands_reached"] = bands_reached(anchors, [b for b in ladder.lengths
                                                 if b not in ladder.exclude], len(trace))
    h = {p.scan: p.height for p in cand}
    out.update(anchors=anchors, quality=q,
               band_height=float(np.median([h[s] for s, _ in anchors])) if anchors else 0.0,
               ok=ladder_fit_ok(q, len(ladder.lengths) - len(ladder.exclude)))
    return out


def detect_ladder(docs, *, matrix=None, kits: Optional[Sequence] = None,
                  channels: Optional[Sequence[int]] = None,
                  base_order: Optional[str] = None,
                  enough: Optional[int] = None, min_docs: int = 8) -> Dict[str, object]:
    """Which size standard a run used, and in which physical channel.

    Every channel of each document is fitted against every kit.  The run's
    answer is the channel and kit that give a clean fit (:func:`ladder_fit_ok`)
    in the most wells.  ROX bleeds into its neighbours, so a neighbouring
    channel can fit the same ladder just as cleanly at a fraction of the
    height: among channels that fit about as often, the one with the tallest
    bands is the ladder's own.  With *enough*, documents are read in turn
    and the search stops once *min_docs* have been read and one combination
    is clean in *enough* of them (pass a lazy iterable to load no more than
    needed).  Returns ``ladder`` (kit key, or None),
    ``channel``, ``wells``, ``wells_ok``, ``rivals`` (other kits clean in at
    least half as many wells in that channel) and ``table`` (every
    combination's count of clean wells and median band height).
    """
    if kits is None:
        kits = auto_kits()
    ladders = {}
    for k in kits:
        try:
            ladders[k if isinstance(k, str) else k.name] = load_ladder(k) if isinstance(k, str) else k
        except Exception:                                  # noqa: BLE001
            continue                   # a broken user ladder must not stop detection
    tally: Dict[Tuple[int, str], List[Tuple[float, float]]] = {}
    n_docs = 0
    n_ch = 4
    focus = None          # once a combination is clean in 2 wells: those seen clean
    for doc in docs:
        order = base_order or getattr(doc, "base_order", None) or DEFAULT_BASE_ORDER
        _, corr = _channels_for_sizing(doc, order, matrix)
        n_docs += 1
        n_ch = max(n_ch, corr.shape[1])
        for c in (channels or range(1, corr.shape[1] + 1)):
            y = corr[:, c - 1]
            pk = detect_peaks(y, min_prominence_frac=DEFAULT_LADDER_PROMINENCE_FRAC,
                              min_height_frac=DEFAULT_LADDER_HEIGHT_FRAC)
            noise = _baseline_noise(y)
            if not pk or noise <= 0 or max(p.height for p in pk) < MIN_LADDER_SNR * noise:
                continue                       # nothing standing clear of the noise
            fits: Dict[tuple, dict] = {}     # kits with identical tables fit once
            for key, lad in ladders.items():
                if focus is not None and (c, key) not in focus:
                    continue
                if (lad.dye or "").upper() == "LIZ" and not (corr.shape[1] >= 5 and c == 5):
                    continue           # LIZ is the fifth dye: channel 5 of a five-dye instrument
                tkey = (lad.lengths, lad.exclude)
                if tkey not in fits:
                    fits[tkey] = fit_ladder_channel(y, lad, pk)
                r = fits[tkey]
                if r["ok"]:
                    bps = [b for _, b in r["anchors"]]
                    tally.setdefault((c, key), []).append(
                        (r["band_height"], float(r["quality"]["rms_error_spacing"]),
                         float(r["quality"]["unexplained_tall"] or 0.0),
                         len(bps), min(bps), max(bps)))
        if enough and n_docs >= min_docs and tally and \
                max(len(v) for v in tally.values()) >= enough:
            break
        if focus is None and tally and max(len(v) for v in tally.values()) >= 2:
            # the rest of the run only needs to vote among what has fitted
            focus = set(tally)
    table = [{"channel": c, "ladder": k, "wells_ok": len(v),
              "band_height": float(np.median([t[0] for t in v])),
              "rms_spacing": float(np.median([t[1] for t in v])),
              "unexplained": float(np.median([t[2] for t in v])),
              "anchors": float(np.median([t[3] for t in v])),
              "bp_range": (float(np.median([t[4] for t in v])), float(np.median([t[5] for t in v])))}
             for (c, k), v in tally.items()]
    res = {"ladder": None, "channel": None, "wells": n_docs, "wells_ok": 0,
           "table": sorted(table, key=lambda t: (-t["wells_ok"], t["rms_spacing"]))}
    # one or two clean wells out of dozens is chance, not a size standard
    table = [t for t in table if t["wells_ok"] >= min(MIN_DETECT_WELLS, n_docs)]
    if not table:
        return res
    best = max(t["wells_ok"] for t in table)
    near = [t for t in table if t["wells_ok"] >= 0.5 * best]
    channel = max(near, key=lambda t: t["band_height"])["channel"]
    # kits sharing many lengths (ET550-R, ILS600, ILS500) can all fit one
    # ladder -- the aligner bends the curve to absorb 20 vs 25 bp steps.  Of
    # the kits clean in at least half as many wells as the most, the right
    # one leaves the fewest tall peaks unexplained and fits tightest.
    here = [t for t in table if t["channel"] == channel]
    most = max(t["wells_ok"] for t in here)
    contenders = [t for t in here if t["wells_ok"] >= 0.5 * most]
    # Ties (kits nested in one another explain the same peaks equally well)
    # go to the kit that matched more bands, then to the one with fewer bands
    # missing: MapMarker 400 rather than SS400 when 60 and 420 are absent.
    # Stepwise: fewest unexplained peaks; of those, the tightest fits (an
    # evenly spaced kit bent onto another's bands fits several times looser);
    # then the most bands matched, counting one band either way as a tie;
    # then the fewest kit bands; then the tightest fit; then the dye.
    least = min(round(t["unexplained"], 2) for t in contenders)
    group = [t for t in contenders if round(t["unexplained"], 2) == least]
    tight = min(t["rms_spacing"] for t in group)
    group = [t for t in group if t["rms_spacing"] <= 3.0 * max(tight, 0.005)]
    most_anch = max(t["anchors"] for t in group)
    group = [t for t in group if t["anchors"] >= most_anch - 1.0]
    pick = min(group, key=lambda t: (len(ladders[t["ladder"]].lengths), t["rms_spacing"],
                                     not _dye_fits_channel(ladders[t["ladder"]], channel, n_ch)))
    res.update(ladder=pick["ladder"], channel=channel, wells_ok=pick["wells_ok"])
    # a kit that explains the channel as well and fits about as tightly is a
    # real doubt on a weak plate: say so rather than pretend -- unless it
    # lists the same lengths over the bands seen, when the sizes are the same
    lo, hi = pick["bp_range"]

    def unexpl_peaks(t):      # unexplained tall peaks, as a count
        return t["unexplained"] * len(ladders[t["ladder"]].lengths)

    def seen(key):
        return [round(b, 3) for b in ladders[key].lengths if lo - 1e-6 <= b <= hi + 1e-6]

    res["rivals"] = [t["ladder"] for t in contenders if t is not pick
                     and unexpl_peaks(t) <= unexpl_peaks(pick) + 0.5
                     and t["rms_spacing"] <= 2.0 * pick["rms_spacing"]
                     and t["anchors"] >= pick["anchors"] - 1
                     and seen(t["ladder"]) != seen(pick["ladder"])]
    res["same_sizes"] = [t["ladder"] for t in contenders if t is not pick
                         and seen(t["ladder"]) == seen(pick["ladder"])]
    return res


def plate_template(results: Sequence["SizingResult"]) -> Optional[Dict[float, float]]:
    """Median scan of each ladder band over the wells whose own fit is clean:
    where a band of this run is expected (``{bp: scan}``), or None when fewer
    than three wells fit."""
    good = [r for r in results if r.anchors and
            ladder_fit_ok(r.quality, len(r.ladder.lengths) - len(r.ladder.exclude))]
    if len(good) < 3:
        return None
    per: Dict[float, List[int]] = {}
    for r in good:
        for s, b in r.anchors:
            per.setdefault(float(b), []).append(int(s))
    need = max(2, len(good) // 2)
    return {b: float(np.median(v)) for b, v in sorted(per.items()) if len(v) >= need}


# a rescued ladder's bands must be this many times the other peaks in their span
RESCUE_HEIGHT_RATIO = 1.5


def align_to_template(scans: Sequence[int], lengths: Sequence[float],
                      template: Dict[float, float],
                      heights: Optional[Dict[int, float]] = None) -> List[Tuple[int, float]]:
    """Find a well's ladder bands from where the plate's other wells have them.

    Capillaries of one run differ in speed, so the well's band ``b`` is
    looked for at ``a * template[b] + c``: every stretch ``a`` from 0.75 to
    1.33 and shift ``c`` is voted on by the peaks it lines up with bands, the
    best is refined with a quadratic on the matches, and the bands missed are
    grown one by one from the trusted ones (:func:`_extend`).  With *heights*
    ({scan: height}) a fit whose bands are no taller than the other peaks
    between them is refused as a chance match."""
    scans = sorted(int(s) for s in scans)
    lengths = [float(b) for b in lengths]
    tb = np.asarray([b for b in lengths if b in template], float)
    if len(scans) < 3 or tb.size < 4:
        return []
    ts = np.asarray([template[b] for b in tb], float)
    sc = np.asarray(scans, float)
    tol0 = 0.3 * float(np.min(np.diff(ts)))

    def match(pred, tol):
        pairs, used = [], set()
        for j in np.argsort([np.min(np.abs(sc - p)) for p in pred]):
            d = np.abs(sc - pred[j])
            i = int(np.argmin(d))
            if d[i] <= tol and i not in used:
                used.add(i)
                pairs.append((i, int(j)))
        return sorted(pairs, key=lambda p: p[1])

    best, best_key = None, None
    for a in np.linspace(0.75, 1.33, 59):
        # vote: each (peak, band) pair proposes a shift; shifts that line up
        # many bands land in the same bin.  Only the best few bins are tried.
        off = (sc[:, None] - a * ts[None, :]).ravel()
        w = tol0 * a
        bins = np.floor(off / w).astype(int)
        u, cnt = np.unique(bins, return_counts=True)
        if u.size == 0:
            continue
        both = {}
        for b, n in zip(u, cnt):
            both[b] = both.get(b, 0) + n
            both[b + 1] = both.get(b + 1, 0) + n     # a shift straddling two bins
        top = sorted(both, key=both.get, reverse=True)[:4]
        for b in top:
            sel = off[(bins == b) | (bins == b - 1)]
            if sel.size == 0:
                continue
            c = float(np.median(sel))
            pairs = match(a * ts + c, w)
            if len(pairs) < 3:
                continue
            err = sum(abs(sc[i] - (a * ts[j] + c)) for i, j in pairs)
            key = (len(pairs), -err)
            if best_key is None or key > best_key:
                best, best_key = (a, c), key
    if best is None:
        return []
    pred = best[0] * ts + best[1]
    pairs = match(pred, tol0 * best[0])
    # chance check: in a channel crowded with small peaks any template finds
    # a peak near every band.  Real bands of one kit stand above the other
    # peaks between them (dye bleed, pull-up, noise); chance matches are just
    # more of the same peaks.  Measured on cell ID plates: rescued ladders
    # 1.5-7x the other peaks in their span, chance matches on noise 1.0x.
    if heights is not None and pairs:
        hs = np.asarray([heights.get(int(sc[i]), 0.0) for i, _ in pairs], float)
        lo, hi = sc[pairs[0][0]], sc[pairs[-1][0]]
        used = {i for i, _ in pairs}
        rest = [heights.get(int(x), 0.0) for k, x in enumerate(sc) if lo <= x <= hi and k not in used]
        if rest and float(np.median(hs)) < RESCUE_HEIGHT_RATIO * float(np.median(rest)):
            return []
    for _ in range(3):                      # follow the capillary's own curvature
        if len(pairs) < 4:
            break
        coef = np.polyfit(ts[[j for _, j in pairs]], sc[[i for i, _ in pairs]], 2)
        pairs = match(np.polyval(coef, ts), 0.6 * tol0 * best[0])
    # back to indices into the kit's own table, then grow from the trusted
    idx = {b: k for k, b in enumerate(lengths)}
    full = [(i, idx[float(tb[j])]) for i, j in pairs]
    full = sorted(full, key=lambda p: p[1])
    mono = []
    for i, j in full:                       # bands come out in order of length
        if not mono or i > mono[-1][0]:
            mono.append((i, j))
    full = mono
    good = _trusted(scans, lengths, full)
    if len(good) >= 4:
        full = _extend(scans, lengths, good)
    return [(scans[i], lengths[j]) for i, j in full]


# --------------------------------------------------------------------------- #
# per-run band correction
# --------------------------------------------------------------------------- #
# A ladder lot can carry a band that migrates off its listed length: in the
# lab's ET550-R runs from September 2012 on, the 125 bp band sits ~2 bp long
# against its neighbours in nearly every well, while 2011-12 runs agree with
# the table (docs/FRAGMENT_SIZING_QC.md).  The sizing curve passes through
# every anchor, so such a band bends the sizes around it.  A run with enough
# clean wells measures each band against its neighbours; a band that is off
# the same way across the run gets the length the run says it has.
BAND_CORRECT_MIN_WELLS = 8      # clean wells measuring a band
BAND_CORRECT_MIN_BP = 0.75      # smaller offsets are left alone (noise ~0.6 bp)
BAND_CORRECT_MAX_BP = 4.0       # larger ones are a mismatch, not a lot quirk
BAND_CORRECT_MAX_BANDS = 3      # isolated bands only
BAND_CORRECT_AGREE = 0.75       # share of wells whose offset has the same sign


def local_band_offsets(anchors: Sequence[Tuple[int, float]], *, k: int = 3,
                       deg: int = 3) -> Dict[float, float]:
    """For each anchor with at least two anchors on either side: the length a
    local cubic (in scan, through up to *k* anchors per side) predicts for it
    from its neighbours, minus its own length."""
    pts = sorted((float(s), float(b)) for s, b in anchors)
    sc = np.array([p[0] for p in pts])
    bp = np.array([p[1] for p in pts])
    out: Dict[float, float] = {}
    for i in range(2, len(pts) - 2):
        nb = [j for j in range(max(0, i - k), min(len(pts), i + k + 1)) if j != i]
        x = (sc[nb] - sc[i]) / 100.0
        coef = np.polyfit(x, bp[nb], deg)
        out[bp[i]] = float(coef[-1] - bp[i])
    return out


def run_band_corrections(results: Sequence["SizingResult"]) -> Dict[float, float]:
    """``{listed bp: length to use}`` for the bands of this run's ladder that
    sit off their listed length the same way across its clean wells (empty
    when none do, or when the run has too few clean wells to say).

    Bands are corrected one at a time, the worst first, because one band off
    pulls its neighbours' offsets the other way.  The result is kept only if
    it makes the run's wells fit their ladder better."""
    good = [r for r in results if r.anchors and
            ladder_fit_ok(r.quality, len(r.ladder.lengths) - len(r.ladder.exclude))]
    if len(good) < BAND_CORRECT_MIN_WELLS:
        return {}
    corr: Dict[float, float] = {}

    def used(b):
        return corr.get(float(b), float(b))

    for _ in range(BAND_CORRECT_MAX_BANDS):
        per: Dict[float, List[float]] = {}
        for r in good:
            back = {used(b): float(b) for _, b in r.anchors}
            for u, off in local_band_offsets([(s, used(b)) for s, b in r.anchors]).items():
                per.setdefault(back[u], []).append(off)
        best = None
        for b, offs in per.items():
            if len(offs) < BAND_CORRECT_MIN_WELLS:
                continue
            med = float(np.median(offs))
            same = np.mean(np.sign(offs) == np.sign(med))
            if abs(med) >= BAND_CORRECT_MIN_BP and same >= BAND_CORRECT_AGREE \
                    and (best is None or abs(med) > abs(best[1])):
                best = (b, med)
        if best is None:
            break
        new = used(best[0]) + best[1]
        if abs(new - best[0]) > BAND_CORRECT_MAX_BP:
            break
        corr[best[0]] = round(new, 2)
    if not corr:
        return {}

    def run_rms(c):
        vals = []
        for r in good:
            q = ladder_fit_quality([s for s, _ in r.anchors],
                                   [c.get(float(b), float(b)) for _, b in r.anchors])
            if q.get("rms_error_bp") is not None:
                vals.append(q["rms_error_bp"])
        return float(np.median(vals)) if vals else float("inf")

    return corr if run_rms(corr) < run_rms({}) else {}


def apply_band_corrections(res: "SizingResult", corr: Dict[float, float]) -> None:
    """Re-size *res* in place with the corrected band lengths *corr*: its
    anchors, sizing curve, sample lengths, fit quality and band table."""
    if not corr or len(res.anchors) < 2:
        return
    fix = lambda b: corr.get(float(b), float(b))  # noqa: E731
    res.anchors = [(s, fix(b)) for s, b in res.anchors]
    q = ladder_fit_quality([s for s, _ in res.anchors], [b for _, b in res.anchors])
    q["unexplained_tall"] = res.quality.get("unexplained_tall")
    q["bands_reached"] = res.quality.get("bands_reached")
    res.quality = q
    curve = SizeCurve([s for s, _ in res.anchors], [b for _, b in res.anchors])
    rms = q.get("rms_error_bp")
    for r in res.rows:
        r["length_bp"] = round(float(curve.bp([r["scan"]])[0]), 2)
        r["ladder_rms_bp"] = None if rms is None else round(float(rms), 3)
    for r in res.standards:
        u = fix(r["band_bp"])
        r["band_used_bp"] = u
        others = [(s, b) for s, b in res.anchors if s != r["scan"]]
        if len(others) >= 2:
            pred = float(SizeCurve([o[0] for o in others], [o[1] for o in others]).bp([r["scan"]])[0])
            r["predicted_bp"] = round(pred, 2)
            r["residual_bp"] = round(pred - u, 2)
    res.band_corrections = dict(corr)


# --------------------------------------------------------------------------- #
# top-level sizing
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# sample peaks: detection, shape and a quality score
# --------------------------------------------------------------------------- #
DEFAULT_MIN_PEAK_HEIGHT = 50.0   # signal units above the baseline (GP 2.2 lists ~50 RFU up)
SAMPLE_SMOOTH = 9                # Savitzky-Golay window after a 3-point median
BASELINE_WINDOW = 201            # scans: running-minimum baseline
PULLUP_RATIO = 3.0               # another dye this many times taller at the same scan


def baseline_corrected(y: np.ndarray) -> np.ndarray:
    """Smoothed trace minus a running-minimum baseline (itself smoothed):
    the baseline correction every sizing height and area is measured on."""
    from scipy.ndimage import minimum_filter1d, uniform_filter1d
    from scipy.signal import medfilt, savgol_filter
    y = np.asarray(y, float)
    if y.size >= SAMPLE_SMOOTH + 2:
        y = savgol_filter(medfilt(y, 3), SAMPLE_SMOOTH, 2, mode="interp")
    w = min(BASELINE_WINDOW, max(3, y.size // 2 * 2 - 1))
    return y - uniform_filter1d(minimum_filter1d(y, w), w)


def _noise_sd(z: np.ndarray) -> float:
    """Robust noise SD of a baseline-corrected trace (MAD of its quiet part)."""
    z = np.asarray(z, float)
    q = z[np.abs(z) <= np.percentile(np.abs(z), 80)]
    if q.size < 10:
        return 0.0
    return float(1.4826 * np.median(np.abs(q - np.median(q))))


def detect_sample_peaks(y, min_height: float = DEFAULT_MIN_PEAK_HEIGHT,
                        distance: int = 3) -> Tuple[List[dict], np.ndarray, float]:
    """Peaks of one (dye-corrected) channel: (peaks, baseline-corrected
    trace, noise SD).

    A peak must stand ``min_height`` above the running baseline (and five
    noise SDs), with a prominence of at least 80 % of that -- the kind of peak
    list Genetic Profiler 2.2 reports.  Each peak carries height, full width
    at half maximum (scans) and area above the baseline."""
    from scipy.signal import find_peaks, peak_widths
    z = baseline_corrected(y)
    noise = _noise_sd(z)
    thr = max(float(min_height), 5.0 * noise)
    idx, _ = find_peaks(z, height=thr, prominence=0.8 * thr, distance=distance)
    if not len(idx):
        return [], z, noise
    import warnings as _w
    with _w.catch_warnings():
        _w.simplefilter("ignore")
        half = peak_widths(z, idx, rel_height=0.5)
        foot = peak_widths(z, idx, rel_height=0.95)
    peaks = []
    for k, i in enumerate(idx):
        lo, hi = int(np.floor(foot[2][k])), int(np.ceil(foot[3][k])) + 1
        peaks.append({"scan": int(i), "height": float(z[i]),
                      "width": float(half[0][k]),
                      "area": float(np.clip(z[lo:hi], 0, None).sum())})
    return peaks, z, noise


def _flat_top(raw: np.ndarray, scan: int) -> bool:
    """Four or more scans pinned at the channel's maximum around *scan*:
    the detector clipped, so height and area are only lower bounds."""
    seg = raw[max(0, scan - 4):scan + 5]
    if not seg.size:
        return False
    top = float(raw.max())
    return top > 0 and float(seg.max()) >= 0.995 * top and \
        int(np.sum(seg >= 0.995 * top)) >= 4


def peak_quality(*, snr: float, width_ratio: Optional[float], pullup_ratio: float,
                 in_range: bool, flat: bool, ladder_rel: Optional[float]):
    """A 0-100 quality score for one sized peak, and the reasons it is not 100.

    The product of: the ladder fit (relative to the band spacing), signal to
    noise, shape (width against the ladder's at that position), dye pull-up,
    position inside the ladder, and clipping."""
    q, flags = 1.0, []
    if ladder_rel is None:
        q *= 0.5
        flags.append("ladder fit unknown")
    else:
        q *= float(np.clip(1.0 - ladder_rel / MAX_LADDER_RMS_SPACING, 0.0, 1.0))
        if ladder_rel > MAX_LADDER_RMS_SPACING:
            flags.append("poor ladder")
    q *= float(np.clip(np.log(max(snr, 1e-9) / 3.0) / np.log(50.0 / 3.0), 0.05, 1.0))
    if snr < 10:
        flags.append("weak")
    if width_ratio is not None:
        if width_ratio > 2.0:
            q *= 0.5
            flags.append("broad")
        elif width_ratio < 0.4:
            q *= 0.4
            flags.append("spike")
    if pullup_ratio >= PULLUP_RATIO:
        q *= 0.3
        flags.append("possible pull-up")
    if not in_range:
        q *= 0.6
        flags.append("outside ladder")
    if flat:
        q *= 0.7
        flags.append("flat top (saturated)")
    return int(round(100 * q)), flags


def _channels_for_sizing(doc, order: str, matrix):
    """(raw, dye-corrected) arrays of every physical channel of *doc*."""
    from analyzer_core import n_channels
    n = n_channels(doc)
    raw = np.column_stack([physical_channel(doc, c, order) for c in range(1, n + 1)])
    if matrix is not None and not matrix.is_identity() and raw.shape[1] >= matrix.n:
        return raw, matrix.correct(raw)
    return raw, raw


def size_trace(doc, ladder, *, ladder_channel: int = 4,
               sample_channel: int = 2, base_order: Optional[str] = None,
               ladder_prominence_frac: float = DEFAULT_LADDER_PROMINENCE_FRAC,
               ladder_height_frac: float = DEFAULT_LADDER_HEIGHT_FRAC,
               sample_prominence_frac: float = DEFAULT_SAMPLE_PROMINENCE_FRAC,
               sample_height_frac: float = DEFAULT_SAMPLE_HEIGHT_FRAC,
               matrix=None, sample_channels: Optional[Sequence[int]] = None,
               min_height: float = DEFAULT_MIN_PEAK_HEIGHT,
               min_bp: Optional[float] = None,
               template: Optional[Dict[float, float]] = None) -> SizingResult:
    """Size the sample peaks of one loaded trace against *ladder*.

    ``ladder_channel`` / ``sample_channel`` are *physical* channels, 1..N --
    five on a LIZ-capable instrument (see
    :func:`analyzer_core.physical_channel`); ``sample_channels`` sizes several
    dyes at once (every channel but the ladder's, as Genetic Profiler does).
    ``matrix`` (:class:`spectral.SpectralMatrix`) removes dye overlap first,
    so a strong peak in one dye is not sized again as pull-up in the others.
    Every channel is baseline-corrected (:func:`baseline_corrected`) before
    peaks are found and measured.  Sample peaks below ``min_bp`` (default:
    half the smallest band -- the primer/dye blob ahead of the ladder) are
    left out.

    Besides the sized sample rows the result lists the ladder bands it found
    (``standards``: band, scan, height, and how well the curve through the
    other bands predicts each), and every sample row carries width,
    signal-to-noise, a 0-100 ``quality`` and the ``flags`` behind it.
    """
    if isinstance(ladder, str):
        ladder = load_ladder(ladder)
    order = (base_order or getattr(doc, "base_order", None) or DEFAULT_BASE_ORDER)
    warnings: List[str] = []
    chans = [int(c) for c in (sample_channels or [sample_channel])]
    if ladder_channel in chans:
        if len(chans) == 1:
            warnings.append("ladder and sample channel are the same; sizing will "
                            "match the sample's own peaks to the ladder")
        else:
            chans = [c for c in chans if c != ladder_channel]
    well = getattr(doc, "well", "") or Path(str(getattr(doc, "path", ""))).stem
    file_name = Path(str(getattr(doc, "path", ""))).name
    mname = "" if matrix is None or matrix.is_identity() else matrix.name

    raw, corr = _channels_for_sizing(doc, order, matrix)
    for c in [ladder_channel] + chans:
        if not 1 <= c <= raw.shape[1]:
            raise ValueError(f"Channel must be 1..{raw.shape[1]}, got {c!r}")
    ladder_trace = corr[:, ladder_channel - 1]

    def empty(why):
        warnings.append(why)
        return SizingResult(ladder=ladder, ladder_channel=ladder_channel,
                            sample_channel=chans[0], base_order=order,
                            anchors=[], quality=ladder_fit_quality([], []),
                            rows=[], well=str(well), file=file_name,
                            warnings=warnings, matrix=mname, sample_channels=chans)

    ladder_peaks = detect_peaks(ladder_trace,
                                min_prominence_frac=ladder_prominence_frac,
                                min_height_frac=ladder_height_frac)
    # Every threshold in ``detect_peaks`` is a fraction of the largest sample in
    # its own channel, which leaves nothing to say "this channel holds no
    # ladder": require the largest band to clear the baseline noise by a wide,
    # measured margin (MegaBACE: noise ~6-21 against bands of 800-2100).
    noise = _baseline_noise(ladder_trace)
    strongest = max((p.height for p in ladder_peaks), default=0.0)
    if not ladder_peaks or noise <= 0.0 or strongest < MIN_LADDER_SNR * noise:
        top = f"{strongest:.0f}" if noise > 0 else "n/a"
        return empty(f"no size standard in the ladder channel (largest peak {top} "
                     f"against noise {noise:.1f}, needs {MIN_LADDER_SNR:g}x); no "
                     "lengths were sized")

    all_peaks = ladder_peaks
    ladder_peaks = ladder_candidates(ladder_peaks, len(ladder.lengths))
    matched = align_peaks([p.scan for p in ladder_peaks], ladder.lengths,
                          n_scans=len(ladder_trace))
    excluded = lambda b: any(abs(b - x) < 1e-6 for x in ladder.exclude)  # noqa: E731
    anchors = [(sc, b) for sc, b in matched if not excluded(b)]
    quality = ladder_fit_quality([sc for sc, _ in anchors], [b for _, b in anchors])
    quality["unexplained_tall"] = unexplained_tall(all_peaks, anchors, len(ladder.lengths))
    kit_bands = [b for b in ladder.lengths if not excluded(b)]
    quality["bands_reached"] = bands_reached(anchors, kit_bands, len(ladder_trace))
    if template is not None and not ladder_fit_ok(quality, len(ladder.lengths)):
        # the well's own fit failed: let the plate's other wells say where the
        # bands should be, and look for them among every peak of the channel
        # (a ladder buried under pull-up is shorter than the junk above it)
        rescued = align_to_template([p.scan for p in all_peaks], ladder.lengths, template,
                                    {p.scan: p.height for p in all_peaks})
        r_anchors = [(sc, b) for sc, b in rescued if not excluded(b)]
        r_quality = ladder_fit_quality([sc for sc, _ in r_anchors], [b for _, b in r_anchors])
        r_quality["bands_reached"] = bands_reached(r_anchors, kit_bands, len(ladder_trace))
        if ladder_fit_ok(r_quality, len(ladder.lengths)):
            matched, anchors, quality = rescued, r_anchors, r_quality
            warnings.append("ladder found with help from the plate's other wells "
                            "(its own peaks alone did not give a clean fit)")

    if len(anchors) < MIN_LADDER_ANCHORS:
        warnings.append(
            f"only {len(anchors)} ladder peak(s) matched the {ladder.name} "
            f"table ({len(ladder_peaks)} peak(s) found); results are not "
            "trustworthy")
    if ladder_peaks and anchors:
        # bands of the table found, not peaks of the channel matched: primer
        # and dye peaks ahead of the ladder are expected and harmless
        reach = quality.get("bands_reached") or len(ladder.lengths)
        if len(anchors) / float(min(reach, len(ladder.lengths))) < 0.6:
            warnings.append(
                f"only {len(anchors)} of the {len(ladder.lengths)} ladder bands "
                "were found; check the ladder channel and the fragment lengths")
    rel = quality.get("rms_error_spacing")
    if rel is not None and rel > MAX_LADDER_RMS_SPACING:
        warnings.append(
            f"poor ladder fit: each band is predicted from the others to only "
            f"{quality['rms_error_bp']:.1f} bp RMS ({rel:.0%} of a band spacing); "
            "the lengths are not reliable -- is there a ladder in this well?")
    unexpl = quality.get("unexplained_tall")
    if unexpl is not None and unexpl > MAX_UNEXPLAINED:
        warnings.append(
            f"the ladder channel holds {unexpl * len(ladder.lengths):.0f} tall peaks "
            f"the {ladder.name} table does not explain; the bands may be sample "
            "peaks that happen to fit -- is this the ladder channel and the right kit?")

    curve = None
    if len(anchors) >= 2:
        curve = SizeCurve([sc for sc, _ in anchors], [b for _, b in anchors])
    elif anchors:
        warnings.append("not enough anchors to build a sizing curve")

    # the ladder bands: position, height, and the curve through the others
    lz = baseline_corrected(ladder_trace)
    std_rows, widths = [], []
    if matched:
        from scipy.signal import peak_widths
        import warnings as _w
        with _w.catch_warnings():
            _w.simplefilter("ignore")     # a clipped band can have no half-width
            w_all = peak_widths(lz, [sc for sc, _ in matched], rel_height=0.5)[0]
        for (sc, b), wd in zip(matched, w_all):
            used = not excluded(b)
            others = [(s2, b2) for s2, b2 in anchors if s2 != sc]
            pred = None
            if len(others) >= 2:
                pred = float(SizeCurve([o[0] for o in others],
                                       [o[1] for o in others]).bp([sc])[0])
            std_rows.append({
                "file": file_name, "well": well, "kind": "standard",
                "channel": ladder_channel,
                "dye": channel_dye(doc, ladder_channel) or ladder.dye,
                "band_bp": b, "band_used_bp": b, "scan": int(sc),
                "height": round(float(lz[sc]), 1),
                "width": round(float(wd), 2),
                "predicted_bp": None if pred is None else round(pred, 2),
                "residual_bp": None if pred is None else round(pred - b, 2),
                "used_in_curve": used})
            if used:
                widths.append((sc, float(wd)))

    def expected_width(scan):
        if len(widths) >= 3:
            k = np.polyfit([w[0] for w in widths], [w[1] for w in widths], 1)
            return max(float(np.polyval(k, scan)), 1.0)
        return float(np.median([w[1] for w in widths])) if widths else None

    if min_bp is None:
        min_bp = 0.5 * float(ladder.lengths[0])
    rms = quality.get("rms_error_bp")
    ladder_rms = None if rms is None else round(float(rms), 3)
    ladder_ok = not warnings
    zall = {c: baseline_corrected(corr[:, c - 1]) for c in range(1, corr.shape[1] + 1)}
    rows: List[dict] = []
    for ch in chans:
        if curve is None:
            break
        peaks, z, nz = detect_sample_peaks(corr[:, ch - 1], min_height=min_height)
        if not peaks:
            continue
        base_letter = order[ch - 1] if ch <= len(order) else ""
        dye = channel_dye(doc, ch)
        scans = [p["scan"] for p in peaks]
        lengths_bp = curve.bp(scans)
        in_range = curve.in_range(scans)
        for pk, bp, rng in zip(peaks, lengths_bp, in_range):
            if bp < min_bp:
                continue                      # the primer / dye blob
            sc = pk["scan"]
            other = max((float(zall[k][max(0, sc - 2):sc + 3].max())
                         for k in zall if k != ch), default=0.0)
            ew = expected_width(sc)
            snr = pk["height"] / nz if nz > 0 else float("inf")
            score, flags = peak_quality(
                snr=snr, width_ratio=None if ew is None else pk["width"] / ew,
                pullup_ratio=other / max(pk["height"], 1e-9), in_range=bool(rng),
                flat=_flat_top(raw[:, ch - 1], sc), ladder_rel=rel)
            rows.append({
                "file": file_name,
                "well": well,
                "kind": "sample",
                "channel": ch,
                "base": base_letter,
                "dye": dye,
                "scan": sc,
                "length_bp": round(float(bp), 2),
                "height": round(pk["height"], 1),
                "width": round(pk["width"], 2),
                "area": round(pk["area"], 1),
                "snr": round(min(snr, 9999.0), 1),
                "quality": score,
                "flags": "; ".join(flags),
                "in_range": bool(rng),
                "ladder_rms_bp": ladder_rms,
                "ladder_ok": ladder_ok,
            })

    return SizingResult(ladder=ladder, ladder_channel=ladder_channel,
                        sample_channel=chans[0], base_order=order,
                        anchors=anchors, quality=quality, rows=rows,
                        well=str(well), file=file_name, warnings=warnings,
                        standards=std_rows, matrix=mname, sample_channels=chans)


def table_rows(results: Sequence[SizingResult]) -> List[dict]:
    """One table for a set of sized wells: per well its ladder bands
    (``kind`` = standard: band, scan, height, predicted bp, residual), then
    its sized sample peaks (``kind`` = sample: channel, dye, scan, size,
    height, width, area, S/N, quality 0-100, flags)."""
    out: List[dict] = []
    for res in results:
        lad = res.ladder.name
        for r in res.standards:
            out.append({"ladder": lad, "matrix": res.matrix, **r})
        for r in res.rows:
            out.append({"ladder": lad, "matrix": res.matrix, **r})
    return out


def plate_matrix_for(paths, base_order: Optional[str] = None, prefer: str = "auto",
                     cache: Optional[dict] = None):
    """The dye-overlap matrix for the run folder of *paths* (see
    :func:`spectral.plate_matrix`): a run's own Genetic Profiler calibration,
    else measured from the single-dye wells of the whole folder."""
    import spectral
    from analyzer_core import SUPPORTED_EXTS, load_trace
    if prefer == "none":
        return spectral.identity()
    folder = Path(paths[0]).parent
    key = (str(folder.resolve()), prefer, base_order)
    if cache is not None and key in cache:
        return cache[key]
    chans = {}
    if prefer in ("auto", "plate"):
        for f in sorted(folder.iterdir()):
            if f.suffix in SUPPORTED_EXTS and f.is_file():
                try:
                    d = load_trace(f, base_order=base_order) if base_order else load_trace(f)
                except Exception:                          # noqa: BLE001
                    continue
                if getattr(d, "channels", None) is not None:
                    chans[f.stem] = np.asarray(d.channels, float)
    m = spectral.plate_matrix(folder, chans, prefer=prefer)
    if cache is not None:
        cache[key] = m
    return m


# --------------------------------------------------------------------------- #
# the run's ladder: found once per run, cached, used for any of its wells
# --------------------------------------------------------------------------- #
# Bump when a change to detection, rescue or band correction would give a
# different profile, so cached profiles from older code are not reused.
LADDER_ENGINE_VERSION = 8


@dataclass
class RunLadder:
    """What a run's size standard is, worked out from all of its wells.

    ``ladder``/``channel``: the kit and physical channel used for every well.
    ``status``: "found" (detected), "given" (named by the user), "uncertain"
    (detected, but another kit fits nearly as well -- ``rivals``),
    "fallback" (a named kit whose channel was not recognised: its usual
    channel is used) or "none" (no size standard recognised; an ET550-R guess
    in channel 1 is used and every size is unreliable).  ``template`` is where
    the run's clean wells have each band (for rescuing weak wells) and
    ``corrections`` the band lengths this run's ladder lot needs.
    """
    ladder: Ladder
    channel: int
    status: str
    note: str = ""
    rivals: List[str] = field(default_factory=list)
    same_sizes: List[str] = field(default_factory=list)
    wells: int = 0                   # wells in the run
    wells_tested: int = 0            # wells detection looked at
    wells_clean: int = 0             # wells whose own ladder fit is clean
    template: Optional[Dict[float, float]] = None
    corrections: Dict[float, float] = field(default_factory=dict)
    from_cache: bool = False

    @property
    def reliable(self) -> bool:
        return self.status in ("found", "given", "uncertain") and self.wells_clean > 0

    def to_json(self) -> dict:
        lad = self.ladder
        return {"engine": LADDER_ENGINE_VERSION,
                "ladder": {"name": lad.name, "lengths": list(lad.lengths), "dye": lad.dye,
                           "exclude": list(lad.exclude), "source": lad.source,
                           "notes": lad.notes, "vendor": lad.vendor, "color": lad.color,
                           "default_channel": lad.default_channel},
                "channel": self.channel, "status": self.status, "note": self.note,
                "rivals": self.rivals, "same_sizes": self.same_sizes, "wells": self.wells,
                "wells_tested": self.wells_tested, "wells_clean": self.wells_clean,
                "template": None if self.template is None else
                [[b, sc] for b, sc in sorted(self.template.items())],
                "corrections": [[b, u] for b, u in sorted(self.corrections.items())]}

    @classmethod
    def from_json(cls, d: dict) -> "RunLadder":
        lad = Ladder(**{k: (tuple(v) if k in ("lengths", "exclude") else v)
                        for k, v in d["ladder"].items()})
        return cls(ladder=lad, channel=int(d["channel"]), status=d["status"],
                   note=d.get("note", ""), rivals=list(d.get("rivals", [])),
                   same_sizes=list(d.get("same_sizes", [])), wells=int(d.get("wells", 0)),
                   wells_tested=int(d.get("wells_tested", 0)),
                   wells_clean=int(d.get("wells_clean", 0)),
                   template=None if d.get("template") is None else
                   {float(b): float(sc) for b, sc in d["template"]},
                   corrections={float(b): float(u) for b, u in d.get("corrections", [])},
                   from_cache=True)


def default_cache_dir() -> Path:
    import os
    base = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    return Path(base) / "limoncello-scorer" / "ladders"


def _profile_key(paths, ladder, ladder_channel, matrix, base_order, plate_rescue,
                 band_correct) -> Optional[str]:
    """A cache key that changes when any file of the run, or any setting
    that shapes the profile, changes.  None when the files are not on disk."""
    import hashlib
    h = hashlib.sha1()
    for p in sorted(Path(x) for x in paths):
        try:
            st = p.stat()
        except OSError:
            return None
        h.update(f"{p.resolve()}|{st.st_size}|{st.st_mtime_ns}\n".encode())
    lad = "auto" if ladder is None or isinstance(ladder, str) and ladder.lower() == "auto" \
        else repr((ladder.name, ladder.lengths, ladder.exclude)) if isinstance(ladder, Ladder) \
        else str(ladder)
    mat = matrix if isinstance(matrix, str) else getattr(matrix, "name", repr(matrix))
    h.update(repr((LADDER_ENGINE_VERSION, lad, ladder_channel, mat, base_order,
                   plate_rescue, band_correct)).encode())
    return h.hexdigest()


def _spread(keys: list) -> list:
    """*keys* reordered to spread over the plate: A01, the middle, then the
    gaps -- an unlucky first few wells (empty, overloaded) must not decide."""
    order, step = [], len(keys)
    while step >= 1 and len(order) < len(keys):
        order += [k for k in keys[::step] if k not in order]
        step //= 2
    return order


def find_run_ladder(paths, ladder="auto", *, ladder_channel: Optional[int] = None,
                    base_order: Optional[str] = None, matrix=None, load=None,
                    docs: Optional[dict] = None, detect_sample: int = 8,
                    plate_rescue: bool = True, band_correct: bool = True,
                    cache_dir="default", on_error=None) -> RunLadder:
    """Work out one run's size standard from its wells (*paths*: every trace
    file of the run, so a selection of a few wells is sized the same way as
    the whole plate).

    1. **Kit and channel** (:func:`detect_ladder`): unless both are given,
       every kit (*ladder* "auto") or the named one is tried in every channel
       (or the named one) of wells spread over the plate, until one
       combination is clean in several wells.
    2. **Template**: where the clean wells have each band, for rescuing weak
       wells (:func:`plate_template`).
    3. **Band corrections** for this ladder lot (:func:`run_band_corrections`).

    The result is cached under *cache_dir* (None: no cache), keyed to the
    files' sizes and times and to these settings."""
    from analyzer_core import load_trace
    auto_kit = ladder is None or (isinstance(ladder, str) and ladder.lower() == "auto")
    if not auto_kit and isinstance(ladder, str):
        ladder = load_ladder(ladder)
    paths = [Path(p) for p in paths]
    cdir = default_cache_dir() if cache_dir == "default" else (Path(cache_dir) if cache_dir else None)
    key = _profile_key(paths, None if auto_kit else ladder, ladder_channel, matrix,
                       base_order, plate_rescue, band_correct) if cdir else None
    if key:
        try:
            d = json.loads((cdir / f"{key}.json").read_text(encoding="utf-8"))
            if d.get("engine") == LADDER_ENGINE_VERSION:
                return RunLadder.from_json(d)
        except (OSError, ValueError, KeyError, TypeError):
            pass
    if load is None:
        def load(p):
            return load_trace(Path(p), base_order=base_order) if base_order \
                else load_trace(Path(p))
    docs = {} if docs is None else docs

    def get(p):
        if p not in docs:
            docs[p] = load(p)
        return docs[p]

    def lazy(keys):
        for k in keys:
            try:
                yield get(k)
            except Exception as exc:                       # noqa: BLE001
                if on_error:
                    on_error(k, exc)

    order = _spread(list(paths))
    run_lad, lch, status, note = ladder, ladder_channel, "given", ""
    rivals, same, tested = [], [], 0
    if auto_kit or lch is None:
        det = detect_ladder(lazy(order[:4 * detect_sample]), matrix=matrix,
                            base_order=base_order,
                            kits=auto_kits() if auto_kit else (ladder,),
                            channels=None if lch is None else [lch],
                            enough=6, min_docs=detect_sample)
        if det["table"] and det["wells_ok"] < 6 and len(order) > det["wells"]:
            # a weak plate: something fitted, but in too few of the test
            # wells to name it (or to tell the ladder's channel from the one
            # its dye bleeds into) -- look at the rest of the run first
            det = detect_ladder(lazy(order), matrix=matrix, base_order=base_order,
                                kits=auto_kits() if auto_kit else (ladder,),
                                channels=None if lch is None else [lch],
                                enough=6, min_docs=detect_sample)
        tested = det["wells"]
        if det["ladder"] is None and not auto_kit:
            lch = lch or ladder.default_channel or 4
            status = "fallback"
            note = (f"ladder channel not recognised automatically in {tested} "
                    f"test wells; using channel {lch}, the {ladder.name}'s usual one")
        elif det["ladder"] is None:
            run_lad = load_ladder("et550r")
            lch = lch or run_lad.default_channel or 1
            status = "none"
            note = (f"no size standard recognised in {tested} test wells of this run "
                    f"(tried {len(auto_kits())} kits in every channel); lengths below rest on "
                    "an ET550-R guess in channel 1 and are not reliable")
        else:
            run_lad = load_ladder(det["ladder"]) if auto_kit else ladder
            lch = det["channel"]
            rivals = list(det.get("rivals") or [])
            same = list(det.get("same_sizes") or [])
            status = "uncertain" if rivals else "found"
            note = (f"ladder {run_lad.name} in channel {lch}, found automatically "
                    f"(clean in {det['wells_ok']} of {tested} test wells)")
            if rivals:
                note += ("; kit uncertain: " + ", ".join(load_ladder(k).name for k in rivals)
                         + " fits nearly as well -- name the kit if you know it")
    lch = lch or run_lad.default_channel or 4

    # every well's own ladder fit: the clean ones make the template and
    # measure the band corrections
    from types import SimpleNamespace
    fits = []
    for doc in lazy(order):
        try:
            order_ = base_order or getattr(doc, "base_order", None) or DEFAULT_BASE_ORDER
            _, corr = _channels_for_sizing(doc, order_, matrix)
            if not 1 <= lch <= corr.shape[1]:
                continue
            f = fit_ladder_channel(corr[:, lch - 1], run_lad)
            fits.append(SimpleNamespace(ladder=run_lad, anchors=f["anchors"],
                                        quality=f["quality"]))
        except Exception:                                  # noqa: BLE001
            continue
    nb = len(run_lad.lengths) - len(run_lad.exclude)
    clean = [f for f in fits if f.anchors and ladder_fit_ok(f.quality, nb)]
    tmpl = plate_template(fits) if plate_rescue else None
    corrections = run_band_corrections(fits) if band_correct and status != "none" else {}
    if corrections:
        note += ("; band lengths corrected for this run's ladder lot: " +
                 ", ".join(f"{b:g} -> {u:g} bp" for b, u in sorted(corrections.items())))
    prof = RunLadder(ladder=run_lad, channel=int(lch), status=status, note=note,
                     rivals=rivals, same_sizes=same, wells=len(paths), wells_tested=tested,
                     wells_clean=len(clean), template=tmpl, corrections=corrections)
    if key:
        try:
            cdir.mkdir(parents=True, exist_ok=True)
            tmp = cdir / f"{key}.json.tmp"
            tmp.write_text(json.dumps(prof.to_json()), encoding="utf-8")
            tmp.replace(cdir / f"{key}.json")
        except OSError:
            pass
    return prof


def _run_files(folder: Path, given: List[Path]) -> List[Path]:
    """Every trace file of the run *given* belongs to (its folder), or just
    *given* when they are not all real files there."""
    from analyzer_core import SUPPORTED_EXTS
    try:
        if not all(p.is_file() for p in given):
            return given
        exts = {p.suffix.lower() for p in given}
        allf = sorted(f for f in folder.iterdir()
                      if f.is_file() and f.suffix.lower() in exts
                      and f.suffix in SUPPORTED_EXTS + tuple(e.lower() for e in SUPPORTED_EXTS))
        have = {f.resolve() for f in allf}
        return allf if all(p.resolve() in have for p in given) else given
    except OSError:
        return given


def size_files(paths, ladder="auto", *, ladder_channel: Optional[int] = None,
               sample_channels: Optional[Sequence[int]] = None,
               base_order: Optional[str] = None, matrix="auto",
               min_height: float = DEFAULT_MIN_PEAK_HEIGHT, load=None,
               on_error=None, plate_rescue: bool = True,
               detect_sample: int = 8, band_correct: bool = True,
               whole_run: bool = True, cache_dir="default",
               matrix_cache: Optional[dict] = None,
               profiles: Optional[dict] = None) -> List[SizingResult]:
    """Size many files, run folder by run folder, each run with its own
    dye-overlap matrix (*matrix*: "auto" / "run" / "plate" / "none", or a
    :class:`spectral.SpectralMatrix`).  *sample_channels* None sizes every
    channel but the ladder's.  *load* (path -> document) defaults to
    :func:`analyzer_core.load_trace`.

    Each run's ladder is worked out once (:func:`find_run_ladder`): *ladder*
    "auto" (or None) finds the kit, a *ladder_channel* of None finds the
    channel.  With *whole_run* the run's profile comes from every trace file
    in its folder, so sizing a few selected wells gives the same answer as
    sizing the plate.  With *plate_rescue*, a well whose own ladder fit fails
    is fitted again guided by where the run's clean wells have their bands;
    with *band_correct*, a band that sits off its listed length the same way
    across the run gets the length the run measures.  *profiles*, if given,
    receives each run folder's :class:`RunLadder`."""
    from analyzer_core import load_trace, n_channels
    if load is None:
        def load(p):
            return load_trace(Path(p), base_order=base_order) if base_order \
                else load_trace(Path(p))
    groups: Dict[str, list] = {}
    for p in paths:
        groups.setdefault(str(Path(p).parent), []).append(Path(p))
    cache: dict = {} if matrix_cache is None else matrix_cache
    out: List[SizingResult] = []
    for folder, ps in groups.items():
        run = _run_files(Path(folder), ps) if whole_run else ps
        if isinstance(matrix, str):
            try:
                m = plate_matrix_for(run, base_order=base_order, prefer=matrix,
                                     cache=cache)
            except Exception:                              # noqa: BLE001
                m = None
        else:
            m = matrix
        docs: Dict[Path, object] = {}
        prof = find_run_ladder(run, ladder, ladder_channel=ladder_channel,
                               base_order=base_order, matrix=m, load=load, docs=docs,
                               detect_sample=detect_sample, plate_rescue=plate_rescue,
                               band_correct=band_correct, cache_dir=cache_dir,
                               on_error=on_error)
        if profiles is not None:
            profiles[folder] = prof
        lch = prof.channel
        for p in ps:
            try:
                doc = docs[p] if p in docs else load(p)
                chans = list(sample_channels) if sample_channels else \
                    [c for c in range(1, n_channels(doc) + 1) if c != lch]
                r = size_trace(doc, prof.ladder, ladder_channel=lch, sample_channels=chans,
                               base_order=base_order, matrix=m, min_height=min_height,
                               template=prof.template if plate_rescue else None)
                if prof.corrections:
                    apply_band_corrections(r, prof.corrections)
                r.detection = prof.note
                r.path = str(p)
                if prof.status in ("none", "fallback"):
                    r.warnings.insert(0, prof.note)
                out.append(r)
            except Exception as exc:                       # noqa: BLE001
                if on_error:
                    on_error(p, exc)
    return out
