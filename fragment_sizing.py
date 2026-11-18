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
    sample_channels: List[int] = field(default_factory=list)

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
    inf = float("inf")
    dp = np.full((n + 1, k + 1), inf)
    back = np.zeros((n + 1, k + 1), dtype=np.int8)  # 0 match, 1 skip scan, 2 skip length
    dp[0, 0] = 0.0
    for i in range(n + 1):
        for j in range(k + 1):
            if i == 0 and j == 0:
                continue
            best, choice = inf, 0
            if i and j:
                cost = dp[i - 1, j - 1] + abs(pred[j - 1] - scans[i - 1])
                if cost < best:
                    best, choice = cost, 0
            if i and dp[i - 1, j] + gap < best:
                best, choice = dp[i - 1, j] + gap, 1
            if j and dp[i, j - 1] + gap < best:
                best, choice = dp[i, j - 1] + gap, 2
            dp[i, j], back[i, j] = best, choice
    pairs = []
    i, j = n, k
    while i or j:
        choice = back[i, j]
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
                     pairs: List[Tuple[int, int]]) -> float:
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
    best = float("inf")
    for u in (np.log(bp), bp):
        deg = 2 if len(pairs) >= 4 else 1
        coef = np.polyfit(u, v, deg)
        resid = v - np.polyval(coef, u)
        best = min(best, float(np.sqrt(np.mean(resid ** 2))) / scans_per_bp)
    return best + (len(lengths) - len(pairs))


def align_peaks(scans: Sequence[int], lengths: Sequence[float],
                iterations: int = 5) -> List[Tuple[int, float]]:
    """Match ladder peak scans to the kit's known lengths.

    Bootstraps a migration model (scan linear in log(bp), or in bp), aligns
    with DP, refits the model on the matches, and repeats -- the iteration is what absorbs the curvature
    a single global fit cannot.  One stray peak at either end (a primer peak,
    dye bleed from a sample channel) would wreck a bootstrap taken from the
    outermost peaks, so several end pairings are tried -- each of the first and
    last few peaks against each of the first and last few lengths -- and the
    alignment whose anchors sit best on a smooth curve wins; a repair pass then
    re-aligns from the anchors that alignment got clearly right.
    Returns ``(scan, length)`` anchors.
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
        score = _alignment_score(scans, lengths, pairs)
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
               min_bp: Optional[float] = None) -> SizingResult:
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

    ladder_peaks = ladder_candidates(ladder_peaks, len(ladder.lengths))
    matched = align_peaks([p.scan for p in ladder_peaks], ladder.lengths)
    excluded = lambda b: any(abs(b - x) < 1e-6 for x in ladder.exclude)  # noqa: E731
    anchors = [(sc, b) for sc, b in matched if not excluded(b)]
    quality = ladder_fit_quality([sc for sc, _ in anchors], [b for _, b in anchors])

    if len(anchors) < MIN_LADDER_ANCHORS:
        warnings.append(
            f"only {len(anchors)} ladder peak(s) matched the {ladder.name} "
            f"table ({len(ladder_peaks)} peak(s) found); results are not "
            "trustworthy")
    if ladder_peaks and anchors:
        # bands of the table found, not peaks of the channel matched: primer
        # and dye peaks ahead of the ladder are expected and harmless
        if len(anchors) / float(len(ladder.lengths)) < 0.6:
            warnings.append(
                f"only {len(anchors)} of the {len(ladder.lengths)} ladder bands "
                "were found; check the ladder channel and the fragment lengths")
    rel = quality.get("rms_error_spacing")
    if rel is not None and rel > MAX_LADDER_RMS_SPACING:
        warnings.append(
            f"poor ladder fit: each band is predicted from the others to only "
            f"{quality['rms_error_bp']:.1f} bp RMS ({rel:.0%} of a band spacing); "
            "the lengths are not reliable -- is there a ladder in this well?")

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
                "band_bp": b, "scan": int(sc), "height": round(float(lz[sc]), 1),
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


def size_files(paths, ladder, *, ladder_channel: Optional[int] = None,
               sample_channels: Optional[Sequence[int]] = None,
               base_order: Optional[str] = None, matrix="auto",
               min_height: float = DEFAULT_MIN_PEAK_HEIGHT, load=None,
               on_error=None) -> List[SizingResult]:
    """Size many files, run folder by run folder, each run with its own
    dye-overlap matrix (*matrix*: "auto" / "run" / "plate" / "none", or a
    :class:`spectral.SpectralMatrix`).  *sample_channels* None sizes every
    channel but the ladder's.  *load* (path -> document) defaults to
    :func:`analyzer_core.load_trace`."""
    from analyzer_core import load_trace, n_channels
    if isinstance(ladder, str):
        ladder = load_ladder(ladder)
    lch = ladder_channel or ladder.default_channel or 4
    if load is None:
        def load(p):
            return load_trace(Path(p), base_order=base_order) if base_order \
                else load_trace(Path(p))
    groups: Dict[str, list] = {}
    for p in paths:
        groups.setdefault(str(Path(p).parent), []).append(Path(p))
    cache: dict = {}
    out: List[SizingResult] = []
    for folder, ps in groups.items():
        if isinstance(matrix, str):
            try:
                m = plate_matrix_for(ps, base_order=base_order, prefer=matrix,
                                     cache=cache)
            except Exception:                              # noqa: BLE001
                m = None
        else:
            m = matrix
        for p in ps:
            try:
                doc = load(p)
                chans = list(sample_channels) if sample_channels else \
                    [c for c in range(1, n_channels(doc) + 1) if c != lch]
                out.append(size_trace(doc, ladder, ladder_channel=lch,
                                      sample_channels=chans, base_order=base_order,
                                      matrix=m, min_height=min_height))
            except Exception as exc:                       # noqa: BLE001
                if on_error:
                    on_error(p, exc)
    return out
