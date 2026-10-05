"""Make the flat-layout modules importable when pytest runs from the repo.

The modules here are flat, not a package: ``genotyping``, ``scorer`` and
``sequence_analyzer`` are siblings of this directory, and ``analyzer_core`` is
declared in ``pyproject.toml`` as one of this project's ``py-modules``, so it
is a sibling too.  Without this the imports resolve only by luck of the working
directory.

The repository root is added as well, because the basecaller packages
(``cimarron_basecaller``) and the sibling copies of these same modules live
there.  The subproject root stays in front so its own modules keep shadowing the
top-level copies, as they always have.
"""
import sys
from pathlib import Path

SUBPROJECT = Path(__file__).resolve().parent.parent
REPO_ROOT = SUBPROJECT.parent

for path in (SUBPROJECT, REPO_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))
# Subproject first: it is the later insert, so it ends up at the front.
if sys.path[0] != str(SUBPROJECT):
    sys.path.remove(str(SUBPROJECT))
    sys.path.insert(0, str(SUBPROJECT))
