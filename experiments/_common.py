"""import shim: makes `legacy/*` and `src/aim3` importable from any experiment directory."""
import os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
for p in (os.path.join(ROOT, "legacy"), os.path.join(ROOT, "src")):
    if p not in sys.path:
        sys.path.insert(0, p)
RESULTS = os.path.join(ROOT, "results")
