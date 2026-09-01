#!/usr/bin/env python3
"""Zero-dependency test runner (no pytest needed): runs every `test_*` function of every
`tests/test_*.py` module.

Usage: python3 tests/run_tests.py        -> prints one line per test, exits 1 if any failed.
"""
import glob, importlib, os, sys, time, traceback

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
HERE = os.path.join(ROOT, "tests")
sys.path.insert(0, HERE)


def main():
    mods = sorted(os.path.basename(p)[:-3] for p in glob.glob(os.path.join(HERE, "test_*.py")))
    total = failed = 0
    for mname in mods:
        mod = importlib.import_module(mname)
        for name in sorted(n for n in dir(mod) if n.startswith("test_")):
            total += 1
            t0 = time.time()
            try:
                getattr(mod, name)()
                print(f"PASS {mname}.{name} ({time.time()-t0:.2f}s)")
            except Exception:
                failed += 1
                print(f"FAIL {mname}.{name} ({time.time()-t0:.2f}s)")
                traceback.print_exc()
    print(f"{total-failed}/{total} tests passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
