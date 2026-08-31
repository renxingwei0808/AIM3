#!/usr/bin/env python3
"""Zero-dependency test runner (no pytest needed): runs every `test_*` function of tests/test_core.py.

Usage: python3 tests/run_tests.py        -> prints one line per test, exits 1 on the first failure count.
"""
import os, sys, time, traceback

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "tests"))
import test_core

def main():
    names = sorted(n for n in dir(test_core) if n.startswith("test_"))
    failed = 0
    for name in names:
        t0 = time.time()
        try:
            getattr(test_core, name)()
            print(f"PASS {name} ({time.time()-t0:.2f}s)")
        except Exception:
            failed += 1
            print(f"FAIL {name} ({time.time()-t0:.2f}s)")
            traceback.print_exc()
    print(f"{len(names)-failed}/{len(names)} tests passed")
    return 1 if failed else 0

if __name__ == "__main__":
    sys.exit(main())
