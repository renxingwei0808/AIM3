#!/usr/bin/env bash
# Fast end-to-end check (~1 min): reference impl self-test, structures, E02, E03 small, E06 small, E08, unit tests.
# NOTE: the `| head` pipes are followed by `|| true` because `set -o pipefail` otherwise turns the
# SIGPIPE/BrokenPipeError of the producer into a script failure.
set -euo pipefail
cd "$(dirname "$0")/.."
python3 -c "import sys; sys.path.insert(0,'src'); import random; from aim3 import AIM3; A=AIM3(11,(2,3,4)); print('toy selftest', A.selftest(random.Random(1),3))"
python3 -c "import sys; sys.path.insert(0,'src'); import random; from aim3.structures import SboxStructures; print('structures n=128', SboxStructures(128,11).check_random(random.Random(2),3))"
python3 experiments/E02_official_baseline/table4_reproduce.py | head -4 || true
python3 experiments/E03_syzygy_and_degree_fall/sbox_degree_fall.py 9 4 --no-mac4
python3 experiments/E06_kernel_bilinear_model/kernel_bilinear_toy.py 7 2
python3 experiments/E08_parameter_audit/exponent_audit.py 128 | head -2 || true
python3 tests/run_tests.py
echo "smoke test OK"
