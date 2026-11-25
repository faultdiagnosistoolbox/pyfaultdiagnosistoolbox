"""Minimal example: MSO() crashes on an unknown that occurs in no equation.

The crash needs both
  * a zero column in X, i.e. an unknown that occurs in no equation, and
  * no overdetermined part.

    X = [[1, 0]]            # e1 contains x1; x2 occurs nowhere

Then model.MSO() ends in SIGSEGV/SIGTRAP inside dmpermlib.MSO. The crash
is memory corruption, so without help it shows up only now and then, often
first after repeated calls. With the macOS malloc debug variables
MallocScribble=1 MallocGuardEdges=1 it happens on the first call.

Reference cases that do not crash:
    [[1]]               no overdetermined part, no zero column -> []
    [[1, 0], [1, 0]]    zero column, overdetermined part       -> [[1, 0]]
    [[1], [0]]          zero row                                -> [[1]]

The only check for unused unknowns is in Lint(), which prints
"Variable x2 is not included in model" as a warning and changes nothing.

Probable cause (not verified with a rebuilt extension): MSOAlg::MSO calls
SM.Plus(), which keeps the whole matrix when the redundancy of the last
Dulmage-Mendelsohn block is <= 0. FindMSO then recurses on a matrix with
more columns than rows (fi = rows - cols < 0, never 1). When there is an
overdetermined part, Plus() keeps only M+, and the zero column is removed.
Possible fixes: return no MSO sets when M+ is empty (in MSOAlg::MSO or in
the Python wrapper), and/or remove zero columns before calling dmpermlib.

Run (each case in its own process so the crash does not stop the script):
    PYTHONPATH=src python src/lab/mso_zero_column.py
"""
import json
import os
import subprocess
import sys

CASES = [
    ([[1, 0]], "zero column, no overdetermined part"),
    ([[0, 1]], "zero column first, no overdetermined part"),
    ([[1]], "no zero column, no overdetermined part"),
    ([[1, 0], [1, 0]], "zero column, overdetermined part"),
    ([[1], [0]], "zero row"),
]


def child(X):
    import faultdiagnosistoolbox as fdt

    model = fdt.DiagnosisModel({"type": "MatrixStruc", "X": X, "F": [], "Z": []})
    print([list(map(int, m)) for m in model.MSO()], flush=True)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        child(json.loads(sys.argv[1]))
        sys.exit(0)
    env = dict(os.environ, MallocScribble="1", MallocGuardEdges="1")
    for X, what in CASES:
        r = subprocess.run([sys.executable, __file__, json.dumps(X)], env=env, capture_output=True, text=True)
        out = r.stdout.strip() if r.returncode == 0 else f"CRASH (exit code {r.returncode})"
        print(f"{str(X):18s} {what:42s} MSO() -> {out}")
