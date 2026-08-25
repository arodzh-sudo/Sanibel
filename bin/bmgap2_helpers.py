#!/usr/bin/env python3
"""
bmgap2_helpers.py - shared boilerplate for the run_bmgap2_*.py host scripts.

The three BMGAP2 steps (AMR, LocusExtractor, BMScan) read the MLST scheme to decide
whether to run, skipping (exit 0) any sample that is not Neisseria meningitidis or
H. influenzae.
"""

import subprocess
import sys

MENINGITIS_SCHEMES = ('neisseria', 'hinfluenzae')
ON_FAIL_MODES = ('fail', 'continue')


def run_tool(cmd, tool_tag, cwd=None, on_fail='fail'):
    if on_fail not in ON_FAIL_MODES:
        raise ValueError(
            f"{tool_tag}: unknown on_fail {on_fail!r}, expected one of {ON_FAIL_MODES}")

    try:
        result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, check=False)
    except Exception as e:
        print(f"{tool_tag}: Unexpected error - {e}", file=sys.stderr)
        sys.exit(1)

    # Replay before the exit branches or a failing subtool's log is discarded
    if result.stdout:
        print(result.stdout)

    if result.returncode != 0:
        print(f"{tool_tag}: subtool exited with code {result.returncode}", file=sys.stderr)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        if on_fail == 'fail':
            sys.exit(1)
        print(f"{tool_tag}: continuing after failure", file=sys.stderr)

    return result


def read_scheme(mlst_file, sample_name, tool_tag):
    try:
        with open(mlst_file) as f:
            fields = f.readline().strip().split()
    except Exception as e:
        print(f"{tool_tag}: Error reading MLST file - {e}", file=sys.stderr)
        sys.exit(1)

    if len(fields) < 2:
        print(f"{tool_tag}: Error - {mlst_file} carries no scheme field", file=sys.stderr)
        sys.exit(1)

    scheme = fields[1]
    if scheme not in MENINGITIS_SCHEMES:
        print(f"{tool_tag}: Skipping {sample_name} "
              f"- not a meningitis species (scheme={scheme})")
        sys.exit(0)

    return scheme
