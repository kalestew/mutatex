#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""generate_mutatex_position_list_from_mutinfo.py

Create a MutateX-style position list (e.g. "AA25") directly from a *mutinfo.txt*
file produced by Rosetta Flex ddG saturation mutagenesis pipelines.

Each line of *mutinfo.txt* is expected to follow the standard four-column format
used by RosettaDDG scripts, for example::

    A.A.25.A,A-A25A,A25A,A25
    A.A.25.C,A-A25C,A25C,A25

The first column encodes the mutation in the form
"<chain>.<wildtype AA 1-letter>.<resnum>.<mutant AA 1-letter>".  This script
extracts the chain, wild-type amino-acid, and residue number to build a position
identifier accepted by MutateX ("<WT AA><chain><resnum>").

Usage
-----
$ python generate_mutatex_position_list_from_mutinfo.py \
        -m flexddg/mutinfo.txt \
        -o position_list.txt

The output file will contain unique, sorted position identifiers – one per line.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List, Set, Tuple


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate MutateX position list from mutinfo.txt")
    p.add_argument(
        "-m",
        "--mutinfo",
        default="flexddg/mutinfo.txt",
        help="Path to mutinfo.txt (default: flexddg/mutinfo.txt)",
    )
    p.add_argument(
        "-o",
        "--output",
        default="position_list.txt",
        help="Filename for the generated position list (default: position_list.txt)",
    )
    return p.parse_args()


def _extract_position(token: str):
    """Return (wt_aa, chain, resnum_str, resnum_int) parsed from the mutinfo first-field token.

    The expected token format is '<chain>.<wt_aa>.<resnum>.<mut_aa>'.  Minimal
    validation is performed – the function returns *None* for malformed tokens.
    """
    parts = token.split(".")
    if len(parts) < 3:
        return None  # type: ignore[return-value]

    chain, wt_aa, resnum_str = parts[0], parts[1], parts[2]

    # Validate wild-type amino-acid single letter code
    if not re.fullmatch(r"[A-Z]", wt_aa):
        return None  # type: ignore[return-value]

    # Residue numbers may sometimes include insertion codes (e.g., '25A').  We
    # split any trailing letters to obtain the numeric portion for sorting, but
    # keep the original string for later reconstruction.
    match = re.fullmatch(r"(\d+)([A-Z]?)", resnum_str)
    if match is None:
        return None  # type: ignore[return-value]

    resnum_int = int(match.group(1))
    return wt_aa, chain, resnum_str, resnum_int


def _position_identifier(wt_aa: str, chain: str, resnum_str: str) -> str:
    """Return MutateX position identifier string (e.g., 'AA25')."""
    return f"{wt_aa}{chain}{resnum_str}"


def generate_position_list(mutinfo_path: Path) -> List[str]:
    """Parse *mutinfo.txt* and return a sorted list of unique position IDs."""
    positions: Set[Tuple[str, str, str, int]] = set()

    with mutinfo_path.open() as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            first_field = line.split(",", 1)[0]
            parsed = _extract_position(first_field)
            if parsed is None:
                print(f"[WARN] Skipping unrecognised line: {line}")
                continue
            positions.add(parsed)

    # Sort by chain then residue number for readability
    sorted_positions = sorted(positions, key=lambda t: (t[1], t[3]))
    return [_position_identifier(t[0], t[1], t[2]) for t in sorted_positions]


def main() -> None:  # noqa: D401
    args = _parse_args()
    mutinfo_path = Path(args.mutinfo).expanduser().resolve()

    if not mutinfo_path.is_file():
        raise SystemExit(f"❌  mutinfo file not found: {mutinfo_path}")

    position_list = generate_position_list(mutinfo_path)

    output_path = Path(args.output).expanduser().resolve()
    with output_path.open("w") as out_f:
        out_f.write("\n".join(position_list) + "\n")

    print(f"✅  Wrote {len(position_list)} position(s) to {output_path}")


if __name__ == "__main__":
    main() 