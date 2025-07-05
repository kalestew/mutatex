#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ddg2heatmap_wrapper.py – helper around MutateX `ddg2heatmap`

This wrapper reads a MutateX *position list* file, removes any residue
identifiers for which the corresponding energy file is missing in the data
directory, then invokes `ddg2heatmap` with the filtered list.

Motivation
----------
When running `rosetta_ddg_aggregate --mutatex-convert`, some residues may lack
energy files (e.g. because the original Flex ddG replicas failed). Passing those
positions to `ddg2heatmap` causes it to abort with an "Couldn't open energy
file ..." error.  This wrapper automatically excludes the missing entries so
you still get a heat-map for the residues that *are* present.

Usage example
-------------
    python ddg2heatmap_wrapper.py \
        -p 41D1_rosettaPreped.pdb \
        -d mutatex_compatible \
        -l residues.txt \
        -q position_list.txt \
        -o heatmap.pdf

All unknown options after the recognised ones are forwarded verbatim to
`ddg2heatmap`. Use this to pass flags such as `-c` (cmap) or `-t` (transpose).
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Tuple


def _parse_args() -> Tuple[argparse.Namespace, List[str]]:
    p = argparse.ArgumentParser(add_help=False)  # we supply -h via ddg2heatmap

    required = p.add_argument_group("required arguments")
    required.add_argument("-p", "--pdb", required=True, help="Input PDB file")
    required.add_argument("-d", "--data-directory", required=True, help="Directory containing MutateX energy files")
    required.add_argument("-l", "--mutation-list", required=True, help="Mutation list file (single-letter residues)")
    required.add_argument("-q", "--position-list", required=True, help="Position list file (e.g. CA23)")

    optional = p.add_argument_group("optional arguments")
    optional.add_argument("-o", "--output", default="heatmap.pdf", help="Output filename passed to ddg2heatmap (-o)")
    optional.add_argument("--keep-temp", action="store_true", help="Keep the filtered position list file instead of deleting it")
    optional.add_argument("-h", "--help", action="help", help="Show this help and exit")

    # allow every other arg; they will be forwarded to ddg2heatmap
    return p.parse_known_args()


def _filter_positions(position_file: Path, data_dir: Path) -> Tuple[Path, int]:
    """Return path to filtered position list (tmp file) and number of removed entries."""
    positions: List[str] = [line.strip() for line in position_file.read_text().splitlines() if line.strip()]
    available: List[str] = []
    missing: List[str] = []

    for pos in positions:
        energy_file = data_dir / pos
        if energy_file.is_file():
            available.append(pos)
        else:
            missing.append(pos)

    if missing:
        sys.stderr.write(f"[INFO] Excluding {len(missing)} residue(s) with missing energy files:\n")
        for m in missing:
            sys.stderr.write(f"  • {m}\n")
    else:
        sys.stderr.write("[INFO] All positions have energy files – no filtering needed.\n")

    # write available positions to tmp file
    tmp_handle = tempfile.NamedTemporaryFile("w", prefix="filtered_pos_", suffix=".txt", delete=False)
    for pos in available:
        tmp_handle.write(pos + "\n")
    tmp_handle.flush()
    tmp_handle.close()
    return Path(tmp_handle.name), len(missing)


def main() -> None:  # noqa: D401
    args, extra = _parse_args()

    pdb_path = Path(args.pdb).expanduser().resolve()
    data_dir = Path(args.data_directory).expanduser().resolve()
    mutlist_path = Path(args.mutation_list).expanduser().resolve()
    poslist_path = Path(args.position_list).expanduser().resolve()

    for path in (pdb_path, mutlist_path, poslist_path):
        if not path.is_file():
            sys.exit(f"❌  File not found: {path}")

    if not data_dir.is_dir():
        sys.exit(f"❌  Data directory not found: {data_dir}")

    filtered_pos_path, n_removed = _filter_positions(poslist_path, data_dir)

    ddg_cmd = [
        "ddg2heatmap",
        "-p",
        str(pdb_path),
        "-d",
        str(data_dir),
        "-l",
        str(mutlist_path),
        "-q",
        str(filtered_pos_path),
        "-o",
        args.output,
    ] + extra

    sys.stderr.write("[INFO] Running: " + " ".join(ddg_cmd) + "\n")

    try:
        subprocess.run(ddg_cmd, check=True)
    finally:
        if not args.keep_temp:
            try:
                os.unlink(filtered_pos_path)
            except OSError:
                pass

    sys.stderr.write("[INFO] ddg2heatmap completed.\n")
    if n_removed:
        sys.stderr.write(f"[INFO] Filtered out {n_removed} missing position(s).\n")


if __name__ == "__main__":
    main() 