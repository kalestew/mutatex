#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from Bio import PDB
from Bio.Data.IUPACData import protein_letters_3to1

def parse_residue_spans(span_list):
    """Parses residue span strings like ['A:30-37', 'B:50-55']."""
    residues = []
    for span in span_list:
        chain, rng = span.split(":")
        start, end = map(int, rng.split("-"))
        residues.append((chain, start, end))
    return residues

def generate_position_list(structure, spans, output_file):
    """Generates a list of residues in the format: [1-letter WT residue][Chain ID][Residue Number]."""
    model = structure[0]
    pos_entries = []

    for chain_id, start, end in spans:
        if chain_id not in model:
            print(f"Warning: Chain {chain_id} not found in structure.")
            continue

        chain = model[chain_id]
        for res in chain.get_residues():
            if res.id[0] != ' ':
                continue  # Skip heteroatoms/waters
            resseq = res.id[1]
            if start <= resseq <= end:
                res_3letter = res.resname.capitalize()
                if res_3letter not in protein_letters_3to1:
                    print(f"Skipping unknown residue {res.resname} at {chain_id}{resseq}")
                    continue
                one_letter = protein_letters_3to1[res_3letter]
                pos_entry = f"{one_letter}{chain_id}{resseq}"
                pos_entries.append(pos_entry)

    with open(output_file, 'w') as f:
        for entry in pos_entries:
            f.write(entry + "\n")
    print(f"Generated position list: {output_file} with {len(pos_entries)} entries.")

def main():
    parser = argparse.ArgumentParser(description="Generate a position list for mutagenesis protocol.")
    parser.add_argument("-p", "--pdb", required=True, help="Input PDB file")
    parser.add_argument("-s", "--spans", nargs='+', required=True,
                        help="Residue spans to include, e.g., A:30-37 B:50-60")
    parser.add_argument("-o", "--output", default="position_list.txt",
                        help="Output file name for position list")
    args = parser.parse_args()

    pdb_parser = PDB.PDBParser(QUIET=True)
    structure = pdb_parser.get_structure("input", args.pdb)
    spans = parse_residue_spans(args.spans)

    generate_position_list(structure, spans, args.output)

if __name__ == "__main__":
    main()
