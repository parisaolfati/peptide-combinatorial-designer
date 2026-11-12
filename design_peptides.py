# -*- coding: utf-8 -*-
import argparse
from pathlib import Path
import pandas as pd
import itertools


def generate_peptides(input_file: Path, output_file: Path, wild_peptide: str):
    """Load best mutations, generate single and combination peptide mutants, remove duplicates, and export."""
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    # Load dataset
    if input_file.suffix.lower() == ".xlsx":
        df = pd.read_excel(input_file)
    elif input_file.suffix.lower() in [".csv", ".txt"]:
        df = pd.read_csv(input_file)
    else:
        raise ValueError(f"Unsupported file format: {input_file.suffix}")

    # Validate required columns
    required_cols = {"Position", "Mutation"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in dataset: {missing}")

    # Extract mutated amino acid code from Mutation column
    def extract_mutant_aa(mutation_name):
        return str(mutation_name).split("_")[0].upper()

    df["Mutant_AA"] = df["Mutation"].apply(extract_mutant_aa)

    # List of mutations: [(Position, AminoAcid)]
    mutations = list(zip(df["Position"], df["Mutant_AA"]))

    def mutate(seq, position, new_aa):
        """Substitute amino acid at specified 1-based position."""
        aa1_letter = {
            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
            "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
            "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
        }
        seq_list = list(seq)
        idx = int(position) - 1
        if 0 <= idx < len(seq_list):
            seq_list[idx] = aa1_letter.get(new_aa, "X")
        return "".join(seq_list)

    # Generate single point mutants
    single_mutants = []
    for pos, aa in mutations:
        mutated_seq = mutate(wild_peptide, pos, aa)
        single_mutants.append({
            "Type": "Single",
            "Mutations": f"{aa}{pos}",
            "Sequence": mutated_seq
        })

    # Generate combination mutants (2 and 3 point mutations)
    combination_mutants = []
    for r in [2, 3]:
        for combo in itertools.combinations(mutations, r):
            seq = wild_peptide
            mut_desc = []
            for pos, aa in combo:
                seq = mutate(seq, pos, aa)
                mut_desc.append(f"{aa}{pos}")
            combination_mutants.append({
                "Type": f"{r}-Combo",
                "Mutations": "+".join(mut_desc),
                "Sequence": seq
            })

    # Combine results and remove duplicate sequences
    all_mutants = single_mutants + combination_mutants
    df_out = pd.DataFrame(all_mutants)
    df_out = df_out.drop_duplicates(subset=["Sequence"]).reset_index(drop=True)

    # Save output
    if output_file.suffix.lower() == ".xlsx":
        df_out.to_excel(output_file, index=False, engine="openpyxl")
    else:
        df_out.to_csv(output_file, index=False)

    print(f"File created successfully: {output_file.resolve()}\n")
    print(df_out.head(10).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate combinatorial peptide mutants from screening mutation results.")
    parser.add_argument(
        "input_file",
        nargs="?",
        default="best_mutations_from_all.xlsx",
        type=Path,
        help="Path to the input best mutations file (default: best_mutations_from_all.xlsx)."
    )
    parser.add_argument(
        "-o", "--output",
        default="designed_peptides.xlsx",
        type=Path,
        help="Path for the output designed peptides file (default: designed_peptides.xlsx)."
    )
    parser.add_argument(
        "--wild",
        default="LIHADPPGVGLTGF",
        help="Wild-type peptide sequence (default: LIHADPPGVGLTGF)."
    )
    args = parser.parse_args()
    
    generate_peptides(args.input_file, args.output, args.wild)
