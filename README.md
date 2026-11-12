# Peptide Combinatorial Designer

A Python script to generate single and combinatorial (2-point and 3-point) mutant peptide sequences based on screened optimal mutation positions, filter duplicates, and export structured Excel files.

## Features

- Takes optimized mutation positions as input.
- Automatically handles 3-letter to 1-letter amino acid conversion.
- Generates single point mutants as well as 2-combo and 3-combo combinatorial variations of the wild-type peptide.
- Deduplicates sequences and exports an organized report.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

Run with default settings:
```bash
python design_peptides.py
```

Or specify custom inputs, outputs, and wild-type sequences:
```bash
python design_peptides.py best_mutations_from_all.xlsx -o custom_designed.xlsx --wild LIHADPPGVGLTGF
```
