#!/usr/bin/env python3
"""
RN_Categorization.py

Categorize HIV-1 sequenced FASTA reads by the number of
R-site(RBEIII) and N-site(NF-KB) motif occurrences they contain (e.g. RN2, RN3, RN4,
R2N2, R2N3, R2N4, or "others").

Each input FASTA record's sequence is scanned for all 26 single-base
variants of the canonical R-site and all 300+ single-base variants of
the canonical N-site. The read is then binned into a category based on how many R-site and N-site
matches it contains.

Usage
-----
    python RN_Categorization.py --input-dir data/fasta --output results/RN_Counts.txt

    # Process a single file instead of a directory
    python RN_Categorization.py --input-dir data/sample1.fasta --output results/RN_Counts.txt

Requirements
------------
    biopython (Bio.SeqIO)

See README.md for details on the motif definitions and category logic.
"""

import argparse
import glob
import os
from Bio import SeqIO

# Canonical R-site and N-site sequences, plus every single-base-substitution
# variant of each (accounts for one sequencing/PCR error per site).
R_SITES = [
    'ACTGCTGA', 'ATCGCTGA', 'TCTGCTGA', 'GCTGCTGA', 'CCTGCTGA', 'AATGCTGA',
    'ATTGCTGA', 'AGTGCTGA', 'ACAGCTGA', 'ACGGCTGA', 'ACCGCTGA', 'ACTACTGA',
    'ACTTCTGA', 'ACTCCTGA', 'ACTGATGA', 'ACTGTTGA', 'ACTGGTGA', 'ACTGCAGA',
    'ACTGCGGA', 'ACTGCCGA', 'ACTGCTAA', 'ACTGCTTA', 'ACTGCTCA', 'ACTGCTGT',
    'ACTGCTGG', 'ACTGCTGC',
]

N_SITES = [
    'GGGACTAAG', 'AGGACTTTCC', 'TGGACTTTCC', 'GGGACTTTCC', 'CGGACTTTCC',
    'AAGACTTTCC', 'ATGACTTTCC', 'ACGACTTTCC', 'AGAACTTTCC', 'AGTACTTTCC',
    'AGCACTTTCC', 'AGGTCTTTCC', 'AGGGCTTTCC', 'AGGCCTTTCC', 'AGGAATTTCC',
    'AGGATTTTCC', 'AGGAGTTTCC', 'AGGACATTCC', 'AGGACGTTCC', 'AGGACCTTCC',
    'AGGACTATCC', 'AGGACTGTCC', 'AGGACTCTCC', 'AGGACTTACC', 'AGGACTTGCC',
    'AGGACTTCCC', 'AGGACTTTAC', 'AGGACTTTTC', 'AGGACTTTGC', 'AGGACTTTCA',
    'AGGACTTTCT', 'AGGACTTTCG', 'TAGACTTTCC', 'TTGACTTTCC', 'TCGACTTTCC',
    'TGAACTTTCC', 'TGTACTTTCC', 'TGCACTTTCC', 'TGGTCTTTCC', 'TGGGCTTTCC',
    'TGGCCTTTCC', 'TGGAATTTCC', 'TGGATTTTCC', 'TGGAGTTTCC', 'TGGACATTCC',
    'TGGACGTTCC', 'TGGACCTTCC', 'TGGACTATCC', 'TGGACTGTCC', 'TGGACTCTCC',
    'TGGACTTACC', 'TGGACTTGCC', 'TGGACTTCCC', 'TGGACTTTAC', 'TGGACTTTTC',
    'TGGACTTTGC', 'TGGACTTTCA', 'TGGACTTTCT', 'TGGACTTTCG', 'GAGACTTTCC',
    'GTGACTTTCC', 'GCGACTTTCC', 'GGAACTTTCC', 'GGTACTTTCC', 'GGCACTTTCC',
    'GGGTCTTTCC', 'GGGGCTTTCC', 'GGGCCTTTCC', 'GGGAATTTCC', 'GGGATTTTCC',
    'GGGAGTTTCC', 'GGGACATTCC', 'GGGACGTTCC', 'GGGACCTTCC', 'GGGACTATCC',
    'GGGACTGTCC', 'GGGACTCTCC', 'GGGACTTACC', 'GGGACTTGCC', 'GGGACTTCCC',
    'GGGACTTTAC', 'GGGACTTTTC', 'GGGACTTTGC', 'GGGACTTTCA', 'GGGACTTTCT',
    'GGGACTTTCG', 'CAGACTTTCC', 'CTGACTTTCC', 'CCGACTTTCC', 'CGAACTTTCC',
    'CGTACTTTCC', 'CGCACTTTCC', 'CGGTCTTTCC', 'CGGGCTTTCC', 'CGGCCTTTCC',
    'CGGAATTTCC', 'CGGATTTTCC', 'CGGAGTTTCC', 'CGGACATTCC', 'CGGACGTTCC',
    'CGGACCTTCC', 'CGGACTATCC', 'CGGACTGTCC', 'CGGACTCTCC', 'CGGACTTACC',
    'CGGACTTGCC', 'CGGACTTCCC', 'CGGACTTTAC', 'CGGACTTTTC', 'CGGACTTTGC',
    'CGGACTTTCA', 'CGGACTTTCT', 'CGGACTTTCG', 'GAAACTTTCC', 'GATACTTTCC',
    'GACACTTTCC', 'GAGTCTTTCC', 'GAGGCTTTCC', 'GAGCCTTTCC', 'GAGAATTTCC',
    'GAGATTTTCC', 'GAGAGTTTCC', 'GAGACATTCC', 'GAGACGTTCC', 'GAGACCTTCC',
    'GAGACTATCC', 'GAGACTGTCC', 'GAGACTCTCC', 'GAGACTTACC', 'GAGACTTGCC',
    'GAGACTTCCC', 'GAGACTTTAC', 'GAGACTTTTC', 'GAGACTTTGC', 'GAGACTTTCA',
    'GAGACTTTCT', 'GAGACTTTCG', 'GTAACTTTCC', 'GTTACTTTCC', 'GTCACTTTCC',
    'GTGTCTTTCC', 'GTGGCTTTCC', 'GTGCCTTTCC', 'GTGAATTTCC', 'GTGATTTTCC',
    'GTGAGTTTCC', 'GTGACATTCC', 'GTGACGTTCC', 'GTGACCTTCC', 'GTGACTATCC',
    'GTGACTGTCC', 'GTGACTCTCC', 'GTGACTTACC', 'GTGACTTGCC', 'GTGACTTCCC',
    'GTGACTTTAC', 'GTGACTTTTC', 'GTGACTTTGC', 'GTGACTTTCA', 'GTGACTTTCT',
    'GTGACTTTCG', 'GCAACTTTCC', 'GCTACTTTCC', 'GCCACTTTCC', 'GCGTCTTTCC',
    'GCGGCTTTCC', 'GCGCCTTTCC', 'GCGAATTTCC', 'GCGATTTTCC', 'GCGAGTTTCC',
    'GCGACATTCC', 'GCGACGTTCC', 'GCGACCTTCC', 'GCGACTATCC', 'GCGACTGTCC',
    'GCGACTCTCC', 'GCGACTTACC', 'GCGACTTGCC', 'GCGACTTCCC', 'GCGACTTTAC',
    'GCGACTTTTC', 'GCGACTTTGC', 'GCGACTTTCA', 'GCGACTTTCT', 'GCGACTTTCG',
    'GGATCTTTCC', 'GGAGCTTTCC', 'GGACCTTTCC', 'GGAAATTTCC', 'GGAATTTTCC',
    'GGAAGTTTCC', 'GGAACATTCC', 'GGAACGTTCC', 'GGAACCTTCC', 'GGAACTATCC',
    'GGAACTGTCC', 'GGAACTCTCC', 'GGAACTTACC', 'GGAACTTGCC', 'GGAACTTCCC',
    'GGAACTTTAC', 'GGAACTTTTC', 'GGAACTTTGC', 'GGAACTTTCA', 'GGAACTTTCT',
    'GGAACTTTCG', 'GGTTCTTTCC', 'GGTGCTTTCC', 'GGTCCTTTCC', 'GGTAATTTCC',
    'GGTATTTTCC', 'GGTAGTTTCC', 'GGTACATTCC', 'GGTACGTTCC', 'GGTACCTTCC',
    'GGTACTATCC', 'GGTACTGTCC', 'GGTACTCTCC', 'GGTACTTACC', 'GGTACTTGCC',
    'GGTACTTCCC', 'GGTACTTTAC', 'GGTACTTTTC', 'GGTACTTTGC', 'GGTACTTTCA',
    'GGTACTTTCT', 'GGTACTTTCG', 'GGCTCTTTCC', 'GGCGCTTTCC', 'GGCCCTTTCC',
    'GGCAATTTCC', 'GGCATTTTCC', 'GGCAGTTTCC', 'GGCACATTCC', 'GGCACGTTCC',
    'GGCACCTTCC', 'GGCACTATCC', 'GGCACTGTCC', 'GGCACTCTCC', 'GGCACTTACC',
    'GGCACTTGCC', 'GGCACTTCCC', 'GGCACTTTAC', 'GGCACTTTTC', 'GGCACTTTGC',
    'GGCACTTTCA', 'GGCACTTTCT', 'GGCACTTTCG', 'GGGTATTTCC', 'GGGTTTTTCC',
    'GGGTGTTTCC', 'GGGTCATTCC', 'GGGTCGTTCC', 'GGGTCCTTCC', 'GGGTCTATCC',
    'GGGTCTGTCC', 'GGGTCTCTCC', 'GGGTCTTACC', 'GGGTCTTGCC', 'GGGTCTTCCC',
    'GGGTCTTTAC', 'GGGTCTTTTC', 'GGGTCTTTGC', 'GGGTCTTTCA', 'GGGTCTTTCT',
    'GGGTCTTTCG', 'GGGGATTTCC', 'GGGGTTTTCC', 'GGGGGTTTCC', 'GGGGCATTCC',
    'GGGGCCTTCC', 'GGGGCTATCC', 'GGGGCTGTCC', 'GGGGCTCTCC', 'GGGGCTTACC',
    'GGGGCTTGCC', 'GGGGCTTCCC', 'GGGGCTTTAC', 'GGGGCTTTTC', 'GGGGCTTTGC',
    'GGGGCTTTCA', 'GGGGCTTTCT', 'GGGGCTTTCG', 'GGGCATTTCC', 'GGGCTTTTCC',
    'GGGCGTTTCC', 'GGGCCATTCC', 'GGGCCGTTCC', 'GGGCCCTTCC', 'GGGCCTATCC',
    'GGGCCTGTCC', 'GGGCCTCTCC', 'GGGCCTTACC', 'GGGCCTTGCC', 'GGGCCTTCCC',
    'GGGCCTTTAC', 'GGGCCTTTTC', 'GGGCCTTTGC', 'GGGCCTTTCA', 'GGGCCTTTCT',
    'GGGCCTTTCG', 'GGGAAATTCC', 'GGGAAGTTCC', 'GGGAACTTCC', 'GGGAATATCC',
    'GGGAATGTCC', 'GGGAATCTCC', 'GGGAATTACC', 'GGGAATTGCC', 'GGGAATTCCC',
    'GGGAATTTAC', 'GGGAATTTTC', 'GGGAATTTGC', 'GGGAATTTCA', 'GGGAATTTCT',
    'GGGAATTTCG', 'GGGATATTCC', 'GGGATGTTCC', 'GGGATCTTCC', 'GGGATTATCC',
    'GGGATTGTCC', 'GGGATTCTCC', 'GGGATTTACC', 'GGGATTTGCC', 'GGGATTTCCC',
    'GGGATTTTAC', 'GGGATTTTTC', 'GGGATTTTGC', 'GGGATTTTCA', 'GGGATTTTCT',
    'GGGATTTTCG', 'GGGAGATTCC', 'GGGAGGTTCC', 'GGGAGCTTCC', 'GGGAGTATCC',
    'GGGAGTGTCC', 'GGGAGTCTCC', 'GGGAGTTACC', 'GGGAGTTGCC', 'GGGAGTTCCC',
    'GGGAGTTTAC', 'GGGAGTTTTC', 'GGGAGTTTGC', 'GGGAGTTTCA', 'GGGAGTTTCT',
    'GGGAGTTTCG', 'GGGACAATCC', 'GGGACAGTCC', 'GGGACACTCC', 'GGGACATACC',
    'GGGACATGCC', 'GGGACATCCC', 'GGGACATTAC', 'GGGACATTTC', 'GGGACATTGC',
    'GGGACATTCA', 'GGGACATTCT', 'GGGACATTCG', 'GGGACGATCC', 'GGGACGGTCC',
    'GGGACGCTCC', 'GGGACGTACC', 'GGGACGTGCC', 'GGGACGTCCC', 'GGGACGTTAC',
    'GGGACGTTTC', 'GGGACGTTGC', 'GGGACGTTCA', 'GGGACGTTCT', 'GGGACGTTCG',
    'GGGACCATCC', 'GGGACCGTCC', 'GGGACCCTCC', 'GGGACCTACC', 'GGGACCTGCC',
    'GGGACCTCCC', 'GGGACCTTAC', 'GGGACCTTTC', 'GGGACCTTGC', 'GGGACCTTCA',
    'GGGACCTTCT', 'GGGACCTTCG', 'GGGACTAACC', 'GGGACTAGCC', 'GGGACTACCC',
    'GGGACTATAC', 'GGGACTATTC', 'GGGACTATGC', 'GGGACTATCA', 'GGGACTATCT',
    'GGGACTATCG', 'GGGACTGACC', 'GGGACTGGCC', 'GGGACTGCCC', 'GGGACTGTAC',
    'GGGACTGTTC', 'GGGACTGTGC', 'GGGACTGTCA', 'GGGACTGTCT', 'GGGACTGTCG',
    'GGGACTCACC', 'GGGACTCGCC', 'GGGACTCCCC', 'GGGACTCTAC', 'GGGACTCTTC',
    'GGGACTCTGC', 'GGGACTCTCA', 'GGGACTCTCT', 'GGGACTCTCG', 'GGGACTTAAC',
    'GGGACTTATC', 'GGGACTTAGC', 'GGGACTTACA', 'GGGACTTACT', 'GGGACTTACG',
    'GGGACTTGAC', 'GGGACTTGTC', 'GGGACTTGGC', 'GGGACTTGCA', 'GGGACTTGCT',
    'GGGACTTGCG', 'GGGACTTCAC', 'GGGACTTCTC', 'GGGACTTCGC', 'GGGACTTCCA',
    'GGGACTTCCT', 'GGGACTTCCG', 'GGGACTTTAA', 'GGGACTTTAT', 'GGGACTTTAG',
    'GGGACTTTTA', 'GGGACTTTTT', 'GGGACTTTTG', 'GGGACTTTGA', 'GGGACTTTGT',
    'GGGACTTTGG', 'GGGGCGTTCC', 'AGGGCGTTCC', 'TGGGCGTTCC', 'CGGGCGTTCC',
    'GAGGCGTTCC', 'GTGGCGTTCC', 'GCGGCGTTCC', 'GGAGCGTTCC', 'GGTGCGTTCC',
    'GGCGCGTTCC', 'GGGACGTTCC', 'GGGTCGTTCC', 'GGGCCGTTCC', 'GGGGAGTTCC',
    'GGGGTGTTCC', 'GGGGGGTTCC', 'GGGGCATTCC', 'GGGGCTTTCC', 'GGGGCCTTCC',
    'GGGGCGATCC', 'GGGGCGGTCC', 'GGGGCGCTCC', 'GGGGCGTACC', 'GGGGCGTGCC',
    'GGGGCGTCCC', 'GGGGCGTTAC', 'GGGGCGTTTC', 'GGGGCGTTGC', 'GGGGCGTTCA',
    'GGGGCGTTCT', 'GGGGCGTTCG',
]

# Maps (R-site count, N-site count) pairs to a category.
CATEGORY_MAP = {
    (1, 2): "RN2",
    (1, 3): "RN3",
    (1, 4): "RN4",
    (2, 2): "R2N2",
    (2, 3): "R2N3",
    (2, 4): "R2N4",
}


def categorize_sequence(sequence):
    """Return the category for a single sequence based on R/N site counts."""
    count_r_sites = sum(sequence.count(site) for site in R_SITES)
    count_n_sites = sum(sequence.count(site) for site in N_SITES)
    return CATEGORY_MAP.get((count_r_sites, count_n_sites), "others")


def process_file(input_file):
    """Return a dict of category -> read count for one FASTA file."""
    category_counts = {
        "RN2": 0, "RN3": 0, "RN4": 0,
        "R2N2": 0, "R2N3": 0, "R2N4": 0,
        "others": 0,
    }
    with open(input_file, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            category = categorize_sequence(str(record.seq))
            category_counts[category] += 1
    return category_counts


def main():
    parser = argparse.ArgumentParser(
        description="Categorize FASTA reads by R-site/N-site motif counts."
    )
    parser.add_argument(
        "--input-dir", required=True,
        help="Directory containing input .fasta files, or a path/glob to a single file "
             "(e.g. 'data/*.fasta')."
    )
    parser.add_argument(
        "--output", required=True,
        help="Path to the output counts text file."
    )
    args = parser.parse_args()

    if os.path.isdir(args.input_dir):
        input_files = sorted(glob.glob(os.path.join(args.input_dir, "*.fasta")))
    else:
        input_files = sorted(glob.glob(args.input_dir))

    if not input_files:
        raise SystemExit(f"No .fasta files found for input: {args.input_dir}")

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    with open(args.output, "w") as out_handle:
        for input_file in input_files:
            counts = process_file(input_file)
            out_handle.write(f"Counts for file: {input_file}\n")
            for category, count in counts.items():
                out_handle.write(f"{category}: {count}\n")
            out_handle.write("\n")
            print(f"Processed {input_file}")

    print(f"Done. Results written to {args.output}")


if __name__ == "__main__":
    main()
