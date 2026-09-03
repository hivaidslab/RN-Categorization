# RN-Categorization

Categorizes HIV-1 sequencing reads (FASTA) by how many R-site(RBEIII) and N-site(NF-KB) motif occurrences they contain, and writes per-file counts
of each category.

## Background

Each read is scanned for the canonical **R-site** and **N-site** motifs, plus
every single-base-substitution variant of each. A read is then binned by how many R-site and
N-site hits it has:

| R-sites | N-sites | Category |
|---------|---------|----------|
| 1       | 2       | RN2      |
| 1       | 3       | RN3      |
| 1       | 4       | RN4      |
| 2       | 2       | R2N2     |
| 2       | 3       | R2N3     |
| 2       | 4       | R2N4     |
| any other combination | | others |

## Requirements

- Python 3.8+
- [Biopython](https://biopython.org/)

Install with:

```bash
pip install -r requirements.txt
```

## Usage

```bash
python RN_Categorization.py --input-dir path/to/fasta_folder --output results/RN_Counts.txt
```

- `--input-dir` — a folder containing one or more `.fasta` files, or a glob
  pattern to a specific set of files (e.g. `"data/*.fasta"`).
- `--output` — path to the text file the per-file category counts will be
  written to.

### Example

```bash
python RN_Categorization.py --input-dir example_data/ --output results/RN_Counts.txt
```

Output (`RN_Counts.txt`) looks like:

```
Counts for file: example_data/sample1.fasta
RN2: 12
RN3: 4
RN4: 1
R2N2: 0
R2N3: 0
R2N4: 0
others: 3
```

## Input format

Standard FASTA files (`.fasta`), one record per read.

## Reproducing on your own data

1. Clone this repo and install requirements.
2. Place your `.fasta` files in a folder.
3. Run the command above, pointing `--input-dir` at that folder.

No paths need to be edited in the script itself — everything is passed via
command-line arguments.
