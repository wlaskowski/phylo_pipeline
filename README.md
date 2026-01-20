# Phylogenomic pipeline

This repository contains a Snakemake-based phylogenomic pipeline developed for species tree reconstruction using both 1:1 orthologues and paralogous gene families.

## Summary

The pipeline:
- downloads proteomes from NCBI Datasets,
- clusters proteins using MMseqs2,
- builds MSAs and gene trees for orthologues and paralogues
- computes species trees using:
	- IQTREE for majority-ruled consensus (case: 1:1 orthologues),
	- ASTRAL for species supertree (case: 1:1 orthologues),
	- ASTRAL-Pro (case: paralogues).

The workflow was developed for a specific comparative genomics project.

## Requirements

The pipeline requires the following software:
- Snakemake
- Conda
- Java (required for ASTRAL / ASTRAL-Pro)  

All other dependencies (Python packages and external tools) are installed automatically via Conda environments defined in `env/pipeline.yaml` and managed by Snakemake.  

Pipeline requires the following to be in the same working directory:
- `Snakefile`,
- `env/pipeline.yaml`,
- `config.yaml`,
- `species.txt`
- `scripts/` directory with `rename_proteomes.sh`, `msa_utils.py`, `make_msa_input_ortologs.py`, `make_msa_input_paralogs.py` 

## Usage

Run:
	```bash
	snakemake -p -j <N> --use-conda   

Options:
- -p - print shell commands executed by Snakemake,
- -j <N> - specify number of parallel jobs according to available CPU cores,
- --use-conda - automatically create and use Conda environments used in the workflow.   

Optionally, if you want to try using the pipeline on a different set of proteomes:
- Edit 'species.txt' to include new genome accessions.


