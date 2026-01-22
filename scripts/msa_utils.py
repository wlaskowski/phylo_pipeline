"""
Utilities for preparing pre-cluster FAST files for MSA step.

Module contains two helper functions used by Snakemake rules for:
- parsing MMseqs2 cluster TSV output and filtering by accessions,
- loading protein sequences from a directory of renamed proteomes (FASTA)
"""

from collections import defaultdict
from pathlib import Path
from Bio import SeqIO

def find_all_clusters(infile, ACCESSIONS):
    """
    Parse MMseqs2 cluster TSV and return members for each cluster representative.
    """
    with open(infile) as f:
        all_clusters = defaultdict(list)

        for line in f:
            line = line.rstrip().split()
            cluster_rep = line[0]
            genome_member = line[1].split("|")[0]
            genome_member_gene_id = line[1]
            if genome_member in ACCESSIONS:
                all_clusters[cluster_rep].append(genome_member_gene_id)
    
    return all_clusters


def load_sequences(data_dir):
    """
    Load all protein sequences from a directory containing .faa files.
    """
    data_dir = Path(data_dir)
    seqs = {}

    for faa_file in data_dir.glob("*.faa"):
        for rec in SeqIO.parse(str(faa_file), "fasta"):
            seqs[rec.id] = str(rec.seq)

    return seqs