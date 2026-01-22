"""
Snakemake helper script that creates per-cluster FASTA files for paralogue analysis stage.

Script:
    - reads MMseqs2 clustering results using functions from msa_utils module
    - filters by genomes listed in 'species.txt'
    - writes one FASTA per cluster containing all copies (including paralogs)
    writes output to 'snakemake.output.msa_dir.'

Thsi script is excecuted by Snakemake.
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path

from msa_utils import find_all_clusters, load_sequences


def save_fasta_for_msa_all_copies(seqs_all, clusters, outdir, min_size=2):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    n = 0
    for cluster_rep in sorted(clusters):
        members = clusters[cluster_rep]
        if len(members) < min_size:
            continue

        outpath = outdir / f"cluster_paralogues_{n}.faa"
        n += 1

        records = []
        for gene_prot_id in members:
            seq = seqs_all.get(gene_prot_id)
            if seq is None:
                continue

            unique_id = gene_prot_id.replace("|", "__")

            records.append(SeqRecord(Seq(seq), id=unique_id, description=""))
        SeqIO.write(records, outpath, "fasta")


all_clusters_file = snakemake.input.clusters
species_file = snakemake.input.species
proteomes_dir = snakemake.input.proteomes_dir
msa_outdir = snakemake.output.msa_dir

ACCESSIONS = {a.strip().split()[0] for a in open(species_file) if a.strip()}

clusters = find_all_clusters(
    infile=all_clusters_file, 
    ACCESSIONS=ACCESSIONS
    )

seqs = load_sequences(proteomes_dir)
save_fasta_for_msa_all_copies(seqs, clusters, msa_outdir, min_size=2)



