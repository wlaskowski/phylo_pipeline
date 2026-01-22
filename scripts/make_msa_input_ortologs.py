"""
Snakemake helper script which creates per-cluster FASTA files for 1:1 orthologues.

Script:
    - reads MMseqs2 clustering results
    - filters clusters to keep only genomes listed in 'species.txt'
    - selects clusters that contain exactly one gene per genome
    - writes one FASTA per cluster to 'snakemake.output.msa_dir'

Script is executed by Snakemake.
"""

from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from msa_utils import find_all_clusters, load_sequences


def find_1_1_orthologues(clusters, num_all_genomes):
    n = 0
    orthologues_1_1 = {}
    for cluster_rep in sorted(clusters):
        genomes_in_cluster = clusters[cluster_rep]

        genomes_all_unique = []
        for genome_and_gene in genomes_in_cluster:
            genome = genome_and_gene.split("|")[0]
            genomes_all_unique.append(genome)

        if len(set(genomes_all_unique)) == num_all_genomes and len(genomes_all_unique) == num_all_genomes:
            cluster_id = f"cluster{n}"
            orthologues_1_1[cluster_id] = genomes_in_cluster
            n += 1

    return orthologues_1_1


def save_fasta_for_msa(seqs_all, orthologues_1_1, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for cluster_id, members in orthologues_1_1.items():
        outpath = outdir/f"{cluster_id}.faa"
        records = []

        for gene_prot_id in members:
            genome = gene_prot_id.split("|")[0]
            seq = seqs_all[gene_prot_id]
            rec = SeqRecord(
                Seq(seq),
                id=genome,
                description=""
            )
            records.append(rec)
    
        SeqIO.write(records, outpath, "fasta")


all_clusters_file = snakemake.input.clusters
species_file      = snakemake.input.species
proteomes_dir     = snakemake.input.proteomes_dir
msa_outdir        = snakemake.output.msa_dir

ACCESSIONS = {a.strip() for a in open(species_file) if a.strip()}
num_all_genomes = len(ACCESSIONS)


all_clusters = find_all_clusters(
                infile=all_clusters_file, 
                ACCESSIONS=ACCESSIONS
                )
orthologues_1_1 = find_1_1_orthologues(
                    clusters=all_clusters,
                    num_all_genomes=num_all_genomes
                    )
seq_dict = load_sequences(
            data_dir=proteomes_dir
            )
save_fasta_for_msa(
    seqs_all=seq_dict, 
    orthologues_1_1=orthologues_1_1,
    outdir=msa_outdir
    )







    


