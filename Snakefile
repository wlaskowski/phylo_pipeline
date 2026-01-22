# Phylogenetic pipeline for inferring majority-rule consensus and supertree from given proteomes

# Overview:
# 1) Download proteomes from NCBI Datasets for RefSeq accessions listed in species.txt file
# 2) Prefix FASTE headers with <ACCESSION|> to make unique IDs across all genomes
# 3) Cluster all proteins with MMseqs2
# 4) Build clusters for two separate workflows:
# 		- 1:1 orthologues:
#			a) IQTREE2 gene trees
#			b) consensus/supertree
# 		- paralogs (all copies):
#			a) FastTree gene trees
#			b) supertree 

import os
from snakemake.io import glob_wildcards

configfile: "config.yaml"
SPECIES_FILE = config["species_file"]

ACCESSIONS = [a.strip() for a in open(SPECIES_FILE) if a.strip()]

rule all:
	input:
		"results/consensus/consensus.contree",
		"results/astral/astral_species_tree.tre",
		"results/astral_pro/astral_pro_species_tree.tre"


rule download_proteomes_zip:
	input:
		SPECIES_FILE
	output:
		"data/proteomes.zip"
	conda: "env/pipeline.yaml"
	shell:
		r"""
		mkdir -p data
		datasets download genome accession \
			--inputfile {input} \
			--include protein \
			--assembly-source refseq \
			--filename ./data/proteomes.zip
		"""


rule unpack_proteomes:
	input:
		"data/proteomes.zip"
	output:
		expand(
			"data/proteomes/ncbi_dataset/data/{acc}/protein.faa",
			acc=ACCESSIONS
		)
	shell:
		r"""
		rm -rf data/proteomes
		mkdir -p data/proteomes
		unzip {input} -d data/proteomes
		"""


rule rename_proteomes:
	input:
		expand(
			"data/proteomes/ncbi_dataset/data/{acc}/protein.faa",
			acc=ACCESSIONS
		)
	output:
		directory("data/proteomes/renamed"),
		expand("data/proteomes/renamed/{acc}.faa", acc=ACCESSIONS)
	shell:
		"bash scripts/rename_proteomes.sh"


rule concat_proteomes:
	input:
		expand(
			"data/proteomes/renamed/{acc}.faa",
			acc=ACCESSIONS
		)
	output:
		"results/proteomes_all.faa"
	shell:
		r"""
		mkdir -p results
		cat {input} > {output}
		"""

rule make_clusters:
	input:
		"results/proteomes_all.faa"
	output:
		fasta = "results/mmseqs/mmseqs_clusters_all_seqs.fasta",
		tsv = "results/mmseqs/mmseqs_clusters_cluster.tsv"
	conda: "env/pipeline.yaml"
	threads: 12
	shell:
		r"""
		mkdir -p results/mmseqs
		mkdir -p results/mmseqs/mmseqs_tmp
		mmseqs easy-cluster \
			{input} \
			results/mmseqs/mmseqs_clusters \
			results/mmseqs/mmseqs_tmp \
			--min-seq-id 0.5 \
			-c 0.8 \
			--cov-mode 1 \
			--threads {threads}
		"""


checkpoint make_msa_input_ortologs:
	input:
		clusters="results/mmseqs/mmseqs_clusters_cluster.tsv",
		species=SPECIES_FILE,
		proteomes_dir="data/proteomes/renamed"
	output:
		msa_dir=directory("results/msa/msa_input_ortologs")
	conda: "env/pipeline.yaml"
	script:
		"scripts/make_msa_input_ortologs.py"


checkpoint make_msa_input_paralogs:
	input:
		clusters="results/mmseqs/mmseqs_clusters_cluster.tsv",
		species=SPECIES_FILE,
		proteomes_dir="data/proteomes/renamed"
	output:
		msa_dir=directory("results/msa/msa_input_paralogs")
	conda: "env/pipeline.yaml"
	script:
		"scripts/make_msa_input_paralogs.py"


def clusters_from_msa_ortologs(wc):
	ck = checkpoints.make_msa_input_ortologs.get(**wc)
	return glob_wildcards(os.path.join(ck.output.msa_dir, "{cluster}.faa")).cluster

def clusters_from_msa_paralogs(wc):
	ck = checkpoints.make_msa_input_paralogs.get(**wc)
	return glob_wildcards(os.path.join(ck.output.msa_dir, "{cluster}.faa")).cluster


rule run_mafft_ortologs:
	input:
		msa_faa = "results/msa/msa_input_ortologs/{cluster}.faa"
	output:
		"results/msa/msa_aligned_ortologs/{cluster}.aln.faa"
	threads: 4
	conda: "env/pipeline.yaml"
	shell:
		r"""
		mkdir -p results/msa/msa_aligned_ortologs
		mafft --auto --thread {threads} {input.msa_faa} > {output}
		"""


rule run_trimal_ortologs:
	input:
		"results/msa/msa_aligned_ortologs/{cluster}.aln.faa"
	output:
		"results/msa/msa_trimmed_ortologs/{cluster}.trimmed.aln.faa"
	conda: "env/pipeline.yaml"
	shell:
		r"""
		mkdir -p results/msa/msa_trimmed_ortologs
		trimal -in {input} -out {output} -automated1 || cp {input} {output}
		"""

rule run_iqtree:
	input:
		aln="results/msa/msa_trimmed_ortologs/{cluster}.trimmed.aln.faa"
	output:
		tree="results/iqtree/{cluster}.treefile"
	threads: 2
	conda: "env/pipeline.yaml"	
	shell:
		r"""
		mkdir -p results/iqtree
		iqtree2 \
			-s {input.aln} \
			-pre results/iqtree/{wildcards.cluster} \
			-m MFP \
			-B 1000 \
			-nt {threads}
		"""


rule run_mafft_paralogs:
	input:
		msa_faa="results/msa/msa_input_paralogs/{cluster}.faa"
	output:
		"results/msa/msa_aligned_paralogs/{cluster}.aln.faa"
	threads: 2
	conda: "env/pipeline.yaml"
	shell:
		r"""
		mkdir -p results/msa/msa_aligned_paralogs
		mafft --auto --thread {threads} {input.msa_faa} > {output}
		"""


rule run_trimal_paralogs:
	input:
		"results/msa/msa_aligned_paralogs/{cluster}.aln.faa"
	output:
		"results/msa/msa_trimmed_paralogs/{cluster}.trimmed.aln.faa"
	conda: "env/pipeline.yaml"
	shell:
		r"""
		mkdir -p results/msa/msa_trimmed_paralogs
		trimal -in {input} -out {output} -automated1 || true
		if [ ! -s "{output}" ]; then
			cp "{input}" "{output}"
		fi
		"""


rule run_fasttree_paralogs:
	input:
		aln="results/msa/msa_trimmed_paralogs/{cluster}.trimmed.aln.faa"
	output:
		tree="results/fasttree_paralogs/{cluster}.tre"
	conda: "env/pipeline.yaml"
	shell:
		r"""
		mkdir -p results/fasttree_paralogs
		FastTree -lg -gamma {input.aln} > {output.tree}
		"""


rule concat_paralog_trees:
	input:
		lambda wc: expand(
			"results/fasttree_paralogs/{cluster}.tre",
			cluster=clusters_from_msa_paralogs(wc)
		)
	output:
		"results/fasttree_paralogs/all_paralog_trees.nwk"
	shell:
		r"""
		mkdir -p results/fasttree_paralogs
		cat {input} > {output}
		"""


rule make_astral_pro_mapping:
	input:
		"results/fasttree_paralogs/all_paralog_trees.nwk"
	output:
		"results/astral_pro/mapping.tsv"
	shell:
		r"""
		mkdir -p results/astral_pro
		grep -oE '[^(),:;]+' {input} | sort -u \
		| awk -F'__' 'NF>=2 {{print $0 "\t" $1}}' \
		> {output}
		"""

rule concat_ortolog_trees:
	input:
		lambda wc: expand(
			"results/iqtree/{cluster}.treefile",
			cluster=clusters_from_msa_ortologs(wc)
		)
	output:
		"results/iqtree/all_ortolog_trees.nwk"
	shell:
		r"""
		cat {input} > {output}
		"""

rule run_consensus:
	input:
		"results/iqtree/all_ortolog_trees.nwk"
	output:
		contree="results/consensus/consensus.contree"
	threads: 1
	conda: "env/pipeline.yaml"
	shell:
		r"""
		mkdir -p results/consensus
		iqtree2 -con {input} \
			-pre results/consensus/consensus \
			-nt {threads}
		test -s {output.contree}
		"""


rule run_astral_ortologs:
	input:
		"results/iqtree/all_ortolog_trees.nwk"
	output:
		"results/astral/astral_species_tree.tre"
	threads: 1
	conda: "env/pipeline.yaml"
	shell:
		r"""
		mkdir -p results/astral
		astral -i {input} -o {output}
		"""


rule run_astral_pro_paralogs:
	input:
		trees="results/fasttree_paralogs/all_paralog_trees.nwk",
		mapping="results/astral_pro/mapping.tsv"
	output:
		"results/astral_pro/astral_pro_species_tree.tre"
	conda: "env/pipeline.yaml"
	shell:
		r"""
		mkdir -p results/astral_pro
		export JAVA_TOOL_OPTIONS="-Xmx8g"
		astral-pro3 -i {input.trees} -a {input.mapping} -o {output}
		"""


