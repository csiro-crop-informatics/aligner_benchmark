#!/bin/bash


# Dataset name
DATASET="dataset_human_hg19_RefSeq_t1r1-40adapter"


# Genome
GENOME_PATH="/project/itmatlab/aligner_benchmark/dataset/human/genome"
GENOME_FILE="ucsc.hg19.fa"
GENOME_FILES_LIST="chr10.fa,chr19.fa,chr5.fa,chr11.fa,chr12.fa,chr1.fa,chr13.fa,chr2.fa,chr6.fa,chr14.fa,chr15.fa,chr20.fa,chr16.fa,chr21.fa,chr17.fa,chr22.fa,chr7.fa,chr3.fa,chr8.fa,chrX.fa,chrY.fa,chr4.fa,chr18.fa,chr9.fa"
GENOME_NAME="ucsc.hg19"


# Annotation
GTF_PATH="/project/itmatlab/aligner_benchmark/dataset/human/annotation"
GTF_FILE="hg19_RefSeq_genes_ucsc.gtf"


# Stats
INNER_MATE_MEAN=21
INNER_MATE_SD=41
FRAGMENT_LENGTH_MEAN=221
FRAGMENT_LENGTH_SD=$INNER_MATE_SD
FRAGMENT_LENGTH_MIN=100
FRAGMENT_LENGTH_MAX=500


# Reads
READS_PATH="/project/itmatlab/aligner_benchmark/dataset/human/dataset_t1r1"
READS_FILE="40.adapter.simulated_reads_HG19t1r1.fa"
READS_FILE_CHANGE_ID="40.adapter.simulated_reads_HG19t1r1.changeID.fa"
READ_1_FILE="40.adapter.simulated_reads_HG19t1r1.forward.fa"
READ_1_FILE_CHANGE_ID="40.adapter.simulated_reads_HG19t1r1.forward.changeID.fa"
READ_2_FILE="40.adapter.simulated_reads_HG19t1r1.reverse.fa"
READ_2_FILE_CHANGE_ID="40.adapter.simulated_reads_HG19t1r1.reverse.changeID.fa"


# Olego junctions file
OLEGO_JUNCTION_FILE="hg19.intron.hmr.bed"

# Olego regression model file
OLEGO_REGRESSION_MODEL_FILE="hg.cfg"

# Olego regression model file prefix
OLEGO_REGRESSION_MODEL_FILE_PREFIX=""

# RUM index name
RUM_INDEX_NAME="hg19"

# USEQ annotation file
USEQ_ANNOTATION_PATH="/project/itmatlab/aligner_benchmark/dataset/human/annotation"
USEQ_ANNOTATION_FILE="hg19_RefSeq_refFlat_ucsc.txt"
USEQ_ANNOTATION_PREFIX="hg19_RefSeq_refFlat_ucsc"

# Compare to truth
CIG_FILE="40.adapter.simulated_reads_HG19t1r1.cig"
TRANSCRIPTS="simulated_reads_transcripts_HG19t1r1.txt"
JUNCTIONS_CROSSED="simulated_reads_junctions-crossed_HG19t1r1.txt"

