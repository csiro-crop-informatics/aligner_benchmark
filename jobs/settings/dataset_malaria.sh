#!/bin/bash


# Dataset name
DATASET="dataset_malaria_pfal"


# Genome
GENOME_PATH="/project/itmatlab/aligner_benchmark/dataset/malaria/genome"
GENOME_FILE="pfal.fa"
GENOME_FILES_LIST="chr10.fa,chr5.fa,chr11.fa,chr12.fa,chr1.fa,chr13.fa,chr2.fa,chr6.fa,chr14.fa,chr7.fa,chr3.fa,chr8.fa,chrX.fa,chrM.fa,chr4.fa,chr9.fa"
GENOME_NAME="pfal"


# Annotation
GTF_PATH="/project/itmatlab/aligner_benchmark/dataset/malaria/annotation"
GTF_FILE="simulator_config_geneinfo_pfal_GTF.gtf"
BED_PATH="/project/itmatlab/aligner_benchmark/dataset/malaria/annotation"
BED_FILE="simulator_config_geneinfo_pfal_BED.bed"


# Olego junctions file
OLEGO_JUNCTION_FILE="pfal.junctions.bed"

# Olego regression model file
OLEGO_REGRESSION_MODEL_FILE="pfal.regression_model.cfg"

# Olego	regression model file prefix
OLEGO_REGRESSION_MODEL_FILE_PREFIX="pfal.regression_model"

# USEQ annotation file
USEQ_ANNOTATION_PATH="/project/itmatlab/aligner_benchmark/dataset/malaria/annotation"
USEQ_ANNOTATION_FILE="simulator_config_geneinfo_pfal_refFlat.txt"
USEQ_ANNOTATION_PREFIX="simulator_config_geneinfo_pfal_refFlat"

# RUM index name
RUM_INDEX_NAME="pfalciparum"
