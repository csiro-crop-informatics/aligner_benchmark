#!/bin/bash
#BSUB -J hisat2-index-denovo_job            # LSF job name
#BSUB -o hisat2-index-denovo_job.%J.out     # Name of the job output file
#BSUB -e hisat2-index-denovo_job.%J.error   # Name of the job error file

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="hisat2"

# define library
LIBRARY=$DATASET


# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME-denovo
mkdir -p $OUTPUT_INDEX_PATH
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY

# build index
$HISAT2_PATH/$HISAT2_INDEX_CMD \
	-p $NUM_THREAD \
	$GENOME_PATH/$GENOME_FILE \
	$OUTPUT_INDEX_PATH/$GENOME_NAME
