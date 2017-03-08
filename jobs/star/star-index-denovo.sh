#!/bin/bash
#BSUB -J star-index-denovo_job            # LSF job name
#BSUB -o star-index-denovo_job.%J.out     # Name of the job output file
#BSUB -e star-index-denovo_job.%J.error   # Name of the job error file

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="star"

# define library
LIBRARY=$DATASET


# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME-denovo
mkdir -p $OUTPUT_INDEX_PATH
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-denovo

# build index
$STAR_PATH/$STAR_INDEX_CMD \
        --runThreadN $NUM_THREAD \
        --runMode genomeGenerate \
        --genomeDir $OUTPUT_INDEX_PATH \
        --genomeFastaFiles $GENOME_PATH/$GENOME_FILE
