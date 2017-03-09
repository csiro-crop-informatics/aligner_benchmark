#!/bin/bash
#BSUB -J subread-index_job            # LSF job name
#BSUB -o subread-index_job.%J.out     # Name of the job output file
#BSUB -e subread-index_job.%J.error   # Name of the job error file

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="subread"

# define library
LIBRARY=$DATASET


# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
mkdir -p $OUTPUT_INDEX_PATH
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY

# build index
$SUBREAD_PATH/$SUBREAD_INDEX_CMD \
	$GENOME_PATH/$GENOME_FILE \
	-B \
	-o $OUTPUT_INDEX_PATH/$GENOME_NAME

	
#-B 	do not split index (more memory but faster alignment)
