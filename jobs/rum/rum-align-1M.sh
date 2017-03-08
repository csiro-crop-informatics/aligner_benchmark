#!/bin/bash
#BSUB -J rum-align_job            # LSF job name
#BSUB -o rum-align_job.%J.out     # Name of the job output file
#BSUB -e rum-align_job.%J.error   # Name of the job error file


# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="rum"

# define library
LIBRARY=$DATASET


# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-1M
mkdir -p $OUTPUT_ALIGNMENT_PATH


# run rum
$RUM_PATH/$RUM_ALIGN_CMD align \
	--index-dir $OUTPUT_INDEX_PATH \
	--name ${LIBRARY}_job \
	--output $OUTPUT_ALIGNMENT_PATH \
	--chunks $NUM_THREAD \
	$READS_PATH/1M.$READ_1_FILE $READS_PATH/1M.$READ_2_FILE \
	--verbose \
	--preserve-names

