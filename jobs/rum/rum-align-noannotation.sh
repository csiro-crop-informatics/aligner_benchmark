#!/bin/bash
#BSUB -J rum-align-noannotation_job            # LSF job name
#BSUB -o rum-align-noannotation_job.%J.out     # Name of the job output file
#BSUB -e rum-align-noannotation_job.%J.error   # Name of the job error file


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
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-noannotation
mkdir -p $OUTPUT_ALIGNMENT_PATH


# run rum
$RUM_PATH/$RUM_ALIGN_CMD align \
	--index-dir $OUTPUT_INDEX_PATH \
	--name ${LIBRARY}_job \
	--output $OUTPUT_ALIGNMENT_PATH \
	--chunks $NUM_THREAD \
	$READS_PATH/$READ_1_FILE $READS_PATH/$READ_2_FILE \
	--verbose \
	--genome-only \
	--preserve-names

