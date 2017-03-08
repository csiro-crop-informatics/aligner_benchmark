#!/bin/bash
#BSUB -J rum-align-tune_job            # LSF job name
#BSUB -o rum-align-tune_job.%J.out     # Name of the job output file
#BSUB -e rum-align-tune_job.%J.error   # Name of the job error file


# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="rum"

# define library
LIBRARY=$DATASET

# set parameters to tune
BLAT_MIN_IDENTITY=$3
BLAT_REP_MATCH=$4
BLAT_STEP_SIZE=$5
BLAT_TILE_SIZE=$6

# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-tune/tune-$BLAT_MIN_IDENTITY-$BLAT_REP_MATCH-$BLAT_STEP_SIZE-$BLAT_TILE_SIZE
mkdir -p $OUTPUT_ALIGNMENT_PATH


# run rum
$RUM_PATH/$RUM_ALIGN_CMD align \
	--index-dir $OUTPUT_INDEX_PATH \
	--name ${LIBRARY}_job \
	--output $OUTPUT_ALIGNMENT_PATH \
	--chunks $NUM_THREAD \
	$READS_PATH/1M.$READ_1_FILE $READS_PATH/1M.$READ_2_FILE \
	--verbose \
	--preserve-names \
	--blat-min-identity $BLAT_MIN_IDENTITY \
	--blat-rep-match $BLAT_REP_MATCH \
	--blat-step-size $BLAT_STEP_SIZE \
	--blat-tile-size $BLAT_TILE_SIZE
