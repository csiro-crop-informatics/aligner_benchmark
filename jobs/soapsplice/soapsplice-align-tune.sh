#!/bin/bash
#BSUB -J soapsplice-align-tune_job            # LSF job name
#BSUB -o soapsplice-align-tune_job.%J.out     # Name of the job output file
#BSUB -e soapsplice-align-tune_job.%J.error   # Name of the job error file

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="soapsplice"

# define library
LIBRARY=$DATASET

# set parameters to tune
MISMATCHES=$3
INDEL=$4
TAIL=$5
SHORT_LENGTH=$6


# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-tune/tune-$MISMATCHES-$INDEL-$TAIL-$SHORT_LENGTH
mkdir -p $OUTPUT_ALIGNMENT_PATH


# run soapsplice
$SOAPSPLICE_PATH/$SOAPSPLICE_ALIGN_CMD \
	-d $OUTPUT_INDEX_PATH/$GENOME_NAME.index \
	-1 $READS_PATH/1M.$READ_1_FILE \
	-2 $READS_PATH/1M.$READ_2_FILE \
	-o $OUTPUT_ALIGNMENT_PATH/$GENOME_NAME \
	-p $NUM_THREAD \
        -f 2 \
        -l 0 \
        -I $FRAGMENT_LENGTH_MEAN \
	-m $MISMATCHES \
	-g $INDEL \
	-i $TAIL \
	-a $SHORT_LENGTH
