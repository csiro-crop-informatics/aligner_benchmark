#!/bin/bash
#BSUB -J mapsplice2-align_job            # LSF job name
#BSUB -o mapsplice2-align_job.%J.out     # Name of the job output file
#BSUB -e mapsplice2-align_job.%J.error   # Name of the job error file

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="mapsplice2"

# define library
LIBRARY=$DATASET


# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-1M
mkdir -p $OUTPUT_ALIGNMENT_PATH


# run mapsplice2
python $MAPSPLICE2_PATH/$MAPSPLICE2_ALIGN_CMD \
	--threads $NUM_THREAD \
	--non-canonical-double-anchor \
	--output $OUTPUT_ALIGNMENT_PATH \
	-c $GENOME_PATH/files_nospace \
	-x $OUTPUT_INDEX_PATH/$GENOME_NAME \
	-1 $READS_PATH/1M.$READ_1_FILE_CHANGE_ID \
	-2 $READS_PATH/1M.$READ_2_FILE_CHANGE_ID

