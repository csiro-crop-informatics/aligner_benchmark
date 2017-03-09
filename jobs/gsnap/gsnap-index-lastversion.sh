#!/bin/bash
#BSUB -J gsnap-index_job            # LSF job name
#BSUB -o gsnap-index_job.%J.out     # Name of the job output file
#BSUB -e gsnap-index_job.%J.error   # Name of the job error file

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="gsnap"

# define library
LIBRARY=$DATASET
echo $LIBRARY


# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index-lastversion/$GENOME_NAME
mkdir -p $OUTPUT_INDEX_PATH
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY

# build index
$LASTVERSION_GSNAP_PATH/$GSNAP_INDEX_CMD -D $OUTPUT_INDEX_PATH -d $GENOME_NAME $GENOME_PATH/$GENOME_FILE
