#!/bin/bash
#BSUB -J mapsplice2-index_job            # LSF job name
#BSUB -o mapsplice2-index_job.%J.out     # Name of the job output file
#BSUB -e mapsplice2-index_job.%J.error   # Name of the job error file

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
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index-lastversion/$GENOME_NAME
mkdir -p $OUTPUT_INDEX_PATH
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-1M-lastversion

# build index
#cd $GENOME_PATH/files_nospace/
#$BOWTIE_PATH/$BOWTIE_INDEX_CMD $GENOME_FILES_LIST $OUTPUT_INDEX_PATH/$GENOME_NAME

cd $LASTVERSION_MAPSPLICE2_PATH/indexing

$LASTVERSION_MAPSPLICE2_PATH/indexing/buildIndex $GENOME_PATH/files_nospace/ $OUTPUT_INDEX_PATH
