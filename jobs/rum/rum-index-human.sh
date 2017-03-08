#!/bin/bash
#BSUB -J rum-index_job            # LSF job name
#BSUB -o rum-index_job.%J.out     # Name of the job output file
#BSUB -e rum-index_job.%J.error   # Name of the job error file


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
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index
mkdir -p $OUTPUT_INDEX_PATH
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY

# build index
perl $RUM_PATH/create_indexes_from_ucsc.pl ucsc.hg19_genome.txt ucsc.hg19_info_files
