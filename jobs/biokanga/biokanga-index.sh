#!/bin/bash

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="biokanga"

# define library
LIBRARY=$DATASET

# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index
mkdir -p $OUTPUT_INDEX_PATH
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY

# build index
$BIOKANGA_PATH/$BIOKANGA_INDEX_CMD -T$NUM_THREAD -o$OUTPUT_INDEX_PATH/$GENOME_NAME.sfx -i$GENOME_PATH/$GENOME_FILE -r$GENOME_NAME
