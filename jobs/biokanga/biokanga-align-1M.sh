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
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-1M
mkdir -p $OUTPUT_ALIGNMENT_PATH


# run biokanga
$BIOKANGA_PATH/$BIOKANGA_ALIGN_CMD \
        -I$OUTPUT_INDEX_PATH/$GENOME_NAME.sfx \
        -M5 \
        -U1 \
        -c50 \
        -D10000 \
        -i$READS_PATH/1M.$READ_1_FILE -u$READS_PATH/1M.$READ_2_FILE \
        -o$OUTPUT_ALIGNMENT_PATH/Aligned.out.sam \
        -T$NUM_THREAD
        
