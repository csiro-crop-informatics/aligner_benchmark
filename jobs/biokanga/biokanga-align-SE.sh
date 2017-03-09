#!/bin/bash
# treating the PE reads as if SE at default sensitivity and allowing for microInDels and splice junctions
# allowing 5 subs, 5bp microInDels, chimeric trimming 50%, allocate multimappers with uniques, 50Kbp splice junctions
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
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY
mkdir -p $OUTPUT_ALIGNMENT_PATH

# run biokanga
$BIOKANGA_PATH/$BIOKANGA_ALIGN_CMD -FbiokangaalignSE.log \
	-I$OUTPUT_INDEX_PATH/$GENOME_NAME.sfx \
        -m0 \
        -s10 \
        -M5 \
        -a5 \
        -c50 \
        -r3 \
        -A50000 \
        -i$READS_PATH/$READ_1_FILE -i$READS_PATH/$READ_2_FILE \
	-o$OUTPUT_ALIGNMENT_PATH/Aligned.out.sam \
	-T$NUM_THREAD

