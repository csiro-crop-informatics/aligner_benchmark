#!/bin/bash
# treating the PE reads as if PE at default sensitivity and relying on chimeric trimming and large insert sizes to cross splice junctions
# aligning default sensitivity, output as SAM, 5 subs allowed, 2 indeterminates, chimeric 50%, inserts upto 50Kbp
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
$BIOKANGA_PATH/$BIOKANGA_ALIGN_CMD -Fbiokangaalign_PE.log \
	-I$OUTPUT_INDEX_PATH/$GENOME_NAME.sfx \
        -m0 \
        -s5 \
        -M5 \
        -n2 \
        -c50 \
        -U2 \
        -D50000 \
        -i$READS_PATH/$READ_1_FILE -u$READS_PATH/$READ_2_FILE \
	-o$OUTPUT_ALIGNMENT_PATH/Aligned.out.sam \
	-T$NUM_THREAD

