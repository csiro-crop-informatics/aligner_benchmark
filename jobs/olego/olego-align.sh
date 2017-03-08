#!/bin/bash
#BSUB -J olego-align_job            # LSF job name
#BSUB -o olego-align_job.%J.out     # Name of the job output file
#BSUB -e olego-align_job.%J.error   # Name of the job error file

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="olego"

# define library
LIBRARY=$DATASET


# define output path
OUTPUT_ANNOTATION_PATH=$OUTPUT_ROOT_PATH/$TOOL/annotation
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY
mkdir -p $OUTPUT_ALIGNMENT_PATH


# run olego
$OLEGO_PATH/$OLEGO_ALIGN_CMD \
	--output-file $OUTPUT_ALIGNMENT_PATH/output_1.sam \
	--num-threads $NUM_THREAD \
	--regression-model $OUTPUT_ANNOTATION_PATH/$OLEGO_REGRESSION_MODEL_FILE \
	--verbose \
	--junction-file $OUTPUT_ANNOTATION_PATH/$OLEGO_JUNCTION_FILE \
	$OUTPUT_INDEX_PATH/$GENOME_NAME \
	$READS_PATH/$READ_1_FILE \


$OLEGO_PATH/$OLEGO_ALIGN_CMD \
        --output-file $OUTPUT_ALIGNMENT_PATH/output_2.sam \
        --num-threads $NUM_THREAD \
        --regression-model $OUTPUT_ANNOTATION_PATH/$OLEGO_REGRESSION_MODEL_FILE \
        --verbose \
        --junction-file $OUTPUT_ANNOTATION_PATH/$OLEGO_JUNCTION_FILE \
        $OUTPUT_INDEX_PATH/$GENOME_NAME \
        $READS_PATH/$READ_2_FILE

perl $OLEGO_PATH/mergePEsam.pl -v $OUTPUT_ALIGNMENT_PATH/output_1.sam $OUTPUT_ALIGNMENT_PATH/output_2.sam $OUTPUT_ALIGNMENT_PATH/output.sam

# --non-denovo 		No de novo junction search
# --strand-mode		<value>
# --junction-file <file>

