#!/bin/bash
#BSUB -J olego-align-tune_job            # LSF job name
#BSUB -o olego-align-tune_job.%J.out     # Name of the job output file
#BSUB -e olego-align-tune_job.%J.error   # Name of the job error file

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="olego"

# define library
LIBRARY=$DATASET

# set parameters to tune
TOTAL_DIFF=$3
WORD_SIZE=$4
MAX_WORD_DIFF=$5
WORD_MAX_OVERLAP=$6
REP_ANCHOR=$7
MIN_ANCHOR=$8



# define output path
OUTPUT_ANNOTATION_PATH=$OUTPUT_ROOT_PATH/$TOOL/annotation
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-tune/tune-$TOTAL_DIFF-$WORD_SIZE-$MAX_WORD_DIFF-$WORD_MAX_OVERLAP-$REP_ANCHOR-$MIN_ANCHOR
mkdir -p $OUTPUT_ALIGNMENT_PATH

if [ $REP_ANCHOR != "on" ]
then
	# run olego
	$OLEGO_PATH/$OLEGO_ALIGN_CMD \
	        --output-file $OUTPUT_ALIGNMENT_PATH/output_1.sam \
	        --num-threads $NUM_THREAD \
	        --regression-model $OUTPUT_ANNOTATION_PATH/$OLEGO_REGRESSION_MODEL_FILE \
	        --verbose \
	        --junction-file $OUTPUT_ANNOTATION_PATH/$OLEGO_JUNCTION_FILE \
	        --max-total-diff $TOTAL_DIFF \
	        --word-size $WORD_SIZE \
	        --max-word-diff $MAX_WORD_DIFF \
	        --word-max-overlap $WORD_MAX_OVERLAP \
		--min-anchor $MIN_ANCHOR \
	        $OUTPUT_INDEX_PATH/$GENOME_NAME \
	        $READS_PATH/1M.$READ_1_FILE \

	$OLEGO_PATH/$OLEGO_ALIGN_CMD \
	        --output-file $OUTPUT_ALIGNMENT_PATH/output_2.sam \
	        --num-threads $NUM_THREAD \
	        --regression-model $OUTPUT_ANNOTATION_PATH/$OLEGO_REGRESSION_MODEL_FILE \
	        --verbose \
	        --junction-file $OUTPUT_ANNOTATION_PATH/$OLEGO_JUNCTION_FILE \
	        --max-total-diff $TOTAL_DIFF \
	        --word-size $WORD_SIZE \
	        --max-word-diff $MAX_WORD_DIFF \
	        --word-max-overlap $WORD_MAX_OVERLAP \
                --min-anchor $MIN_ANCHOR \
	        $OUTPUT_INDEX_PATH/$GENOME_NAME \
	        $READS_PATH/1M.$READ_2_FILE
else
    	# run olego
        $OLEGO_PATH/$OLEGO_ALIGN_CMD \
                --output-file $OUTPUT_ALIGNMENT_PATH/output_1.sam \
                --num-threads $NUM_THREAD \
                --regression-model $OUTPUT_ANNOTATION_PATH/$OLEGO_REGRESSION_MODEL_FILE \
                --verbose \
                --junction-file $OUTPUT_ANNOTATION_PATH/$OLEGO_JUNCTION_FILE \
                --max-total-diff $TOTAL_DIFF \
                --word-size $WORD_SIZE \
                --max-word-diff $MAX_WORD_DIFF \
                --word-max-overlap $WORD_MAX_OVERLAP \
                --min-anchor $MIN_ANCHOR \
		--allow-rep-anchor \
                $OUTPUT_INDEX_PATH/$GENOME_NAME \
                $READS_PATH/1M.$READ_1_FILE \

        $OLEGO_PATH/$OLEGO_ALIGN_CMD \
                --output-file $OUTPUT_ALIGNMENT_PATH/output_2.sam \
                --num-threads $NUM_THREAD \
                --regression-model $OUTPUT_ANNOTATION_PATH/$OLEGO_REGRESSION_MODEL_FILE \
                --verbose \
                --junction-file $OUTPUT_ANNOTATION_PATH/$OLEGO_JUNCTION_FILE \
                --max-total-diff $TOTAL_DIFF \
                --word-size $WORD_SIZE \
                --max-word-diff $MAX_WORD_DIFF \
                --word-max-overlap $WORD_MAX_OVERLAP \
                --min-anchor $MIN_ANCHOR \
		--allow-rep-anchor \
                $OUTPUT_INDEX_PATH/$GENOME_NAME \
                $READS_PATH/1M.$READ_2_FILE	

fi

perl $OLEGO_PATH/mergePEsam.pl -v $OUTPUT_ALIGNMENT_PATH/output_1.sam $OUTPUT_ALIGNMENT_PATH/output_2.sam $OUTPUT_ALIGNMENT_PATH/output.sam

# --non-denovo 		No de novo junction search
# --strand-mode		<value>
# --junction-file <file>

