#!/bin/bash
#BSUB -J tophat2-align-tune_job            # LSF job name
#BSUB -o tophat2-align-tune_job.%J.out     # Name of the job output file
#BSUB -e tophat2-align-tune_job.%J.error   # Name of the job error file

if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh
fi
module load samtools-0.1.19
module load bowtie2-2.1.0


# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="tophat2"

# define library
LIBRARY=$DATASET

# set parameters to tune
NUM_MISMATCHES=$3
NUM_GAP_LENGTH=$4
NUM_EDIT_DIST=$5
NUM_REALIGN_EDIT_DIST=$6
NUM_INSERTION_LENGTH=$7
NUM_DELETION_LENGTH=$8
NUM_MULTIHITS=$9


# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-tune/nocoveragesearch-bowtie2sensitive-tune-$NUM_MISMATCHES-$NUM_GAP_LENGTH-$NUM_EDIT_DIST-$NUM_REALIGN_EDIT_DIST-$NUM_INSERTION_LENGTH-$NUM_DELETION_LENGTH-$NUM_MULTIHITS
mkdir -p $OUTPUT_ALIGNMENT_PATH


if [ $NUM_REALIGN_EDIT_DIST != "default" ]
then
	# run tophat2
	$TOPHAT2_PATH/$TOPHAT2_ALIGN_CMD \
        	--output-dir $OUTPUT_ALIGNMENT_PATH \
        	--num-threads $NUM_THREAD \
        	--mate-inner-dist $INNER_MATE_MEAN \
        	--mate-std-dev $INNER_MATE_SD \
        	--b2-very-sensitive \
        	--GTF $GTF_PATH/$GTF_FILE \
        	--read-mismatches $NUM_MISMATCHES \
        	--read-gap-length $NUM_GAP_LENGTH \
        	--read-edit-dist $NUM_EDIT_DIST \
       		--read-realign-edit-dist $NUM_REALIGN_EDIT_DIST \
        	--max-insertion-length $NUM_INSERTION_LENGTH \
        	--max-deletion-length $NUM_DELETION_LENGTH \
        	--max-multihits $NUM_MULTIHITS \
		$OUTPUT_INDEX_PATH/$GENOME_NAME \
		$READS_PATH/$READ_1_FILE $READS_PATH/$READ_2_FILE
else
	# run tophat2
	$TOPHAT2_PATH/$TOPHAT2_ALIGN_CMD \
		--output-dir $OUTPUT_ALIGNMENT_PATH \
                --num-threads $NUM_THREAD \
                --mate-inner-dist $INNER_MATE_MEAN \
                --mate-std-dev $INNER_MATE_SD \
                --b2-very-sensitive \
                --GTF $GTF_PATH/$GTF_FILE \
                --read-mismatches $NUM_MISMATCHES \
                --read-gap-length $NUM_GAP_LENGTH \
                --read-edit-dist $NUM_EDIT_DIST \
                --max-insertion-length $NUM_INSERTION_LENGTH \
                --max-deletion-length $NUM_DELETION_LENGTH \
                --max-multihits $NUM_MULTIHITS \
		$OUTPUT_INDEX_PATH/$GENOME_NAME \
		$READS_PATH/$READ_1_FILE $READS_PATH/$READ_2_FILE
fi

