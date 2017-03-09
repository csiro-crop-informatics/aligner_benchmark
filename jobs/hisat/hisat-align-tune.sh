#!/bin/bash
#BSUB -J hisat-align-tune_job            # LSF job name
#BSUB -o hisat-align-tune_job.%J.out     # Name of the job output file
#BSUB -e hisat-align-tune_job.%J.error   # Name of the job error file

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="hisat"

# define library
LIBRARY=$DATASET

# set parameters to tune
MODE=$3
NUM_MISMATCH=$4
SEED_LENGTH=$5
SEED_INTERVAL=$6
SEED_EXTENSION=$7
RE_SEED=$8
PENALITY_NONCANONICAL=$9
MAX_MISMATCH_PENALITY=${10}
MIN_MISMATCH_PENALITY=${11}

# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-tune/tune-$MODE-$NUM_MISMATCH-$SEED_LENGTH-$SEED_INTERVAL-$SEED_EXTENSION-$RE_SEED-$PENALITY_NONCANONICAL-$MAX_MISMATCH_PENALITY-$MIN_MISMATCH_PENALITY
mkdir -p $OUTPUT_ALIGNMENT_PATH

if [ $MODE != "local" ]
then
	# run hisat (in end-to-end mode)
	$HISAT_PATH/$HISAT_ALIGN_CMD \
		--threads $NUM_THREAD \
		--end-to-end \
		-N $NUM_MISMATCH \
		-L $SEED_LENGTH \
		-i S,1,${SEED_INTERVAL} \
		-D $SEED_EXTENSION \
		-R $RE_SEED \
		--pen-noncansplice $PENALITY_NONCANONICAL \
		--mp $MAX_MISMATCH_PENALITY,$MIN_MISMATCH_PENALITY \
		--time \
		--reorder \
		--known-splicesite-infile $OUTPUT_INDEX_PATH/$GENOME_NAME.splicesites.txt \
		--novel-splicesite-outfile $OUTPUT_ALIGNMENT_PATH/$GENOME_NAME.splicesites.novel.txt \
		--novel-splicesite-infile $OUTPUT_ALIGNMENT_PATH/$GENOME_NAME.splicesites.novel.txt \
		-f \
		-x $OUTPUT_INDEX_PATH/$GENOME_NAME \
		-1 $READS_PATH/1M.$READ_1_FILE \
		-2 $READS_PATH/1M.$READ_2_FILE \
		-S $OUTPUT_ALIGNMENT_PATH/output.sam
else
	# run hisat (in local mode)
        $HISAT_PATH/$HISAT_ALIGN_CMD \
                --threads $NUM_THREAD \
                -N $NUM_MISMATCH \
                -L $SEED_LENGTH \
                -i S,1,${SEED_INTERVAL} \
                -D $SEED_EXTENSION \
                -R $RE_SEED \
                --pen-noncansplice $PENALITY_NONCANONICAL \
                --mp $MAX_MISMATCH_PENALITY,$MIN_MISMATCH_PENALITY \
                --time \
                --reorder \
                --known-splicesite-infile $OUTPUT_INDEX_PATH/$GENOME_NAME.splicesites.txt \
		--novel-splicesite-outfile $OUTPUT_ALIGNMENT_PATH/$GENOME_NAME.splicesites.novel.txt \
                --novel-splicesite-infile $OUTPUT_ALIGNMENT_PATH/$GENOME_NAME.splicesites.novel.txt \
                -f \
                -x $OUTPUT_INDEX_PATH/$GENOME_NAME \
                -1 $READS_PATH/1M.$READ_1_FILE \
                -2 $READS_PATH/1M.$READ_2_FILE \
                -S $OUTPUT_ALIGNMENT_PATH/output.sam
fi

#                --novel-splicesite-outfile $OUTPUT_ALIGNMENT_PATH/$GENOME_NAME.splicesites.novel.txt \

