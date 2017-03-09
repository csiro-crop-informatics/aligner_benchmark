#!/bin/bash
#BSUB -J star-align-tune_job            # LSF job name
#BSUB -o star-align-tune_job.%J.out     # Name of the job output file
#BSUB -e star-align-tune_job.%J.error   # Name of the job error file

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="star"

# define library
LIBRARY=$DATASET

# set parameters to tune
NUM_COLLAPSED_JUNCTIONS=$3
NUM_INSERTED_JUNCTIONS=$4
NUM_MULTIMAPPER=$5
NUM_FILTER_MISMATCHES=$6
RATIO_FILTER_MISMATCHES=$7
SEED_LENGTH=$8
OVERHANG=$9
END_ALIGNMENT_TYPE=${10}


# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-tune/tune-$NUM_COLLAPSED_JUNCTIONS-$NUM_INSERTED_JUNCTIONS-$NUM_MULTIMAPPER-$NUM_FILTER_MISMATCHES-$RATIO_FILTER_MISMATCHES-$SEED_LENGTH-$OVERHANG-$END_ALIGNMENT_TYPE
mkdir -p $OUTPUT_ALIGNMENT_PATH

# run star
$STAR_PATH/$STAR_ALIGN_CMD \
	--runThreadN $NUM_THREAD \
	--genomeDir $OUTPUT_INDEX_PATH \
	--readFilesIn $READS_PATH/1M.$READ_1_FILE $READS_PATH/1M.$READ_2_FILE \
	--outFileNamePrefix $OUTPUT_ALIGNMENT_PATH/$GENOME_NAME \
	--twopassMode Basic \
	--outSAMunmapped Within \
	--limitOutSJcollapsed $NUM_COLLAPSED_JUNCTIONS \
	--limitSjdbInsertNsj $NUM_INSERTED_JUNCTIONS \
	--outFilterMultimapNmax $NUM_MULTIMAPPER \
	--outFilterMismatchNmax $NUM_FILTER_MISMATCHES \
	--outFilterMismatchNoverLmax $RATIO_FILTER_MISMATCHES \
	--seedSearchStartLmax $SEED_LENGTH \
	--alignSJoverhangMin $OVERHANG \
	--alignEndsType $END_ALIGNMENT_TYPE

#	--outSAMunmapped Within \
#	--outSAMorder PairedKeepInputOrder \
#	--twopass1readsN 1
#	--genomeLoad LoadAndRemove

# --genomeLoad		#could speed-up??

