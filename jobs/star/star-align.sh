#!/bin/bash
#BSUB -J star-align_job            # LSF job name
#BSUB -o star-align_job.%J.out     # Name of the job output file
#BSUB -e star-align_job.%J.error   # Name of the job error file

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="star"

# define library
LIBRARY=$DATASET


# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY
mkdir -p $OUTPUT_ALIGNMENT_PATH

# run star
$STAR_PATH/$STAR_ALIGN_CMD \
	--runThreadN $NUM_THREAD \
	--genomeDir $OUTPUT_INDEX_PATH \
	--readFilesIn $READS_PATH/$READ_1_FILE $READS_PATH/$READ_2_FILE \
	--outFileNamePrefix $OUTPUT_ALIGNMENT_PATH/$GENOME_NAME \
	--twopassMode Basic \
	--outSAMunmapped Within


#	--outSAMunmapped Within \
#	--outSAMorder PairedKeepInputOrder \
#	--twopassMode Basic \
#	--twopass1readsN 1

#	--genomeLoad LoadAndRemove

# --genomeLoad		#could speed-up??

