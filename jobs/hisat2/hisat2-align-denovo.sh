#!/bin/bash
#BSUB -J hisat2-align-noannotation_job            # LSF job name
#BSUB -o hisat2-align-noannotation_job.%J.out     # Name of the job output file
#BSUB -e hisat2-align-noannotation_job.%J.error   # Name of the job error file

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="hisat2"

# define library
LIBRARY=$DATASET


# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME-denovo
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-noannotation
mkdir -p $OUTPUT_ALIGNMENT_PATH

# run hisat2
$HISAT2_PATH/$HISAT2_ALIGN_CMD \
	--threads $NUM_THREAD \
	--time \
	--reorder \
	-f \
	-x $OUTPUT_INDEX_PATH/$GENOME_NAME \
	-1 $READS_PATH/$READ_1_FILE \
	-2 $READS_PATH/$READ_2_FILE \
	-S $OUTPUT_ALIGNMENT_PATH/output.sam \
	
#        --time \
#        --reorder \
#        --known-splicesite-infile $OUTPUT_INDEX_PATH/$GENOME_NAME.splicesites.txt \
#        --novel-splicesite-outfile $OUTPUT_ALIGNMENT_PATH/$GENOME_NAME.splicesites.novel.txt \
#        --novel-splicesite-infile $OUTPUT_ALIGNMENT_PATH/$GENOME_NAME.splicesites.novel.txt \

