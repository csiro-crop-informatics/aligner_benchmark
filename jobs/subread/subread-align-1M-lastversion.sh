#!/bin/bash
#BSUB -J subread-align_job            # LSF job name
#BSUB -o subread-align_job.%J.out     # Name of the job output file
#BSUB -e subread-align_job.%J.error   # Name of the job error file

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="subread"

# define library
LIBRARY=$DATASET


# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index-lastversion/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-1M-lastversion
mkdir -p $OUTPUT_ALIGNMENT_PATH


# run subjunc
$LASTVERSION_SUBREAD_PATH/$SUBJUNC_ALIGN_CMD \
	-i $OUTPUT_INDEX_PATH/$GENOME_NAME \
	-r $READS_PATH/1M.$READ_1_FILE \
	-R $READS_PATH/1M.$READ_2_FILE \
	-T $NUM_THREAD \
	--allJunctions \
	--SAMoutput \
	-o $OUTPUT_ALIGNMENT_PATH/$GENOME_NAME

# -S < ff : fr : rf >	strand option
# --allJunctions 	report non-canonical exon-exon junc-tions and structural variants, in addition to the canonical exon-exon junctions (canonical donor/receptor sites detected).
