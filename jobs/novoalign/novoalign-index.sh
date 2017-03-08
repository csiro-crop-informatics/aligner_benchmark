#!/bin/bash
#BSUB -J novoalign-index_job            # LSF job name
#BSUB -o novoalign-index_job.%J.out     # Name of the job output file
#BSUB -e novoalign-index_job.%J.error   # Name of the job error file

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="novoalign"

# define library
LIBRARY=$DATASET


# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
mkdir -p $OUTPUT_INDEX_PATH
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY

# build index
$NOVOALIGN_PATH/$NOVOALIGN_INDEX_CMD \
	-t $NUM_THREAD \
	$OUTPUT_INDEX_PATH/$GENOME_NAME.nix \
	$OUTPUT_INDEX_PATH/${USEQ_ANNOTATION_PREFIX}Rad96Num100kMin10Splices.fasta \
	$OUTPUT_INDEX_PATH/${USEQ_ANNOTATION_PREFIX}Rad96Num100kMin10Transcripts.fasta \
	$OUTPUT_INDEX_PATH/maskedExons/geneMaskedGenome.fasta
