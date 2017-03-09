#!/bin/bash
#BSUB -J contextmap2-index_job            # LSF job name
#BSUB -o contextmap2-index_job.%J.out     # Name of the job output file 
#BSUB -e contextmap2-index_job.%J.error   # Name of the job error file

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="contextmap2"

# define library
LIBRARY=$DATASET

# define output	path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
mkdir -p $OUTPUT_INDEX_PATH
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY

# ALIGNER must be one of "bwa", "bowtie1" or "bowtie2"
ALIGNER="bwa"
ALIGNER_BIN_PATH=""
INDEXER_BIN_PATH=""


# build index
echo "build $ALIGNER index"
if [ "${ALIGNER}" = bowtie1 ]
then
	echo "build Bowtie index"
	$BOWTIE_PATH/$BOWTIE_INDEX_CMD -f $GENOME_PATH/$GENOME_FILE $OUTPUT_INDEX_PATH/$GENOME_NAME
	ALIGNER_BIN_PATH=$BOWTIE_PATH/$BOWTIE_ALIGN_CMD
	INDEXER_BIN_PATH=$BOWTIE_PATH/$BOWTIE_INDEX_CMD
else
	if [ "${ALIGNER}" = bowtie2 ]
	then 
		echo "build Bowtie 2 index"
        	$BOWTIE2_PATH/$BOWTIE2_INDEX_CMD -f $GENOME_PATH/$GENOME_FILE $OUTPUT_INDEX_PATH/$GENOME_NAME
		ALIGNER_BIN_PATH=$BOWTIE2_PATH/$BOWTIE2_ALIGN_CMD
		INDEXER_BIN_PATH=$BOWTIE2_PATH/$BOWTIE2_INDEX_CMD
	else
		if [ "${ALIGNER}" = bwa ]
		then 
			echo "build BWA index"
		        $BWA_PATH/$BWA_INDEX_CMD -a bwtsw -p $OUTPUT_INDEX_PATH/$GENOME_NAME $GENOME_PATH/$GENOME_FILE
			ALIGNER_BIN_PATH=$BWA_PATH/$BWA_ALIGN_CMD
			INDEXER_BIN_PATH=$BWA_PATH/$BWA_INDEX_CMD
		fi
	fi
fi
