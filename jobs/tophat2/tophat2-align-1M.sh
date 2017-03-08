#!/bin/bash
#BSUB -J tophat2-align-1M_job            # LSF job name
#BSUB -o tophat2-align-1M_job.%J.out     # Name of the job output file
#BSUB -e tophat2-align-1M_job.%J.error   # Name of the job error file

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


# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-1M/nocoveragesearch-bowtie2sensitive
mkdir -p $OUTPUT_ALIGNMENT_PATH


# run tophat2
$TOPHAT2_PATH/$TOPHAT2_ALIGN_CMD \
	--output-dir $OUTPUT_ALIGNMENT_PATH \
	--num-threads $NUM_THREAD \
	--mate-inner-dist $INNER_MATE_MEAN \
	--mate-std-dev $INNER_MATE_SD \
	--b2-very-sensitive \
	--GTF $GTF_PATH/$GTF_FILE \
	$OUTPUT_INDEX_PATH/$GENOME_NAME \
	$READS_PATH/1M.$READ_1_FILE $READS_PATH/1M.$READ_2_FILE
	
## --b2-very-sensitive		Option for bowtie 2, increase the sensitivity
