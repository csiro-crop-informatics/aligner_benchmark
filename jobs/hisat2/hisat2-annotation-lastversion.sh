#!/bin/bash
#BSUB -J hisat-annotation_job            # LSF job name
#BSUB -o hisat-annotation_job.%J.out     # Name of the job output file
#BSUB -e hisat-annotation_job.%J.error   # Name of the job error file

if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh
fi
module load python-3.4.2


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
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index-lastversion/$GENOME_NAME
mkdir -p $OUTPUT_INDEX_PATH
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY

$LASTVERSION_HISAT2_PATH/hisat2_extract_splice_sites.py $GTF_PATH/$GTF_FILE > $GENOME_NAME.splicesites.txt
$LASTVERSION_HISAT2_PATH/hisat2_extract_exons.py $GTF_PATH/$GTF_FILE > $GENOME_NAME.exons.txt
mv $GENOME_NAME.splicesites.txt $OUTPUT_INDEX_PATH/$GENOME_NAME.splicesites.txt
mv $GENOME_NAME.exons.txt $OUTPUT_INDEX_PATH/$GENOME_NAME.exons.txt
