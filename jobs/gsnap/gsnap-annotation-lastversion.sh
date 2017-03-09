#!/bin/bash
#BSUB -J gsnap-annotation_job            # LSF job name
#BSUB -o gsnap-annotation_job.%J.out     # Name of the job output file
#BSUB -e gsnap-annotation_job.%J.error   # Name of the job error file

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="gsnap"

# define library
LIBRARY=$DATASET
echo $LIBRARY


# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index-lastversion/$GENOME_NAME
mkdir -p $OUTPUT_INDEX_PATH
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY

cat $GTF_PATH/$GTF_FILE | $LASTVERSION_GSNAP_PATH/gtf_splicesites > foo.splicesites
#cat $GTF_PATH/$GTF_FILE | $LASTVERSION_GSNAP_PATH/gtf_introns > foo.introns

cat foo.splicesites | $LASTVERSION_GSNAP_PATH/iit_store -o $GENOME_NAME.splicesites

mkdir -p $OUTPUT_INDEX_PATH/$GENOME_NAME/$GENOME_NAME.maps/
mv $GENOME_NAME.splicesites.iit $OUTPUT_INDEX_PATH/$GENOME_NAME/$GENOME_NAME.maps/
rm foo.splicesites
