#!/bin/bash
#BSUB -J crac-align_job            # LSF job name
#BSUB -o crac-align_job.%J.out     # Name of the job output file
#BSUB -e crac-align_job.%J.error   # Name of the job error file

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="crac"

# define library
LIBRARY=$DATASET

# load crac library
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CRAC_LIB

# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-1M
mkdir -p $OUTPUT_ALIGNMENT_PATH

# run crac
$CRAC_PATH/$CRAC_ALIGN_CMD \
	-i $OUTPUT_INDEX_PATH/$GENOME_NAME \
	-k 22 \
	-r $READS_PATH/1M.$READ_1_FILE $READS_PATH/1M.$READ_2_FILE \
	--sam $OUTPUT_ALIGNMENT_PATH/output.sam \
	--reads-length 100 \
        --summary $OUTPUT_ALIGNMENT_PATH/summary.txt \
	--nb-threads $NUM_THREAD

#	--reads-length|-m <read_length> (If the read length is fixed, is deeply recommend to specify the read length. CRAC will therefore be much faster)
#	--nb-threads <num_thread>
#	--stranded
