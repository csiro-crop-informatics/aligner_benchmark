#!/bin/bash
#BSUB -J novoalign-align-1M_job            # LSF job name
#BSUB -o novoalign-align-1M_job.%J.out     # Name of the job output file
#BSUB -e novoalign-align-1M_job.%J.error   # Name of the job error file

if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh
fi

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
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index-lastversion/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-1M-lastversion
mkdir -p $OUTPUT_ALIGNMENT_PATH

# run novoalign
#$NOVOALIGN_PATH/$NOVOALIGN_ALIGN_CMD \
#	-d $OUTPUT_INDEX_PATH/$GENOME_NAME.nix \
#	-f $READS_PATH/$READ_1_FILE $READS_PATH/$READ_2_FILE \
#	> $OUTPUT_ALIGNMENT_PATH/output.sam \
#	2> $OUTPUT_ALIGNMENT_PATH/alignment.log

# run novoalign
$LASTVERSION_NOVOALIGN_PATH/$NOVOALIGN_ALIGN_CMD \
	-d $OUTPUT_INDEX_PATH/$GENOME_NAME.nix \
	-f $READS_PATH/1M.$READ_1_FILE $READS_PATH/1M.$READ_2_FILE \
	-F FA \
	-o SAM \
	-r All 10 \
	-i PE $FRAGMENT_LENGTH_MEAN,$FRAGMENT_LENGTH_SD \
	-v 0 70 70 "[>]([^:]*)" \
	> $OUTPUT_ALIGNMENT_PATH/output.tmp.sam \
	2> $OUTPUT_ALIGNMENT_PATH/alignment.log

# -k    Quality calibration. [not clear if it is useful for us, however not included in the unlicensed version]
# --lockidx     (causes crash)

#$SCRIPTS_PATH/fix-missing-qualities.sh $OUTPUT_ALIGNMENT_PATH

#java -Xmx64G -jar $USEQ_PATH/$USEQ_SAM_TRANSCRIPTOME_PARSER_CMD \
#	-f $OUTPUT_ALIGNMENT_PATH/output.tmp.fix.sam \
#	-a 50000 \
#	-n 100 \
#	-u \
#	-s $OUTPUT_ALIGNMENT_PATH/output.sam
