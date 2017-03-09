#!/bin/bash
#BSUB -J subread-align-tune_job            # LSF job name
#BSUB -o subread-align-tune_job.%J.out     # Name of the job output file
#BSUB -e subread-align-tune_job.%J.error   # Name of the job error file

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="subread"

# define library
LIBRARY=$DATASET

# set parameters to tune
MIN_FRAGMENT_LENGTH=$3
INDEL=$4
NUM_HIT_SUBREADS=$5
MISMATCHES=$6
NUM_EXTRACTED_SUBREADS=$7
NUM_HIT_PAIR_SUBREADS=$8
COMPLEX_INDELS=$9

# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-tune/tune-$MIN_FRAGMENT_LENGTH-$INDEL-$NUM_HIT_SUBREADS-$MISMATCHES-$NUM_EXTRACTED_SUBREADS-$NUM_HIT_PAIR_SUBREADS-$COMPLEX_INDELS
mkdir -p $OUTPUT_ALIGNMENT_PATH

if [ $COMPLEX_INDELS != "on" ]
then
        # run subjunc (with no --complexIndels)
        $SUBREAD_PATH/$SUBJUNC_ALIGN_CMD \
                -i $OUTPUT_INDEX_PATH/$GENOME_NAME \
                -r $READS_PATH/1M.$READ_1_FILE \
                -R $READS_PATH/1M.$READ_2_FILE \
                -T $NUM_THREAD \
                --allJunctions \
                -o $OUTPUT_ALIGNMENT_PATH/$GENOME_NAME \
                -d $MIN_FRAGMENT_LENGTH \
                -I $INDEL \
                -m $NUM_HIT_SUBREADS \
                -M $MISMATCHES \
                -n $NUM_EXTRACTED_SUBREADS \
                -p $NUM_HIT_PAIR_SUBREADS
else
	# run subjunc (with --complexIndels)
	$SUBREAD_PATH/$SUBJUNC_ALIGN_CMD \
		-i $OUTPUT_INDEX_PATH/$GENOME_NAME \
		-r $READS_PATH/1M.$READ_1_FILE \
		-R $READS_PATH/1M.$READ_2_FILE \
		-T $NUM_THREAD \
		--allJunctions \
		-o $OUTPUT_ALIGNMENT_PATH/$GENOME_NAME \
		-d $MIN_FRAGMENT_LENGTH \
		-I $INDEL \
		-m $NUM_HIT_SUBREADS \
		-M $MISMATCHES \
		-n $NUM_EXTRACTED_SUBREADS \
		-p $NUM_HIT_PAIR_SUBREADS \
		--complexIndels
fi


# -S < ff : fr : rf >	strand option
# --allJunctions 	report non-canonical exon-exon junc-tions and structural variants, in addition to the canonical exon-exon junctions (canonical donor/receptor sites detected).
