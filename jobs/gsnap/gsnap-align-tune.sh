#!/bin/bash
#BSUB -J gsnap-align-tune_job            # LSF job name
#BSUB -o gsnap-align-tune_job.%J.out     # Name of the job output file
#BSUB -e gsnap-align-tune_job.%J.error   # Name of the job error file

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

# set parameters to tune
MAX_MISMATCHES=$3
INDEL_PENALITY=$4
GMAP_MIN_MATCH_LENGTH=$5
PAIR_EXPECT=$6
PAIR_DEV=$7

# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-tune/tune-$MAX_MISMATCHES-$INDEL_PENALITY-$GMAP_MIN_MATCH_LENGTH-$PAIR_EXPECT-$PAIR_DEV
mkdir -p $OUTPUT_ALIGNMENT_PATH


if [ $PAIR_EXPECT != "default" ]
then
	if [ $PAIR_DEV != "default" ]
	then
		# run gsnap
		$GSNAP_PATH/$GSNAP_ALIGN_CMD \
		        -D $OUTPUT_INDEX_PATH \
		        -d $GENOME_NAME \
	        	-A sam \
		        --max-mismatches $MAX_MISMATCHES \
		        --indel-penalty $INDEL_PENALITY \
		        --gmap-min-match-length $GMAP_MIN_MATCH_LENGTH \
		        --pairexpect $PAIR_EXPECT \
		        --pairdev $PAIR_DEV \
		        --merge-distant-samechr \
		        --ordered \
		        --novelsplicing 1 \
		        --use-splicing $GENOME_NAME.splicesites \
		        --nthreads $NUM_THREAD \
	        	--batch 5 \
		        --expand-offsets 1 \
		        $READS_PATH/$READ_1_FILE $READS_PATH/$READ_2_FILE \
	        	> $OUTPUT_ALIGNMENT_PATH/output.sam	
	else
		# run gsnap
                $GSNAP_PATH/$GSNAP_ALIGN_CMD \
                        -D $OUTPUT_INDEX_PATH \
                        -d $GENOME_NAME \
                       	-A sam \
                        --max-mismatches $MAX_MISMATCHES \
                        --indel-penalty $INDEL_PENALITY \
                        --gmap-min-match-length $GMAP_MIN_MATCH_LENGTH \
                        --pairexpect $PAIR_EXPECT \
                        --merge-distant-samechr \
                        --ordered \
                        --novelsplicing 1 \
                        --use-splicing $GENOME_NAME.splicesites \
                        --nthreads $NUM_THREAD \
                        --batch 5 \
                        --expand-offsets 1 \
                        $READS_PATH/$READ_1_FILE $READS_PATH/$READ_2_FILE \
                        > $OUTPUT_ALIGNMENT_PATH/output.sam

	fi
else
	if [ $PAIR_DEV != "default" ]
       	then
		# run gsnap
                $GSNAP_PATH/$GSNAP_ALIGN_CMD \
                        -D $OUTPUT_INDEX_PATH \
                       	-d $GENOME_NAME \
                        -A sam \
                        --max-mismatches $MAX_MISMATCHES \
                        --indel-penalty $INDEL_PENALITY \
                        --gmap-min-match-length $GMAP_MIN_MATCH_LENGTH \
                        --pairdev $PAIR_DEV \
                        --merge-distant-samechr \
                        --ordered \
                        --novelsplicing 1 \
                        --use-splicing $GENOME_NAME.splicesites \
                        --nthreads $NUM_THREAD \
                        --batch 5 \
                        --expand-offsets 1 \
	                $READS_PATH/$READ_1_FILE $READS_PATH/$READ_2_FILE \
                       	> $OUTPUT_ALIGNMENT_PATH/output.sam	
       	else
		# run gsnap
                $GSNAP_PATH/$GSNAP_ALIGN_CMD \
                        -D $OUTPUT_INDEX_PATH \
                        -d $GENOME_NAME \
                       	-A sam \
                        --max-mismatches $MAX_MISMATCHES \
                        --indel-penalty $INDEL_PENALITY \
                        --gmap-min-match-length $GMAP_MIN_MATCH_LENGTH \
                        --merge-distant-samechr \
                        --ordered \
                        --novelsplicing 1 \
                        --use-splicing $GENOME_NAME.splicesites \
                        --nthreads $NUM_THREAD \
                        --batch 5 \
                        --expand-offsets 1 \
                       	$READS_PATH/$READ_1_FILE $READS_PATH/$READ_2_FILE \
                       	> $OUTPUT_ALIGNMENT_PATH/output.sam
       	fi
fi
