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

# set parameters to tune
K=$3
READ_LENGTH=$4
NO_AMBIGUITY=$5
MAX_LOCS=$6
MIN_PERCENT_SINGLE_LOC=$7
MIN_PERCENT_MULTIPLE_LOC=$8
MIN_REPETITION=$9
MIN_PERCENT_REPETITION_LOC=${10}


# load crac library
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CRAC_LIB

# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-tune/tune-$K-$READ_LENGTH-$NO_AMBIGUITY-$MAX_LOCS-$MIN_PERCENT_SINGLE_LOC-$MIN_PERCENT_MULTIPLE_LOC-$MIN_REPETITION-$MIN_PERCENT_REPETITION_LOC
mkdir -p $OUTPUT_ALIGNMENT_PATH

if [ $NO_AMBIGUITY != "on" ]
then
	if [ $MAX_LOCS != "default" ]
	then
		if [ $MIN_PERCENT_REPETITION_LOC != "default" ]
        	then
			# run crac (without --no-ambiguity, without --max-locs default, without --min-percent-repetition-loc default)
			$CRAC_PATH/$CRAC_ALIGN_CMD \
				-i $OUTPUT_INDEX_PATH/$GENOME_NAME \
				-k $K \
				-r $READS_PATH/1M.$READ_1_FILE $READS_PATH/1M.$READ_2_FILE \
				--sam $OUTPUT_ALIGNMENT_PATH/output.sam \
				--reads-length $READ_LENGTH \
				--max-locs $MAX_LOCS \
				--min-percent-single-loc $MIN_PERCENT_SINGLE_LOC \
				--min-percent-multiple-loc $MIN_PERCENT_MULTIPLE_LOC \
				--min-repetition $MIN_REPETITION \
				--min-percent-repetition-loc $MIN_PERCENT_REPETITION_LOC \
			        --summary $OUTPUT_ALIGNMENT_PATH/summary.txt \
				--nb-threads $NUM_THREAD
		else
			# run crac (without --no-ambiguity, without --max-locs default, with --min-percent-repetition-loc default)
                        $CRAC_PATH/$CRAC_ALIGN_CMD \
                                -i $OUTPUT_INDEX_PATH/$GENOME_NAME \
                               	-k $K \
                                -r $READS_PATH/1M.$READ_1_FILE $READS_PATH/1M.$READ_2_FILE \
                                --sam $OUTPUT_ALIGNMENT_PATH/output.sam \
                                --reads-length $READ_LENGTH \
                                --max-locs $MAX_LOCS \
                                --min-percent-single-loc $MIN_PERCENT_SINGLE_LOC \
                                --min-percent-multiple-loc $MIN_PERCENT_MULTIPLE_LOC \
                                --min-repetition $MIN_REPETITION \
                                --summary $OUTPUT_ALIGNMENT_PATH/summary.txt \
                                --nb-threads $NUM_THREAD
		fi
	else
		if [ $MIN_PERCENT_REPETITION_LOC != "default" ]
	        then
			# run crac (without --no-ambiguity, with --max-locs default, without --min-percent-repetition-loc default)
                	$CRAC_PATH/$CRAC_ALIGN_CMD \
                        	-i $OUTPUT_INDEX_PATH/$GENOME_NAME \
	                       	-k $K \
	                        -r $READS_PATH/1M.$READ_1_FILE $READS_PATH/1M.$READ_2_FILE \
	                        --sam $OUTPUT_ALIGNMENT_PATH/output.sam \
	                        --reads-length $READ_LENGTH \
        	                --min-percent-single-loc $MIN_PERCENT_SINGLE_LOC \
	                        --min-percent-multiple-loc $MIN_PERCENT_MULTIPLE_LOC \
	                        --min-repetition $MIN_REPETITION \
	                        --min-percent-repetition-loc $MIN_PERCENT_REPETITION_LOC \
	                        --summary $OUTPUT_ALIGNMENT_PATH/summary.txt \
	                        --nb-threads $NUM_THREAD
		else
			# run crac (without --no-ambiguity, with --max-locs default, with --min-percent-repetition-loc default)
                        $CRAC_PATH/$CRAC_ALIGN_CMD \
                                -i $OUTPUT_INDEX_PATH/$GENOME_NAME \
                               	-k $K \
                                -r $READS_PATH/1M.$READ_1_FILE $READS_PATH/1M.$READ_2_FILE \
                                --sam $OUTPUT_ALIGNMENT_PATH/output.sam \
                                --reads-length $READ_LENGTH \
                                --min-percent-single-loc $MIN_PERCENT_SINGLE_LOC \
                                --min-percent-multiple-loc $MIN_PERCENT_MULTIPLE_LOC \
                                --min-repetition $MIN_REPETITION \
                                --summary $OUTPUT_ALIGNMENT_PATH/summary.txt \
                                --nb-threads $NUM_THREAD
		fi
	fi
else
	if [ $MAX_LOCS != "default" ]
        then
		if [ $MIN_PERCENT_REPETITION_LOC != "default" ]
        	then
			# run crac (with --no-ambiguity, without --max-locs default, without --min-percent-repetition-loc default)
			$CRAC_PATH/$CRAC_ALIGN_CMD \
			        -i $OUTPUT_INDEX_PATH/$GENOME_NAME \
			       	-k $K \
			        -r $READS_PATH/1M.$READ_1_FILE $READS_PATH/1M.$READ_2_FILE \
			        --sam $OUTPUT_ALIGNMENT_PATH/output.sam \
				--no-ambiguity \
			        --reads-length $READ_LENGTH \
			        --max-locs $MAX_LOCS \
			        --min-percent-single-loc $MIN_PERCENT_SINGLE_LOC \
			        --min-percent-multiple-loc $MIN_PERCENT_MULTIPLE_LOC \
			        --min-repetition $MIN_REPETITION \
			        --min-percent-repetition-loc $MIN_PERCENT_REPETITION_LOC \
			        --summary $OUTPUT_ALIGNMENT_PATH/summary.txt \
			        --nb-threads $NUM_THREAD
		else
			# run crac (with --no-ambiguity, without --max-locs default, with --min-percent-repetition-loc default)
			$CRAC_PATH/$CRAC_ALIGN_CMD \
			        -i $OUTPUT_INDEX_PATH/$GENOME_NAME \
			       	-k $K \
			        -r $READS_PATH/1M.$READ_1_FILE $READS_PATH/1M.$READ_2_FILE \
			        --sam $OUTPUT_ALIGNMENT_PATH/output.sam \
				--no-ambiguity \
			        --reads-length $READ_LENGTH \
			        --max-locs $MAX_LOCS \
			        --min-percent-single-loc $MIN_PERCENT_SINGLE_LOC \
			        --min-percent-multiple-loc $MIN_PERCENT_MULTIPLE_LOC \
			        --min-repetition $MIN_REPETITION \
			        --summary $OUTPUT_ALIGNMENT_PATH/summary.txt \
			        --nb-threads $NUM_THREAD
		fi
	else
		if [ $MIN_PERCENT_REPETITION_LOC != "default" ]
        	then
			# run crac (with --no-ambiguity, with --max-locs default, without --min-percent-repetition-loc default)
		        $CRAC_PATH/$CRAC_ALIGN_CMD \
		                -i $OUTPUT_INDEX_PATH/$GENOME_NAME \
		               	-k $K \
		                -r $READS_PATH/1M.$READ_1_FILE $READS_PATH/1M.$READ_2_FILE \
		                --sam $OUTPUT_ALIGNMENT_PATH/output.sam \
		                --no-ambiguity \
		                --reads-length $READ_LENGTH \
		                --min-percent-single-loc $MIN_PERCENT_SINGLE_LOC \
		                --min-percent-multiple-loc $MIN_PERCENT_MULTIPLE_LOC \
		                --min-repetition $MIN_REPETITION \
		                --min-percent-repetition-loc $MIN_PERCENT_REPETITION_LOC \
		                --summary $OUTPUT_ALIGNMENT_PATH/summary.txt \
		                --nb-threads $NUM_THREAD
		else
			# run crac (with --no-ambiguity, with --max-locs default, with --min-percent-repetition-loc default)
		        $CRAC_PATH/$CRAC_ALIGN_CMD \
		                -i $OUTPUT_INDEX_PATH/$GENOME_NAME \
		               	-k $K \
		                -r $READS_PATH/1M.$READ_1_FILE $READS_PATH/1M.$READ_2_FILE \
		                --sam $OUTPUT_ALIGNMENT_PATH/output.sam \
		                --no-ambiguity \
		                --reads-length $READ_LENGTH \
		                --min-percent-single-loc $MIN_PERCENT_SINGLE_LOC \
		                --min-percent-multiple-loc $MIN_PERCENT_MULTIPLE_LOC \
		                --min-repetition $MIN_REPETITION \
		                --summary $OUTPUT_ALIGNMENT_PATH/summary.txt \
		                --nb-threads $NUM_THREAD
		fi
	fi
fi
