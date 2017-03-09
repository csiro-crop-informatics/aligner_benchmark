#!/bin/bash
#BSUB -J novoalign-align-tune_job            # LSF job name
#BSUB -o novoalign-align-tune_job.%J.out     # Name of the job output file
#BSUB -e novoalign-align-tune_job.%J.error   # Name of the job error file

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

# set parameters to tune
REPEAT_FILTER=$3
PAIREND_DEFAULT=$4
A_SCORE=$5
B_SCORE=$6

# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-tune/tune-$REPEAT_FILTER-$PAIREND_DEFAULT-$A_SCORE-$B_SCORE
mkdir -p $OUTPUT_ALIGNMENT_PATH


if [ $REPEAT_FILTER != "off" ]
then
	if [ $PAIREND_DEFAULT != "on" ]
	then
		# run novoalign
		$NOVOALIGN_PATH/$NOVOALIGN_ALIGN_CMD \
        		-d $OUTPUT_INDEX_PATH/$GENOME_NAME.nix \
	        	-f $READS_PATH/1M.$READ_1_FILE $READS_PATH/1M.$READ_2_FILE \
		        -F FA \
		        -o SAM \
	        	-r All 10 \
			-t $A_SCORE,$B_SCORE \
		        -i PE $FRAGMENT_LENGTH_MEAN,$FRAGMENT_LENGTH_SD \
		        -v 0 70 70 "[>]([^:]*)" \
		        > $OUTPUT_ALIGNMENT_PATH/output.tmp.sam \
	        	2> $OUTPUT_ALIGNMENT_PATH/alignment.log
	else
                # run novoalign
                $NOVOALIGN_PATH/$NOVOALIGN_ALIGN_CMD \
                        -d $OUTPUT_INDEX_PATH/$GENOME_NAME.nix \
                        -f $READS_PATH/1M.$READ_1_FILE $READS_PATH/1M.$READ_2_FILE \
                       	-F FA \
                       	-o SAM \
                        -r All 10 \
                       	-t $A_SCORE,$B_SCORE \
                        -v 0 70 70 "[>]([^:]*)" \
                        > $OUTPUT_ALIGNMENT_PATH/output.tmp.sam \
                        2> $OUTPUT_ALIGNMENT_PATH/alignment.log
	fi
else
    	if [ $PAIREND_DEFAULT != "on" ]
        then
		# run novoalign
		$NOVOALIGN_PATH/$NOVOALIGN_ALIGN_CMD \
			-d $OUTPUT_INDEX_PATH/$GENOME_NAME.nix \
			-f $READS_PATH/1M.$READ_1_FILE $READS_PATH/1M.$READ_2_FILE \
			-F FA \
		        -o SAM \
		        -r All 10 \
                       	-t $A_SCORE,$B_SCORE \
			-h -1 -1 \
		        -i PE $FRAGMENT_LENGTH_MEAN,$FRAGMENT_LENGTH_SD \
		        -v 0 70 70 "[>]([^:]*)" \
		        > $OUTPUT_ALIGNMENT_PATH/output.tmp.sam \
		        2> $OUTPUT_ALIGNMENT_PATH/alignment.log
	else
                $NOVOALIGN_PATH/$NOVOALIGN_ALIGN_CMD \
                        -d $OUTPUT_INDEX_PATH/$GENOME_NAME.nix \
                        -f $READS_PATH/1M.$READ_1_FILE $READS_PATH/1M.$READ_2_FILE \
                       	-F FA \
                       	-o SAM \
                        -r All 10 \
                       	-t $A_SCORE,$B_SCORE \
                        -h -1 -1 \
                        -v 0 70 70 "[>]([^:]*)" \
                        > $OUTPUT_ALIGNMENT_PATH/output.tmp.sam \
                        2> $OUTPUT_ALIGNMENT_PATH/alignment.log	
	fi
fi

echo "alignment completete... i will fix qualities"

$SCRIPTS_PATH/fix-missing-qualities.sh $OUTPUT_ALIGNMENT_PATH

java -Xmx64G -jar $USEQ_PATH/$USEQ_SAM_TRANSCRIPTOME_PARSER_CMD \
	-f $OUTPUT_ALIGNMENT_PATH/output.tmp.fix.sam \
	-a 50000 \
	-n 100 \
	-u \
	-s $OUTPUT_ALIGNMENT_PATH/output.sam
