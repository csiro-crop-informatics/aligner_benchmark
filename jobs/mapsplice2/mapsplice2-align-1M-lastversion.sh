#!/bin/bash
#BSUB -J mapsplice2-align_job            # LSF job name
#BSUB -o mapsplice2-align_job.%J.out     # Name of the job output file
#BSUB -e mapsplice2-align_job.%J.error   # Name of the job error file

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="mapsplice2"

# define library
LIBRARY=$DATASET


# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index-lastversion/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-1M-lastversion
mkdir -p $OUTPUT_ALIGNMENT_PATH


# run mapsplice2
cd $LASTVERSION_MAPSPLICE2_PATH/mapping
$LASTVERSION_MAPSPLICE2_PATH/mapping/mps \
	-1 $READS_PATH/1M.$READ_1_FILE_CHANGE_ID \
	-2 $READS_PATH/1M.$READ_2_FILE_CHANGE_ID \
	-G $OUTPUT_INDEX_PATH \
	-J \
	-T $NUM_THREAD \
	-O $OUTPUT_ALIGNMENT_PATH
	

#Processing Paired-end reads:
#command-line main arguments:
#-1 <string> path_to_read_file_end1
#-2 <string> path_to_read_file_end2
#-G <string> path_to_global_index
#-L <string> path_to_local_index
#-T <int> threads_num
#-O <string> path_to_output_folder
#optional arguments:
#-A Do_phase1_only
#-I annotation file
#-F <string> input_format_fasta_or_fastq
#-H command line help information


#python $MAPSPLICE2_PATH/$MAPSPLICE2_ALIGN_CMD \
#	--threads $NUM_THREAD \
#	--non-canonical-double-anchor \
#	--output $OUTPUT_ALIGNMENT_PATH \
#	-c $GENOME_PATH/files_nospace \
#	-x $OUTPUT_INDEX_PATH/$GENOME_NAME \
#	-1 $READS_PATH/1M.$READ_1_FILE_CHANGE_ID \
#	-2 $READS_PATH/1M.$READ_2_FILE_CHANGE_ID

