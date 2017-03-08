#!/bin/bash
#BSUB -J olego-annotation_job            # LSF job name
#BSUB -o olego-annotation_job.%J.out     # Name of the job output file
#BSUB -e olego-annotation_job.%J.error   # Name of the job error file

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="olego"

# define library
LIBRARY=$DATASET


# define output path
OUTPUT_ANNOTATION_PATH=$OUTPUT_ROOT_PATH/$TOOL/annotation
mkdir -p $OUTPUT_ANNOTATION_PATH
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY


# create regression model
echo "create regression model"
cd $OUTPUT_ANNOTATION_PATH
perl $OLEGO_PATH/regression_model_gen/OLego_regression.pl -g $GENOME_PATH/files -a $BED_PATH/$BED_FILE -o $OLEGO_REGRESSION_MODEL_FILE_PREFIX
cp $OUTPUT_ANNOTATION_PATH/$OLEGO_REGRESSION_MODEL_FILE_PREFIX.cache/$OLEGO_REGRESSION_MODEL_FILE $OUTPUT_ANNOTATION_PATH


# create junctions file
echo "create junctions file"
cd $OUTPUT_ANNOTATION_PATH
perl $OLEGO_PATH/bed2junc.pl $BED_PATH/$BED_FILE $OLEGO_JUNCTION_FILE 


# run olego
#$OLEGO_PATH/$OLEGO_ALIGN_CMD \
#	--output-file $OUTPUT_ALIGNMENT_PATH/output_1.sam \
#	--num-threads $NUM_THREAD \
#	--regression-model $OLEGO_PATH/models/hg.cfg \
#	--verbose \
#	--junction-file $OUTPUT_ROOT_PATH/$TOOL/annotation/$OLEGO_JUNCTION_FILE \
#	$OUTPUT_INDEX_PATH/$GENOME_NAME \
#	$READS_PATH/$READ_1_FILE \


#$OLEGO_PATH/$OLEGO_ALIGN_CMD \
#        --output-file $OUTPUT_ALIGNMENT_PATH/output_2.sam \
#        --num-threads $NUM_THREAD \
#        --regression-model $OLEGO_PATH/models/hg.cfg \
#        --verbose \
#        --junction-file $OUTPUT_ROOT_PATH/$TOOL/annotation/$OLEGO_JUNCTION_FILE	\
#        $OUTPUT_INDEX_PATH/$GENOME_NAME \
#        $READS_PATH/$READ_2_FILE

#perl $OLEGO_PATH/mergePEsam.pl -v $OUTPUT_ALIGNMENT_PATH/output_1.sam $OUTPUT_ALIGNMENT_PATH/output_2.sam $OUTPUT_ALIGNMENT_PATH/output.sam

# --non-denovo 		No de novo junction search
# --strand-mode		<value>
# --junction-file <file>

