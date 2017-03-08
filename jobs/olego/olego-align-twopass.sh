#!/bin/bash
#BSUB -J olego-align-twopass_job            # LSF job name
#BSUB -o olego-align-twopass_job.%J.out     # Name of the job output file
#BSUB -e olego-align-twopass_job.%J.error   # Name of the job error file

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
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY/twopass
mkdir -p $OUTPUT_ALIGNMENT_PATH


# run olego
$OLEGO_PATH/$OLEGO_ALIGN_CMD \
	--output-file $OUTPUT_ALIGNMENT_PATH/output_1.sam \
	--num-threads $NUM_THREAD \
	--regression-model $OUTPUT_ANNOTATION_PATH/$OLEGO_REGRESSION_MODEL_FILE \
	--verbose \
	--junction-file $OUTPUT_ANNOTATION_PATH/$OLEGO_JUNCTION_FILE \
	$OUTPUT_INDEX_PATH/$GENOME_NAME \
	$READS_PATH/$READ_1_FILE


$OLEGO_PATH/$OLEGO_ALIGN_CMD \
        --output-file $OUTPUT_ALIGNMENT_PATH/output_2.sam \
        --num-threads $NUM_THREAD \
        --regression-model $OUTPUT_ANNOTATION_PATH/$OLEGO_REGRESSION_MODEL_FILE \
        --verbose \
        --junction-file $OUTPUT_ANNOTATION_PATH/$OLEGO_JUNCTION_FILE \
        $OUTPUT_INDEX_PATH/$GENOME_NAME \
        $READS_PATH/$READ_2_FILE

# merge 
perl $OLEGO_PATH/mergePEsam.pl -v $OUTPUT_ALIGNMENT_PATH/output_1.sam $OUTPUT_ALIGNMENT_PATH/output_2.sam $OUTPUT_ALIGNMENT_PATH/output.sam
rm $OUTPUT_ALIGNMENT_PATH/output_1.sam
rm $OUTPUT_ALIGNMENT_PATH/output_2.sam 

# convert the SAM file to BED file
perl $OLEGO_PATH/sam2bed.pl -v --use-RNA-strand $OUTPUT_ALIGNMENT_PATH/output.sam $OUTPUT_ALIGNMENT_PATH/output.bed

# find the junctions in the BED file
perl $OLEGO_PATH/bed2junc.pl $OUTPUT_ALIGNMENT_PATH/output.bed $OUTPUT_ALIGNMENT_PATH/output.junc

# second pass
$OLEGO_PATH/$OLEGO_ALIGN_CMD \
        --output-file $OUTPUT_ALIGNMENT_PATH/output_1.twopass.sam \
        --num-threads $NUM_THREAD \
        --regression-model $OUTPUT_ANNOTATION_PATH/$OLEGO_REGRESSION_MODEL_FILE \
        --verbose \
        --junction-file $OUTPUT_ALIGNMENT_PATH/output.junc \
	--non-denovo \
	$OUTPUT_INDEX_PATH/$GENOME_NAME \
        $READS_PATH/$READ_1_FILE \

$OLEGO_PATH/$OLEGO_ALIGN_CMD \
        --output-file $OUTPUT_ALIGNMENT_PATH/output_2.twopass.sam \
        --num-threads $NUM_THREAD \
        --regression-model $OUTPUT_ANNOTATION_PATH/$OLEGO_REGRESSION_MODEL_FILE \
        --verbose \
        --junction-file $OUTPUT_ALIGNMENT_PATH/output.junc \
	--non-denovo \
	$OUTPUT_INDEX_PATH/$GENOME_NAME \
        $READS_PATH/$READ_2_FILE

# merge
perl $OLEGO_PATH/mergePEsam.pl -v $OUTPUT_ALIGNMENT_PATH/output_1.twopass.sam $OUTPUT_ALIGNMENT_PATH/output_2.twopass.sam $OUTPUT_ALIGNMENT_PATH/output.twopass.sam
rm $OUTPUT_ALIGNMENT_PATH/output.sam



# --non-denovo 		No de novo junction search
# --strand-mode		<value>
# --junction-file <file>

