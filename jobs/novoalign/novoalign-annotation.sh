#!/bin/bash
#BSUB -J novoalign-annotation_job            # LSF job name
#BSUB -o novoalign-annotation_job.%J.out     # Name of the job output file
#BSUB -e novoalign-annotation_job.%J.error   # Name of the job error file

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
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
mkdir -p $OUTPUT_INDEX_PATH
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY

# create transcriptome with USEQ
java -Xmx64G -jar $USEQ_PATH/$USEQ_MAKE_TRANSCRIPTOME_CMD \
	-f $GENOME_PATH/files \
	-u $USEQ_ANNOTATION_PATH/$USEQ_ANNOTATION_FILE \
	-r 96

mv $USEQ_ANNOTATION_PATH/${USEQ_ANNOTATION_PREFIX}Rad96*.fasta.gz $OUTPUT_INDEX_PATH/
cd $OUTPUT_INDEX_PATH/
gunzip ${USEQ_ANNOTATION_PREFIX}Rad96*.fasta.gz

# create masked exon with USEQ
java -jar $USEQ_PATH/$USEQ_MASK_EXONS_IN_FASTA_FILES_CMD \
       -f $GENOME_PATH/files \
       -u $USEQ_ANNOTATION_PATH/$USEQ_ANNOTATION_FILE \
       -s $OUTPUT_INDEX_PATH/maskedExons

cd $OUTPUT_INDEX_PATH/maskedExons
for f in chr*.fasta; do (cat "${f}"; echo) >> geneMaskedGenome.fasta; done
