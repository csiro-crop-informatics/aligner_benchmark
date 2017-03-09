#!/bin/bash
#BSUB -J contextmap2-align_job            # LSF job name
#BSUB -o contextmap2-align_job.%J.out     # Name of the job output file 
#BSUB -e contextmap2-align_job.%J.error   # Name of the job error file

if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh
fi
module load java-sdk-1.7.0

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="contextmap2"

# define library
LIBRARY=$DATASET


# define output	path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY
mkdir -p $OUTPUT_ALIGNMENT_PATH

# ALIGNER must be one of "bwa", "bowtie1" or "bowtie2"
ALIGNER="bwa"
ALIGNER_BIN_PATH=""
INDEXER_BIN_PATH=""


# set index variable
echo "build $ALIGNER index"
if [ "${ALIGNER}" = bowtie1 ]
then
	ALIGNER_BIN_PATH=$BOWTIE_PATH/$BOWTIE_ALIGN_CMD
	INDEXER_BIN_PATH=$BOWTIE_PATH/$BOWTIE_INDEX_CMD
else
	if [ "${ALIGNER}" = bowtie2 ]
	then 
		ALIGNER_BIN_PATH=$BOWTIE2_PATH/$BOWTIE2_ALIGN_CMD
		INDEXER_BIN_PATH=$BOWTIE2_PATH/$BOWTIE2_INDEX_CMD
	else
		if [ "${ALIGNER}" = bwa ]
		then 
			ALIGNER_BIN_PATH=$BWA_PATH/$BWA_ALIGN_CMD
			INDEXER_BIN_PATH=$BWA_PATH/$BWA_INDEX_CMD
		fi
	fi
fi


# run contextmap2 with multiple thread
# performance parameters set with the values suggested in the manual, except for the max memory(-Xmx): suggested 12GB, here 64GB)
java -Xms4000M -Xmx64000m -XX:+UseConcMarkSweepGC -XX:NewSize=300M -XX:MaxNewSize=300M -jar $CONTEXTMAP2_PATH/$CONTEXTMAP2_ALIGN_CMD \
	-reads $READS_PATH/$READS_FILE_CHANGE_ID \
	--pairedend \
	-gtf $GTF_PATH/$GTF_FILE \
	--noncanonicaljunctions \
	-aligner_name $ALIGNER \
	-aligner_bin $ALIGNER_BIN_PATH \
	-indexer_bin $INDEXER_BIN_PATH \
	-indices $OUTPUT_INDEX_PATH/$GENOME_NAME \
	-genome $GENOME_PATH/files \
	-o $OUTPUT_ALIGNMENT_PATH \
	--verbose \
	-t $NUM_THREAD

#-gtf $GTF_PATH/$GTF_FILE (not mandatory)
#-t
#--noncanonicaljunctions
#--strandspecific (works only with single end data)


#java -Xms1000m -Xmx4000m -jar $CONTEXTMAP2_PATH/$CONTEXTMAP2_ALIGN_CMD \
#	-reads $READS_PATH/$READS_FILE \
#	--pairedend \
#	-gtf $GTF_PATH/$GTF_FILE \
#	-aligner_name $ALIGNER \
#	-aligner_bin $ALIGNER_BIN_PATH \
#	-indexer_bin $INDEXER_BIN_PATH \
#	-indices $OUTPUT_INDEX_PATH/$GENOME_NAME \
#	-genome $GENOME_PATH \
#	-o $OUTPUT_ALIGNMENT_PATH \
#	--verbose

#java -Xms1000m -Xmx1000m -jar $contextmap_jar mapper -aligner_name $aligner_name -aligner_bin $aligner_bin -indexer_bin $indexer_bin -reads $reads -o $output_dir -indices $aligner_index -genome $genomes_dir >> log.txt 2>&1

#java -Xms1000m -Xmx4000m -XX:+UseConcMarkSweepGC -XX:NewSize=300M -XX:MaxNewSize=300M -jar $contextmap_jar mapper -reads $reads -readlen $readlen -o $output_dir -rmap $rmap -bwtbuild $bwtbuild -bwtindex $bwtindex -genome $genome_dir -t 8

