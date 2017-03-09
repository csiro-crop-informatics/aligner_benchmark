#!/bin/bash
#BSUB -J hisat-align_job            # LSF job name
#BSUB -o hisat-align_job.%J.out     # Name of the job output file
#BSUB -e hisat-align_job.%J.error   # Name of the job error file

# load variables
NUM_THREAD=$1

DATASET_PARAM=$2
source $DATASET_PARAM

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

# define tool
TOOL="hisat"

# define library
LIBRARY=$DATASET


# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY
mkdir -p $OUTPUT_ALIGNMENT_PATH

# run hisat
$HISAT_PATH/$HISAT_ALIGN_CMD \
	--threads $NUM_THREAD \
	--time \
	--reorder \
	--known-splicesite-infile $OUTPUT_INDEX_PATH/$GENOME_NAME.splicesites.txt \
	--novel-splicesite-outfile $OUTPUT_ALIGNMENT_PATH/$GENOME_NAME.splicesites.novel.txt \
	--novel-splicesite-infile $OUTPUT_ALIGNMENT_PATH/$GENOME_NAME.splicesites.novel.txt \
	-f \
	-x $OUTPUT_INDEX_PATH/$GENOME_NAME \
	-1 $READS_PATH/$READ_1_FILE \
	-2 $READS_PATH/$READ_2_FILE \
	-S $OUTPUT_ALIGNMENT_PATH/output.sam \
	

# -p/--threads <num_thread>
# --rna-strandness <string>	single end: F or R ; paired end: RF or FR ; default: unstranded
# --fr/--rf/--ff 		The upstream/downstream mate orientations for a valid paired-end alignment against the forward reference strand.
# -t/--time			Print the wall-clock time required to load the index files and align the reads. This is printed to the "standard error" ("stderr") filehandle. Default: off.
# --reorder			Guarantees that output SAM records are printed in an order corresponding to the order of the reads in the original input file, even when -p is set greater than 1.
# --mm				Use memory-mapped I/O to load the index, rather than typical file I/O. (useful??)

#--known-splicesite-infile <path>	With this mode, you can provide a list of known splice sites, which HISAT makes use of to align reads with small anchors. You can create such a list using python extract_splice_sites.py genes.gtf > splicesites.txt, where extract_splice_sites.py is included in the HISAT package, genes.gtf is a gene annotation file, and splicesites.txt is a list of splice sites with which you provide HISAT in this mode.
#--novel-splicesite-outfile <path>	In this mode, HISAT reports a list of splice sites in the file : chromosome name <tab> genomic position of the flanking base on the left side of an intron <tab> genomic position of the flanking base on the right <tab> strand
#--novel-splicesite-infile <path>	With this mode, you can provide a list of novel splice sites that were generated from the above option "--novel-splicesite-outfile".

