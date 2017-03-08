#!/bin/bash
#BSUB -J gsnap-align_job            # LSF job name
#BSUB -o gsnap-align_job.%J.out     # Name of the job output file
#BSUB -e gsnap-align_job.%J.error   # Name of the job error file

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


# define output path
OUTPUT_INDEX_PATH=$OUTPUT_ROOT_PATH/$TOOL/index-lastversion/$GENOME_NAME
OUTPUT_ALIGNMENT_PATH=$OUTPUT_ROOT_PATH/$TOOL/alignment/$LIBRARY-1M-lastversion
mkdir -p $OUTPUT_ALIGNMENT_PATH


# run gsnap
$LASTVERSION_GSNAP_PATH/$GSNAP_ALIGN_CMD \
        -D $OUTPUT_INDEX_PATH \
        -d $GENOME_NAME \
        -A sam \
        --merge-distant-samechr \
	--ordered \
	--novelsplicing 1 \
        --use-splicing $GENOME_NAME.splicesites \
        --nthreads $NUM_THREAD \
        --batch 5 \
        --expand-offsets 1 \
	$READS_PATH/1M.$READ_1_FILE $READS_PATH/1M.$READ_2_FILE \
        > $OUTPUT_ALIGNMENT_PATH/output.sam


#--#--#   TODO: rimettere "--add-paired-nomappers"  #--#--#





#  -q, --part=INT/INT             Process only the i-th out of every n sequences e.g., 0/100 or 99/100 (useful for distributing jobs to a computer farm).
#  --orientation=STRING           Orientation of paired-end reads.Allowed values: FR (fwd-rev, or typical Illumina; default),RF (rev-fwd, for circularized inserts), or FF (fwd-fwd, same strand)
#  --novelsplicing			detect splices using a probabilistic model	
#  --use-splicing=STRING (or -s)  Look for splicing involving known sites or known introns (in <STRING>.iit), at short or long distances. See README instructions for the distinction between known sites and known introns
#  --clip-overlap                 For paired-end reads whose alignments overlap, clip the overlapping region.
#  --merge-overlap                For paired-end reads whose alignments overlap, merge the two ends into a single end (beta implementation)
#  --merge-distant-samechr	  Report distant splices on the same chromosome as a single splice, if possible. Will produce a single SAM line instead of two SAM lines, which is also done for translocations, inversions, and scramble events

# --expand-offsets=INT           Whether to expand the genomic offsets index. Values: 0 (no, default), or 1 (yes). Expansion gives faster alignment, but requires more memory
# --batch=INT			From "0"(low performance, low memory) to "5" (high performance, high memory)
# --nthreads=INT             Number of worker threads

#  --pairmax-rna=INT              Max total genomic length for RNA-Seq paired reads, or other reads that could have a splice (default 200000).  Used if -N or -s is specified. Should probably match the value for -w, --localsplicedist.
#  --pairexpect=INT               Expected paired-end length, used for calling splices in medial part of paired-end reads (default 200).  Was turned off in previous versions, but reinstated. 
#  --pairdev=INT                  Allowable deviation from expected paired-end length, used for calling splices in medial part of paired-end reads (default 100). Was turned off in previous versions, but reinstated.

#  -O, --ordered                  Print output in same order as input (relevant only if there is more than one worker thread)
