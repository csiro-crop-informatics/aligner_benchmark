#!/bin/bash

# script to remove adapter sequence using CutAdapt
# Usage: 
# remove_adapter.sh read_forward.fa read_reverse.fa

source $HOME/my_python-2.7.5/bin/activate

source /project/itmatlab/aligner_benchmark/jobs/settings/variable.sh

INPUT_FILE_1=$1
INPUT_FILE_2=$2

OUTPUT_FILE_1=trimmed.$1
OUTPUT_FILE_2=trimmed.$2

$CUTADAPT_PATH/$CUTADAPT_CMD -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o $OUTPUT_FILE_1 -p $OUTPUT_FILE_2 $INPUT_FILE_1 $INPUT_FILE_2


# sed 's/b//' input_2.fa  > input_2_new.fa
# sed 's/a//' input_1.fa  > input_1_new.fa
# sed 's/>seq.[0-9]*/&b/' output_2.fa > output_2_new.fa
# sed 's/>seq.[0-9]*/&a/' output_1.fa > output_1_new.fa
