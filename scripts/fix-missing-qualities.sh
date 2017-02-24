#!/bin/sh

# Script to add a quality string to Novoalign output SAM.
# This is necessary for the following postprocessing with USEQ 
# (which otherwise recognize a quality string length different from the sequence length)

if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh
fi
module load samtools-1.1

cd $1

samtools view -H output.tmp.sam > header.sam
samtools view output.tmp.sam > align.txt

num_line=$(wc -l align.txt | cut -d " " -f 1 -)

cut -f 1-10 align.txt > left.part
cut -f 12- align.txt > right.part

num_line_per_file=100000
num_file=$(( num_line / num_line_per_file ))
num_iter=$(( num_file - 1 ))
num_line_remain=$(( num_line % num_line_per_file ))

for i in $(seq 1 ${num_line_per_file}) 
do
	echo "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" >> lines.tmp.txt 
done

touch lines.txt

for i in $(seq 1 ${num_file})
do
	cat lines.txt lines.tmp.txt > tmp.txt
	mv tmp.txt lines.txt
done

for i in $(seq 1 ${num_line_remain}) 
do 
	echo "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" >> lines.txt
done

paste left.part lines.txt right.part > align.fix.sam

cat header.sam align.fix.sam > output.tmp.fix.sam

rm left.part lines.txt right.part align.fix.sam align.txt header.sam lines.tmp.txt 
