# Script to convert the BEERS .BED file into a standard .BED file.
# Usage:
# fromBeersBedtoStandardBed.py -i <_beers_bed_file> -o <standard_bed_file>'
import string, sys, getopt

def main(argv):

	inputfile = ''
	outputfile = ''
	try:
	  opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
	except getopt.GetoptError:
	  print 'fromBeersBedtoStandardBed.py -i <inputfile> -o <outputfile>'
	  sys.exit(2)
	for opt, arg in opts:
	  if opt == '-h':
	     print 'fromBeersBedtoStandardBed.py -i <inputfile> -o <outputfile>'
	     sys.exit()
	  elif opt in ("-i", "--ifile"):
	     inputfile = arg
	  elif opt in ("-o", "--ofile"):
	     outputfile = arg
	print 'Input file is "', inputfile, '"'
	print 'Output file is "', outputfile, '"'


if __name__ == "__main__":
	main(sys.argv[1:])

	# BEERS Bed file fields
	BEERS_BED_FIELD_CHR=0
	BEERS_BED_FIELD_STRAND=1
	BEERS_BED_FIELD_CHR_START=2
	BEERS_BED_FIELD_CHR_END=3
	BEERS_BED_FIELD_BLOCK_COUNT=4
	BEERS_BED_FIELD_BLOCK_STARTS=5
	BEERS_BED_FIELD_BLOCK_ENDS=6
	BEERS_BED_FIELD_FEATURE_NAME=7

	# Standard BED file fields
	STD_BED_FIELD_CHR=0
	STD_BED_FIELD_CHR_START=1
	STD_BED_FIELD_CHR_END=2
	STD_BED_FIELD_FEATURE_NAME=3
	STD_BED_FIELD_SCORE=4
	STD_BED_FIELD_STRAND=5
	STD_BED_FIELD_THICK_START=6
	STD_BED_FIELD_THICK_END=7
	STD_BED_FIELD_ITEM_RGB=8
	STD_BED_FIELD_BLOCK_COUNT=9
	STD_BED_FIELD_BLOCK_SIZES=10
	STD_BED_FIELD_BLOCK_STARTS=11


	default_item_rgb_value = "255,0,0"
	default_score_value = "0"

	# open input file
	input_file = open(sys.argv[2], 'r')

	# create output file
	output_file = open(sys.argv[4], 'w')


	for line in input_file:

		# read BEERS BED line
		beers_bed_entry = line.split()

		# read block end positions
		ends = string.split(beers_bed_entry[BEERS_BED_FIELD_BLOCK_ENDS], ",")
		ends.pop()
		ends = map(int, ends)

		# read block start positions
		starts = string.split(beers_bed_entry[BEERS_BED_FIELD_BLOCK_STARTS], ",")
		starts.pop()
		starts = map(int, starts)

		# find "chromStart" and "chromEnd" for standard BED file
		chr_start = min(starts)
		chr_end = max(ends)

		# BED specification:
		# all "blockStart" positions must be calculated relative to "chromStart"
		# the first "blockStart" value must be 0, so that the first block begins at "chromStart". 
		# the final "blockStart" position plus the final "blockSize" value must equal "chromEnd"
		# blocks may not overlap.

		# remove "chromStart" offset 
		starts[:] = [start - chr_start for start in starts]
		ends[:] = [end - chr_start for end in ends]

		# calculate block sizes
		sizes = [end-start for end,start in zip(ends,starts)]

		# convert to string
		sizes = map(str, sizes)
		starts = map(str, starts)	


		# prepare standard BED entry
		std_bed_entry = list()

		std_bed_entry.insert (STD_BED_FIELD_CHR, 			beers_bed_entry[BEERS_BED_FIELD_CHR])
		std_bed_entry.insert (STD_BED_FIELD_CHR_START, 		str(chr_start))
		std_bed_entry.insert (STD_BED_FIELD_CHR_END, 		str(chr_end))
		std_bed_entry.insert (STD_BED_FIELD_FEATURE_NAME,	beers_bed_entry[BEERS_BED_FIELD_FEATURE_NAME])
		std_bed_entry.insert (STD_BED_FIELD_SCORE,			default_score_value)
		std_bed_entry.insert (STD_BED_FIELD_STRAND,			beers_bed_entry[BEERS_BED_FIELD_STRAND])
		std_bed_entry.insert (STD_BED_FIELD_THICK_START,	str(chr_start))
		std_bed_entry.insert (STD_BED_FIELD_THICK_END,		str(chr_end))
		std_bed_entry.insert (STD_BED_FIELD_ITEM_RGB,		default_item_rgb_value)
		std_bed_entry.insert (STD_BED_FIELD_BLOCK_COUNT,	beers_bed_entry[BEERS_BED_FIELD_BLOCK_COUNT])
		std_bed_entry.insert (STD_BED_FIELD_BLOCK_SIZES, 	','.join(sizes)+",")
		std_bed_entry.insert (STD_BED_FIELD_BLOCK_STARTS,	','.join(starts))

		tab = "\t"
		
		# write standard BED entry
		output_file.write (tab.join(std_bed_entry)+"\n")
		
