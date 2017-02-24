# Script to convert a .BED file into a RefFlat format.
# Usage:
# fromBedtoRefFlat.py -i <bed_file> -o <ref_flat_file>'
import string, sys, getopt

def main(argv):

	inputfile = ''
	outputfile = ''
	try:
	  opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
	except getopt.GetoptError:
	  print 'fromBedtoRefFlat.py -i <inputfile> -o <outputfile>'
	  sys.exit(2)
	for opt, arg in opts:
	  if opt == '-h':
	     print 'fromBedtoRefFlat.py -i <inputfile> -o <outputfile>'
	     sys.exit()
	  elif opt in ("-i", "--ifile"):
	     inputfile = arg
	  elif opt in ("-o", "--ofile"):
	     outputfile = arg
	print 'Input file is "', inputfile, '"'
	print 'Output file is "', outputfile, '"'


if __name__ == "__main__":
	main(sys.argv[1:])

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


	# refFlat file fields
	REFFLAT_FIELD_GENE_NAME=0
	REFFLAT_FIELD_GENE=1
	REFFLAT_FIELD_CHR=2
	REFFLAT_FIELD_STRAND=3
	REFFLAT_FIELD_TRANSCRIPTION_START=4
	REFFLAT_FIELD_TRANSCRIPTION_END=5
	REFFLAT_FIELD_CDS_START=6
	REFFLAT_FIELD_CDS_END=7
	REFFLAT_FIELD_EXON_COUNT=8
	REFFLAT_FIELD_EXON_STARTS=9
	REFFLAT_FIELD_EXON_ENDS=10

	default_item_rgb_value = "255,0,0"
	default_score_value = "0"

	# open input file
	input_file = open(sys.argv[2], 'r')

	# create output file
	output_file = open(sys.argv[4], 'w')


	for line in input_file:

		# read BED line
		bed_entry = line.split()

		# read chromosome start position
		chr_start = int(bed_entry[STD_BED_FIELD_CHR_START])

		# read block sizes
		sizes = string.split(bed_entry[STD_BED_FIELD_BLOCK_SIZES], ",")
		sizes.pop()
		sizes = map(int, sizes)

		# read block start positions
		starts = string.split(bed_entry[STD_BED_FIELD_BLOCK_STARTS], ",")
		starts = map(int, starts)
		starts = [start+chr_start for start in starts]

		

		# calculate block ends
		ends = [start+size for start,size in zip(starts,sizes)]


		# find starts and ends of transcription and CDS
		cds_start = min(starts)
		cds_end = max(ends)
		transcription_start = cds_start
		transcription_end = cds_end


		# convert to string
		starts = map(str, starts)
		ends = map(str, ends)



		# prepare refFlat entry
		refFlat_entry = list()

		refFlat_entry.insert (REFFLAT_FIELD_GENE_NAME, 				bed_entry[STD_BED_FIELD_FEATURE_NAME])
		refFlat_entry.insert (REFFLAT_FIELD_GENE, 					bed_entry[STD_BED_FIELD_FEATURE_NAME])
		refFlat_entry.insert (REFFLAT_FIELD_CHR, 					bed_entry[STD_BED_FIELD_CHR])
		refFlat_entry.insert (REFFLAT_FIELD_STRAND,					bed_entry[STD_BED_FIELD_STRAND])
		refFlat_entry.insert (REFFLAT_FIELD_TRANSCRIPTION_START,	str(transcription_start))
		refFlat_entry.insert (REFFLAT_FIELD_TRANSCRIPTION_END,		str(transcription_end))
		refFlat_entry.insert (REFFLAT_FIELD_CDS_START,				str(cds_start))
		refFlat_entry.insert (REFFLAT_FIELD_CDS_END,				str(cds_end))
		refFlat_entry.insert (REFFLAT_FIELD_EXON_COUNT,				bed_entry[STD_BED_FIELD_BLOCK_COUNT])
		refFlat_entry.insert (REFFLAT_FIELD_EXON_STARTS,			','.join(starts)+",")
		refFlat_entry.insert (REFFLAT_FIELD_EXON_ENDS, 				','.join(ends)+",")

		tab = "\t"
		
		# write standard BED entry
		output_file.write (tab.join(refFlat_entry)+"\n")
		
