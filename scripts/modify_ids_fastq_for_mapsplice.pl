# Script to change the read id 
# from <seq_id>.a into <seq_id>.1 
# from <seq_id>.b into <seq_id>.2

# Usage:
# modify_ids_fastq_for_mapsplice.pl <read_file> > <new_read_file>

open(INFILE, $ARGV[0]);
while($line = <INFILE>) {
    $line =~ s!a!/1!;
    $line =~ s!b!/2!;
    print $line;
}
close(INFILE);

