# Script to add adapter tails on read sequences.
# The script will add adapter sequences of length 10,20,30,40,50 bp on the first 500000 read pairs (5 chunks of 100000 read pairs)
# The Illumina Universal Adapter as template sequence
# Usage:
# add_adapter2fasta_V2.pl read_forward.fa read_reverse.fa

# Illumina Universal Adapter
$adapter = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT";
@ADAPTER = split(//,$adapter);

# reverse complement of the adapter sequence (for reverse reads)
$adapter_rev = reversecomplement($adapter);
@ADAPTER_REV = split(//,$adapter_rev);

$f10 = "";
$f20 = "";
$f30 = "";
$f40 = "";
$f50 = "";
for($i=0; $i<10; $i++) {
    $f10 = $f10 . $ADAPTER[$i];
}
for($i=0; $i<20; $i++) {
    $f20 = $f20 . $ADAPTER[$i];
}
for($i=0; $i<30; $i++) {
    $f30 = $f30 . $ADAPTER[$i];
}
for($i=0; $i<40; $i++) {
    $f40 = $f40 . $ADAPTER[$i];
}
for($i=0; $i<50; $i++) {
    $f50 = $f50 . $ADAPTER[$i];
}
$r10 = "";
$r20 = "";
$r30 = "";
$r40 = "";
$r50 = "";
for($i=0; $i<10; $i++) {
    $r10 = $r10 . $ADAPTER_REV[$i];
}
for($i=0; $i<20; $i++) {
    $r20 = $r20 . $ADAPTER_REV[$i];
}
for($i=0; $i<30; $i++) {
    $r30 = $r30 . $ADAPTER_REV[$i];
}
for($i=0; $i<40; $i++) {
    $r40 = $r40 . $ADAPTER_REV[$i];
}
for($i=0; $i<50; $i++) {
    $r50 = $r50 . $ADAPTER_REV[$i];
}

$forward_out = $ARGV[2];
$reverse_out = $ARGV[3];
open(FORWARD, ">$forward_out");

open(INFILE, $ARGV[0]); # the forward reads fa file
$flag = 0;
$count = 0;
while($line = <INFILE>) {
    if($flag == 0) {
	print FORWARD $line;
	$flag = 1;
	next;
    } 
    chomp($line);
    $flag = 0;
    if($count < 100000) {
	$line =~ s/.{10}$/$f10/;
	$r = int(rand(200));
	if($r < 10) {
	    @a = split(//,$line);
	    $x = getrandbase();
	    while($x eq $a[99-$r]) {
		$x = getrandbase();
	    }
	    $a[99-$r] = $x;
	    $line = "";
	    for($i=0; $i<100; $i++) {
		$line = $line . $a[$i];
	    }
	}
    } elsif($count < 200000) {
	$line =~ s/.{20}$/$f20/;
	$r = int(rand(200));
	if($r < 20) {
	    @a = split(//,$line);
	    $x = getrandbase();
	    while($x eq $a[99-$r]) {
		$x = getrandbase();
	    }
	    $a[99-$r] = $x;
	    $line = "";
	    for($i=0; $i<100; $i++) {
		$line = $line . $a[$i];
	    }
	}
    } elsif($count < 300000) {
	$line =~ s/.{30}$/$f30/;
	$r = int(rand(200));
	if($r < 30) {
	    @a = split(//,$line);
	    $x = getrandbase();
	    while($x eq $a[99-$r]) {
		$x = getrandbase();
	    }
	    $a[99-$r] = $x;
	    $line = "";
	    for($i=0; $i<100; $i++) {
		$line = $line . $a[$i];
	    }
	}
    } elsif($count < 400000) {
	$line =~ s/.{40}$/$f40/;
	$r = int(rand(200));
	if($r < 40) {
	    @a = split(//,$line);
	    $x = getrandbase();
	    while($x eq $a[99-$r]) {
		$x = getrandbase();
	    }
	    $a[99-$r] = $x;
	    $line = "";
	    for($i=0; $i<100; $i++) {
		$line = $line . $a[$i];
	    }
	}
    } elsif($count < 500000) {
	$line =~ s/.{50}$/$f50/;
	$r = int(rand(200));
	if($r < 50) {
	    @a = split(//,$line);
	    $x = getrandbase();
	    while($x eq $a[99-$r]) {
		$x = getrandbase();
	    }
	    $a[99-$r] = $x;
	    $line = "";
	    for($i=0; $i<100; $i++) {
		$line = $line . $a[$i];
	    }
	}
    }
    print FORWARD "$line\n";
    $count++;
}
close(INFILE);
close(FORWARD);

open(REVERSE, ">$reverse_out");
open(INFILE, $ARGV[1]); # the reverse reads fa file
$flag = 0;
$count = 0;
while($line = <INFILE>) {
    if($flag == 0) {
	print REVERSE $line;
	$flag = 1;
	next;
    } 
    chomp($line);
    $flag = 0;
    if($count < 100000) {
	$line =~ s/.{10}$/$r10/;
	$r = int(rand(200));
	if($r < 10) {
	    @a = split(//,$line);
	    $x = getrandbase();
	    while($x eq $a[99-$r]) {
		$x = getrandbase();
	    }
	    $a[99-$r] = $x;
	    $line = "";
	    for($i=0; $i<100; $i++) {
		$line = $line . $a[$i];
	    }
	}
    } elsif($count < 200000) {
	$line =~ s/.{20}$/$r20/;
	$r = int(rand(200));
	if($r < 20) {
	    @a = split(//,$line);
	    $x = getrandbase();
	    while($x eq $a[99-$r]) {
		$x = getrandbase();
	    }
	    $a[99-$r] = $x;
	    $line = "";
	    for($i=0; $i<100; $i++) {
		$line = $line . $a[$i];
	    }
	}
    } elsif($count < 300000) {
	$line =~ s/.{30}$/$r30/;
	$r = int(rand(200));
	if($r < 30) {
	    @a = split(//,$line);
	    $x = getrandbase();
	    while($x eq $a[99-$r]) {
		$x = getrandbase();
	    }
	    $a[99-$r] = $x;
	    $line = "";
	    for($i=0; $i<100; $i++) {
		$line = $line . $a[$i];
	    }
	}
    } elsif($count < 400000) {
	$line =~ s/.{40}$/$r40/;
	$r = int(rand(200));
	if($r < 40) {
	    @a = split(//,$line);
	    $x = getrandbase();
	    while($x eq $a[99-$r]) {
		$x = getrandbase();
	    }
	    $a[99-$r] = $x;
	    $line = "";
	    for($i=0; $i<100; $i++) {
		$line = $line . $a[$i];
	    }
	}
    } elsif($count < 500000) {
	$line =~ s/.{50}$/$r50/;
	$r = int(rand(200));
	if($r < 50) {
	    @a = split(//,$line);
	    $x = getrandbase();
	    while($x eq $a[99-$r]) {
		$x = getrandbase();
	    }
	    $a[99-$r] = $x;
	    $line = "";
	    for($i=0; $i<100; $i++) {
		$line = $line . $a[$i];
	    }
	}
    }
    print REVERSE "$line\n";
    $count++;
}


sub getrandbase () {
    $x = int(rand(4));
    if($x == 0) {
        return "A";
    }
    if($x == 1) {
        return "C";
    }
    if($x == 2) {
        return "G";
    }
    if($x == 3) {
        return "T";
    }
}

sub reversecomplement () {
    ($sq) = @_;
    @A = split(//,$sq);
    $rev = "";
    for($i=@A-1; $i>=0; $i--) {
        $flag = 0;
        if($A[$i] eq 'A') {
            $rev = $rev . "T";
            $flag = 1;
        }
        if($A[$i] eq 'T') {
            $rev = $rev . "A";
            $flag = 1;
        }
        if($A[$i] eq 'C') {
            $rev = $rev . "G";
            $flag = 1;
        }
        if($A[$i] eq 'G') {
            $rev = $rev . "C";
            $flag = 1;
        }
        if($flag == 0) {
            $rev = $rev . $A[$i];
        }
    }

    return $rev;
}

