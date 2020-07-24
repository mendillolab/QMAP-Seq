#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;

my $seqfile = $ARGV[0];

my $prefix = "ACACCG";
my $sgRNAlen = 20;
my $interval = 102;
my $cellLineLen = 8;

# This script takes a set of fastq-formatted reads in an uncompressed file.
# It parses the reads using the parameters above, looking for the presence
# of the prefix, then a string of length $sgRNAlen, then a string of length
# $interval, and finally a string of length $cellLineLen.  It stores the
# identity of the $sgRNAlen string and the $cellLineLen string, and counts
# the number of times each particular combination of sgRNA barcode and cell
# line barcode are observed.  After running over all of the reads, it outputs
# to standard out all of the pairs of barcodes observed and their total
# counts.

if ($#ARGV <0){
    die "Usage:  $0 <fastqfile>\n";
}
$seqfile =~ s/.fq$//g;
$seqfile =~ s/\_trimmed$//g;
my $lineCount = 5;
my ($header,$seq,$qual,$header2);
my %counts;


open (FQ,$seqfile);
print STDERR "Sequence file: $seqfile\n";

my $lost = 0;
my $noAdaptor = 1;
my %sgRNAs;
my %cellLines;

# Read in the Fastq file.
while(<FQ>){
    chomp $_;
    if (($_ =~ /^\@/)&& ($lineCount == 5)){ $header = $_; $lineCount = 1; $noAdaptor = 1;}
    elsif (($lineCount == 2)){ $seq = $_;}
    elsif (($lineCount == 3)){ $header2 = $_}
    elsif ($lineCount == 4){
	# If you are on a sequence line, parse it for barcodes.
	$qual = $_;
	my ($sgRNA,$cellLine) = "";
	if ($seq =~ /\w*$prefix([ACTGN]{$sgRNAlen})[ACTGN]{$interval}([ACTGN]{$cellLineLen})[ACTGN]*/) {
	    $noAdaptor = 0;
	    $sgRNA = $1;
	    $cellLine = $2;
	    # Increment the counter for that combination of bar counts.
	    if (exists($counts{$sgRNA}{$cellLine})){
		$counts{$sgRNA}{$cellLine} ++;
	    } else {
		$counts{$sgRNA}{$cellLine} =1;
		$sgRNAs{$sgRNA} = 1;
		$cellLines{$cellLine} = 1;
	    }
	} 
    }
    if ($lineCount ==4 && $noAdaptor == 1){ 
#	print STDERR "$header lost\n"; 
	$lost++;
    } elsif ( $lineCount == 4 && $noAdaptor == 0) { 
	#print STDERR "$header found!\n";
    }
    $lineCount++;
}
close(FQ);

print STDERR "$lost Reads could not be assigned to sgRNA.\n";

# Print out all of the combinations of sgRNA and cell line barcodes found
# and their counts, tab-separated.
foreach my $sgRNA (keys(%sgRNAs)){
    foreach my $cellLine (keys(%cellLines)){
	if (exists($counts{$sgRNA}{$cellLine})){
	    print "$counts{$sgRNA}{$cellLine}\t$sgRNA\t$cellLine\n";
	}
    }
}
