#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use List::Util qw( min max );

my $countTable = $ARGV[0];
my $adaptorFile = $ARGV[1];
my $maxCount = $ARGV[2];

if ($#ARGV <1){
    die "Usage:  $0 <countTable> <adaptorFile> <maxMismatch>\n";
}

my %clBarcodes;
my %sgBarcodes;
my %labeledTable;
open(IN,$adaptorFile);
while(<IN>){
    chomp $_;
    if ($_ !~ /^Cell Line/){
	my ($cellLine,$sgRNA,$sgBarcode,$clBarcode) = split(/\t/,$_);
	print STDERR "CL:$cellLine\tSG:$sgRNA\tSGB:$sgBarcode\tCLB:$clBarcode\n";
	$clBarcodes{$clBarcode} = $cellLine;
	$sgBarcodes{$sgBarcode} = $sgRNA;
    }
}
my @knownClBarcodes = keys(%clBarcodes);
my @knownSgBarcodes = keys(%sgBarcodes);
print STDERR "Cell line Barcodes:\n@knownClBarcodes\n";
print STDERR "sgRNA Barcodes:\n@knownSgBarcodes\n";

my %results;
my $ambiguous = 0;
my $totalCount = 0;
open(TB,$countTable);
while(<TB>){
    chomp $_;
    my ($pairCount,$sgBarcode,$clBarcode) = split(/\t/,$_);
    my $sgLabel = "";
    my $clLabel = "";
    if (exists($clBarcodes{$clBarcode})){
	print STDERR "Found Cell line $clBarcodes{$clBarcode}\n";
	$clLabel = $clBarcodes{$clBarcode};
    } else {
	print STDERR "Did not find Cell line $clBarcode\n";
	my $mmCount = 0;
	my $matches = 0;
	# Now looking for mismatched barcodes
	foreach my $knownBarcode (@knownClBarcodes){
	    $mmCount = mismatch_count($knownBarcode,$clBarcode);
	    if ($mmCount <= $maxCount){
		print STDERR "Cell line code $clBarcode is $mmCount mismatches from $knownBarcode, $clBarcodes{$knownBarcode}\n";
		$clLabel = $clBarcodes{$knownBarcode};
		$matches++;
	    }
	}
	if ($matches > 1){ 
	    print STDERR "ERROR: Hit more than one Cell line barcode with $clBarcode\n";
	    $clLabel = "";
	    $ambiguous++;
	}
    }
    if (exists($sgBarcodes{$sgBarcode})){
	print STDERR "Found sgRNA $sgBarcodes{$sgBarcode}\n";
	$sgLabel = $sgBarcodes{$sgBarcode};
    } else {
	print STDERR "Did not find sgRNA $sgBarcode\n";
	my $mmCount = 0;
	my $matches = 0;
	foreach my $knownBarcode (@knownSgBarcodes){
	    my $mmCount = mismatch_count($knownBarcode,$sgBarcode);
	    if ($mmCount <= $maxCount){
		print STDERR "Cell line code $sgBarcode is $mmCount mismatches from $knownBarcode, $sgBarcodes{$knownBarcode}\n";
		$sgLabel = $sgBarcodes{$sgBarcode};
		$matches ++;
	    }
	}
	if ($matches > 1){ 
	    print STDERR "ERROR: Hit more than one sgRNA barcode with $sgBarcode\n";
	    $sgLabel = "";
	    $ambiguous++;
	}
    }
    if (($sgLabel ne "") && ($clLabel ne "")){
	if (exists($labeledTable{$sgLabel}{$clLabel})){
	    $labeledTable{$sgLabel}{$clLabel} += $pairCount;
	} else {
	    $labeledTable{$sgLabel}{$clLabel} = $pairCount;
	}
    }
}
@knownClBarcodes = sort { $a cmp $b } @knownClBarcodes;
@knownSgBarcodes = sort { $a cmp $b } @knownSgBarcodes;

foreach my $knownClBarcode (@knownClBarcodes){
    my $clLabel = $clBarcodes{$knownClBarcode};
    foreach my $knownSgBarcode (@knownSgBarcodes){
	my $sgLabel = $sgBarcodes{$knownSgBarcode};
	if (exists($labeledTable{$sgLabel}{$clLabel})){
	    print "$knownClBarcode\t$knownSgBarcode\t$clLabel\t$sgLabel\t$labeledTable{$sgLabel}{$clLabel}\n";
	    $totalCount+=$labeledTable{$sgLabel}{$clLabel};
	} else {
	    print "$knownClBarcode\t$knownSgBarcode\t$clLabel\t$sgLabel\t0\n";
	}
    }
}

print STDERR "Ambiguous barcodes = $ambiguous\n";
print STDERR "Counted barcode pairs = $totalCount\n";


sub mismatch_count($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }
    
