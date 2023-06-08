#!/usr/bin/perl

use strict;
use warnings;
 
## Script that reads a FASTA alignment and converts it into 
## PHYLIP format.
##
## Usage:
## <path_to_script>/FASTAtoPHYL.pl <path_to_aln>
##
## NOTE: This code snippet assumes that you have a PHYLIP file generated
##       with in-house PERL script `FASTAtoPHYL.pl`.
##       In other words, it expectes 6 spaces between taxa and sequence.
##
## Contact information: <sandra.ac93@gmail.com>

#IMPORTANT! CREATE ".TXT" FILES IN UNIX, NOT WINDOWS!!!!! IT THEN ADDS THE "\r" CARRIAGE AND IT MAKES THE SCRIPT CRASH...
		
## Open the input files
open (INFILE1, "<$ARGV[0]") or die "Cannot open $ARGV[0] file: $!";
## Get second argument to set output file name
my $outname = $ARGV[0];
chomp($outname);
$outname =~ s/\.phy//;
$outname =~ s/..*\///;
print "Parsing alignment ".$outname." ... ... \n";
my $outname2 = $outname.".fasta";
## Open the output file to save the alignment file 
open(OUT, ">$outname2") or die "Cannot create the output file: $!";
open(OUT2, ">log_lenseq.txt") or die "Cannot create log file: $!";

## Get lines of input file
my @aln1 = <INFILE1>;

## Create variables 
my $count = 0;
my $lenseq = 0;
my $species = "";
my $species = "";
my $seq = "";

## Loop over all the lines in input file 
foreach my $line (@aln1){

	chomp($line);
	
	# Ommit header
	if( $line =~ /      / ){
		$species = $line;
		$seq = $line;
		# Keep track of length
		$lenseq = length($seq);
		# Get taxon name and sequence
		$species =~ s/      ..*//;
		$seq =~ s/..*      //;
		$count += 1;
		print OUT2 "Taxa: ".$species."\t"."Length of sequence: ".$lenseq."\n";
		print OUT ">".$species."\n".$seq."\n"
	}

}
					
print "Total no. species parsed: ".$count."\nLast seq length: ". $lenseq."\n";

## Close files
close(OUT);
close(OUT2);
close(INFILE1);