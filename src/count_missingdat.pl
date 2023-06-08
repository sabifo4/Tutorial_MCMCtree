#!/usr/bin/perl

use strict;
use warnings;
use Cwd qw(getcwd);
use List::MoreUtils qw(uniq);

## Script that reads a PHYLIP alignment and counts the missing 
## taxa
##
## Usage:
## <path_to_script>/count_missingdat.pl <path_to_aln> 
##
## NOTE: If the user wants to run this script in a loop for several
##       directories, an example is given below: 
##
##       ```  
##       # Run from main directory where other subdirectories
##       # have the alignments
##       cd <subdirectory>
##       for i in `seq <START> <END>`
##       do
##           cd $i
##         	 perl <path_to_script>/count_missingdat.pl *aln
##           cd ..
##       done
##       ```
##
##       This code snippet assumes that you have a PHYLIP file as an input 
##       file.
##
## Contact information: <sandra.ac93@gmail.com>

#IMPORTANT! CREATE ".TXT" FILES IN UNIX, NOT WINDOWS!!!!! IT THEN ADDS THE "\r" CARRIAGE AND IT MAKES THE SCRIPT CRASH...
		
## Open the input files
open (INFILE1, "<$ARGV[0]") or die "Cannot open $ARGV[0] file: $!";
## Get second argument to set output file name
my $outname = $ARGV[0];
chomp($outname);
$outname =~ s/\.aln//;
$outname =~ s/..*\///;
print "Parsing alignment ".$outname." ... ... \n";
my $outname2 = $outname."_missdatpersp.txt";
my $outname3 = $outname."_avgmissdata.txt";
my $outdir = "out_count_NA";
## Create the output dir to save results
if( ! -d $outdir ){
	mkdir $outdir ;
}
open(OUT, ">$outdir/$outname2") or die "Cannot create the output file: $!";
open(OUT2, ">$outdir/$outname3") or die "Cannot create the output file: $!";

## Print header for output 1 file 
print OUT "Taxon\tLen_aln\tNum_gaps\tMissing_taxon\tperc_missdat\n";

## Get lines of input file
my @aln1 = <INFILE1>;

## Create global variables
my $count_sp = 0;
my $len_aln = 0;
my $num_sp = 0;
my $tmp_sp = "";
my $str_gaps = "";
my $tmp_gaps = "";
my $tmp_missdat = 0;

my $is_missing = "N";
my $count_missing = 0;
my $count_no_missing = 0;
my $str_gaps_no_missingtaxa = "";
my $str_gaps_missingtaxa = "";

my $count_lines = 0;
## Loop over all the lines in input file 
foreach my $line (@aln1){

	chomp($line);
	$count_lines += 1;
	# Skip blank lines!
	if( $line =~ /^$/ ){
		next;
	}
	# Have added this condition because it seems the second line 
	# has a weird character that is not detected by first condition 
	# and screws the rest of the parsing
	elsif ( $count_lines != 2 ){
		# Find header
		if( $line =~ /^[0-9]/ && $line !~ /[A-Za-z]/ ){
			$len_aln = $line;
			$len_aln =~ s/..* //; # This is modified from other version of count missdat
								  # There is only one space in this type of aln!
			$len_aln =~ s/\r//;
			$len_aln =~ s/\n//;
			chomp( $len_aln );
			$num_sp = $line;
			$num_sp =~ s/\s..*//;
			$num_sp =~ s/\r//;
			$num_sp =~ s/\n//;
			chomp( $num_sp );
			print "Length alignment: ".$len_aln."\n";
			print "Number of taxa: ".$num_sp."\n";
		}else{
			# Keep counting taxa 
			$count_sp += 1;
			#print "Line read: ".$count_lines."\n";
			# Get taxon name + sequence 	
			my @tmp_line = split( / {2,}/, $line );
			$tmp_sp = $tmp_line[0];
			$str_gaps = $tmp_line[1];
			$str_gaps =~ s/[ATGCatgc]//g;
			$str_gaps =~ s/\r//;
			$str_gaps =~ s/\n//;
			chomp( $str_gaps );
			# If $str_gaps eq $len_aln, flag as missing taxon 
			if ( length( $str_gaps ) == $len_aln ){
				$is_missing = "Y";
				$count_missing += 1;
				$str_gaps_missingtaxa .= $str_gaps;
			}else{
				# Now add only gaps that are not
				# part of missing taxa
				$is_missing = "N";
				$count_no_missing += 1;
				$str_gaps_no_missingtaxa .= $str_gaps;
			}
			# print "seq: ".$str_gaps."\n";
			# Calculate missing data per taxon
			$tmp_missdat = length( $str_gaps )*100/$len_aln;
			print "Taxon: ".$tmp_sp." | Num gaps: ".length( $str_gaps )." | Missing taxon?: ".$is_missing." | % Missing data: ".$tmp_missdat."\n";
			print OUT $tmp_sp."\t".$len_aln."\t".length( $str_gaps )."\t".$is_missing."\t".$tmp_missdat."\n";

		}
	}
}
			
## Check if there has been errors or missed 
## a line
if( $count_sp != $num_sp ){
	print "You must have missed reading one line!\n";
	print "The two numbers counted in two different vars\n";
	print "are not equal and this is a mistake!\n";
}
my $check_total = $count_missing + $count_no_missing;
if( $check_total != $num_sp ){
	print "You must have missed reading one line!\n";
	print "The two numbers counted in two different vars\n";
	print "are not equal and this is a mistake!\n";
}
			
# Compute average missing data
# Before removing taxa: 
#  num_gaps / (num_taxa * length_aln )
# Now, do the same but use only num_gaps that are not part 
# of missing taxa and num_taxa replace but only those taxa 
# that do not have all gaps !
print "\nNumber of missing taxa = ".$count_missing."\n";
print "Number of present taxa = ".$count_no_missing."\n";
print "Total number of taxa = ".$num_sp."\n";

print "\nEq: num_gaps / ( num_taxa * len_aln ) --> "."(".length( $str_gaps_missingtaxa )."+".length( $str_gaps_no_missingtaxa ).")/((".$count_missing."+".$count_no_missing.")*".$len_aln.")\n";
my $missing_taxa_before =  (length( $str_gaps_missingtaxa )+length( $str_gaps_no_missingtaxa ))/(($count_missing+$count_no_missing)*$len_aln);
print "Avg missing data before removing missing taxa = ".$missing_taxa_before." = ".($missing_taxa_before*100)."%\n";
print "Eq: num_gaps / ( num_taxa * len_aln ) --> ".length( $str_gaps_no_missingtaxa )."/(".$count_no_missing."*".$len_aln.")\n";
my $missing_taxa_after = length( $str_gaps_no_missingtaxa )/($count_no_missing*$len_aln);
print "Avg missing data after removing missing taxa = ".$missing_taxa_after." = ".($missing_taxa_after*100)."%\n\n";

print OUT2 "Num_missing_taxa\tNum_present_taxa\tNum_total_taxa\tAvg_missdat_with_missing_taxa\tAvg_missdat_without_missing_taxa\n";
print OUT2 $count_missing."\t".$count_no_missing."\t".$num_sp."\t".$missing_taxa_before."\t".$missing_taxa_after."\n";
print OUT2 "\nEq. before removing missing taxa: num_gaps / ( num_taxa * len_aln ) --> "."(".length( $str_gaps_missingtaxa )."+".length( $str_gaps_no_missingtaxa ).")/((".$count_missing."+".$count_no_missing.")*".$len_aln.")\n";
print OUT2 "Eq. after removing missing taxa:  num_gaps / ( num_taxa * len_aln ) --> ".length( $str_gaps_no_missingtaxa )."/(".$count_no_missing."*".$len_aln.")\n";
print OUT2 "\nMissing data before (total taxa): ".sprintf( "%.2f", ($missing_taxa_before*100) )." (".$num_sp.")\n";
print OUT2 "\nMissing data after (present taxa): ".sprintf( "%.2f", ($missing_taxa_after*100) )." (".$count_no_missing.")\n";

## Close files
close(OUT);
close(OUT2);
close(INFILE1);