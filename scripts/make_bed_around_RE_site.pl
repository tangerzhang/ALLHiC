#!/usr/bin/perl -w
use strict;


# make_bed_around_restriction_site.pl: Make a BED file representing the regions around all occurrences of a restriction site.
#
# For syntax, run with no arguments.
#
# The output BED file is designed for use with bedtools intersect, as follows:
# bedtools intersect -abam [SRR.bam] -b [$BED_out] > [SRR.REduced.bam]
# samtools view -h [SRR.REduced.bam] > [SRR.REduced.sam]
# This restricts a SAM/BAM file to only include reads close to a restriction site, which is a good way to filter Hi-C data, according to Fig. 1b of this paper:
# http://www.nature.com/ng/journal/v43/n11/full/ng.947.html
# Also see PreprocessSAM.pl, which uses the output file.
#
# Josh Burton
# April 2013




if ( scalar @ARGV != 3 ) {
    
    # Report syntax.
    print "\nmake_bed_around_RE_site.pl\n\n";
    print "Find all occurrences of a motif in a genome.  Make a 'POS' file listing these occurrences, and also a BED file representing the regions around these occurrences.\n\n";
    print "SYNTAX:\tmake_bed_around_RE_site.pl <fasta> <motif> <range>\n";
    print "fasta:\tA fasta file representing a genome (reference or draft assembly.)\n";
    print "motif:\tA motif, typically a restriction site sequence (e.g., HindIII = AAGCTT, NcoI = CCATGG, Dpn1 = GATC).\n";
    print "range:\tA number representing how many bp around the sequence to include.  Recommend 500 based on Yaffe & Tanay, Nat. Genetics 2011.\n\n";
    print "OUTPUT FILES:\n";
    print "<fasta>.near_<motif>.<range>.bed\n";
    print "<fasta>.near_pos_of_<motif>.txt\n";
    print "\n";
    exit;
}



# Get command-line arguments.
my ( $FASTA_in, $motif_seq, $range ) = @ARGV;

my $verbose = 0;

# Convert the motif from a string into a regex.  Unroll the IUPAC codes from single letters into Perl-parseable regular expressions.
my $motif_regex = $motif_seq;
$motif_regex =~ s/R/\[AG\]/g;
$motif_regex =~ s/Y/\[CT\]/g;
$motif_regex =~ s/S/\[CG\]/g;
$motif_regex =~ s/W/\[AT\]/g;
$motif_regex =~ s/K/\[GT\]/g;
$motif_regex =~ s/M/\[AC\]/g;
$motif_regex =~ s/B/\[CGT\]/g;
$motif_regex =~ s/D/\[AGT\]/g;
$motif_regex =~ s/H/\[ACT\]/g;
$motif_regex =~ s/V/\[ACG\]/g;
$motif_regex =~ s/N/\[ACGT\]/g;




# Derive an output filename.
my $BED_out = "$FASTA_in.near_$motif_seq.$range.bed";
my $POS_out = "$FASTA_in.pos_of_$motif_seq.txt";


# Determine how many letters needed to be added to each line in order to find instances of the sequence that bridge lines in the fasta.
my $N_prev_chars = length($motif_seq) - 1;


my $contig_name = '';
my $offset = 0;
my $prev_chars;
my @motif_positions;
my $N_motifs_found = 0;


# Open the input fasta file and read through it line-by-line.
print localtime() . ": Reading file $FASTA_in...\n";
open IN, '<', $FASTA_in or die "Can't find file `$FASTA_in'";
open BED, '>', $BED_out or die;
open POS, '>', $POS_out or die;

while (<IN>) {
    my $line = $_;
    chomp $line;
    
    # If this is a header line, we're done with this contig/chromosome (unless we just started), and start a new contig/chromosome.
    if ( $line =~ /^\>(\S+)/ ) {
	
	# The hash %motif_positions contains all positions on the (now complete) old contig at which this motif appears.
	# Convert this list of positions to a set of BED lines, as necessary.
	my ( $prev_start, $prev_end ) = (-1,-1);
	foreach my $pos ( @motif_positions ) {
	    if ( $prev_end == -1 ) {
		$prev_start = $pos;
		$prev_end   = $pos;
	    }
	    if ( $prev_end + 2*$range < $pos ) {
		$prev_start =         $range if $prev_start < $range;
		$prev_end = $offset - $range if $prev_end > $offset - $range; # prevent overflow past the end of the contig/chromosome
		print BED "$contig_name\t", $prev_start - $range, "\t", $prev_end + $range, "\n";
		$prev_start = $pos;
	    }
	    #print "pos = $pos\n";
	    $prev_end = $pos;
	}
	
	# Print the final BED line for this contig/chromosome.
	if (@motif_positions) {
	    $prev_start =         $range if $prev_start < $range;
	    $prev_end = $offset - $range if $prev_end > $offset - $range; # prevent overflow past the end of the contig/chromosome
	    print BED "$contig_name\t", $prev_start - $range, "\t", $prev_end + $range, "\n";
	}
	
	# Get the new contig's name.
	$contig_name = $1;
	print localtime() . ": $contig_name\n" if $verbose;
	print POS ">$contig_name\n";
	
	# Reset other contig-related variables.
	$offset = 0;
	$prev_chars = '';
	@motif_positions = ();
    }
    
    # Otherwise, read through this contig/chromosome.
    else {
	if ( $offset != 0 ) { die unless $prev_chars; }
	
	my $verbose = 0;
	
	# Look for instances of this motif in this line of the fasta (including the overlap characters from the previous line, tacked on at the beginning.)
	my $motif_loc = -1;
	my $target_str = "$prev_chars" . uc $line;
	
	my @matches;
	while ($target_str =~ /$motif_regex/g ) {
	    
	    # Every iteration in this loop represents a new match to the motif regex in the terget string.
	    my $motif_loc = $-[0];
	    
	    # Adjust the location so it properly describes the 0-indexed motif position in this contig.
	    # Then add it to the list of contig positions at which the motif has been seen.
	    $N_motifs_found++;
	    my $true_motif_loc = $motif_loc + $offset - length $prev_chars; # adjust index so it properly describes the 0-indexed motif position in this contig
	    push @motif_positions, $true_motif_loc;
	    
	    print "$contig_name\t$offset\t$prev_chars\t->\t$motif_loc\n" if $verbose;
	    print POS "$true_motif_loc\n";
	}
	
	
	# TODO: remove
	while (0) {
	    $motif_loc = index "$prev_chars$line", $motif_seq, $motif_loc + 1;
	    last if ( $motif_loc == -1 ); # no more instances found
	    
	    # Found a motif!  Add its index to the list of contig positions at which the motif has been seen.
	    $N_motifs_found++;
	    my $true_motif_loc = $motif_loc + $offset - length $prev_chars; # adjust index so it properly describes the 0-indexed motif position in this contig
	    push @motif_positions, $true_motif_loc;
	    
	    print "$contig_name\t$offset\t$prev_chars\t->\t$motif_loc\n" if $verbose;
	    print POS "$true_motif_loc\n";
	}
	
	
	# Save the last few characters of this line, so that they can be appended onto the next line in a search for the sequence.
	my $line_len = length $line;
	$prev_chars = substr( $line, $line_len - $N_prev_chars );
	$offset += $line_len;
    }
    
    
}

close IN;
close BED;
close POS;


print localtime() . ": Done!  Found $N_motifs_found total instances of motif $motif_seq.  Created files:\n";
print "$BED_out\n$POS_out\n";
