#!/usr/bin/perl -w
use strict;



# PreprocessSAMs.pl
#
# Syntax: PreprocessSAMs.pl <sam or bam filename> <draft assembly fasta>
#
# This Perl script prepares a SAM/BAM file for use with Lachesis.
# Specifically, it pre-processes the file with bedtools, samtools, picard to remove redundant, chimeric, and/or uninformative read pairs.
# This creates a dataset of Hi-C links with as strong a signal as possible, and it's also as small as possible, so as to reduce I/O runtime in Lachesis.
# (NOTE: As of August 24, 2013, I'm no longer removing PCR duplicates.  Picard's MarkDuplicates is extremely slow and resource-intensive - far more so than
# the runtime benefit in Lachesis of having fewer reads.  I don't think it's removing PCR duplicates properly, nor do I think PCR duplicate removal is even
# necessary - http://seqanswers.com/forums/showthread.php?t=6854).
#
# This script will determine whether the file is a SAM or a BAM file, and then run the following commands:
#
# COMMAND                                OUTPUT FILENAME                               WHAT THE COMMAND DOES
# make_bed_around_RE_site.pl             <fasta>.near_<RE>.<range>.bed                 Prepare the bed file for bedtools intersect (next command)
# bedtools intersect                     <head>.REduced.bam                            Remove all reads that aren't within 500 bp of a restriction site
### picard SortSam.jar                     <head>.REduced.sort_coord.bam                 Sort the file in coordinate order so PCR duplicates can be removed
### picard MarkDuplicates.jar              <head>.REduced.sort_coord.nodups.bam          Remove PCR duplicates
### picard SortSam.jar                     <head>.REduced.nodups.bam                     Sort the file in query-name order so Lachesis can read it
# samtools view -F12                     <head>.REduced.nodups.paired_only.bam         Filter out all pairs in which both reads are not aligned
# samtools flagstat                      <head>.REduced.nodups.paired_only.flagstat    Make a flagstat file that describes the contents of the BAM file
#
#
# The final output file will be <head>.REduced.paired_only.bam.  This is what should be entered into the Lachesis INI file under the key "SAM_FILES".
#
# To pre-process several SAM/BAM files in parallel, use the script PreprocessSAMs.sh, which can be submitted to a cluster via qsub.
#
# Josh Burton
# July 2013



################################
#                              #
#   USER-DEFINED PARAMETERS    #
#                              #
################################


my $dry_run = 0; # if true, just print the commands to be run - don't actually run them
#my $RE_site = 'AAGCTT'; # the restriction enzyme site at which the DNA was cut for the Hi-C experiment

# Paths to the necessary scripts and software packages.
my $make_bed_around_RE_site_pl = 'make_bed_around_RE_site.pl';
my $bedtools = 'bedtools';
my $samtools = 'samtools';
#my $mem = "16G";
#my $picard_head = "java -d64 -Xmx$mem -jar /net/shendure/vol10/jnburton/extern/picard-tools-1.50/";



################################
#                              #
#         SUBROUTINES          #
#                              #
################################

# Print and then run a command in bash (unless $dry_run, in which case just print it.)
# First argument: the command to run.
# Second argument (optional): the file to redirect stdout to.
sub run_cmd(@) {
    
    my ($cmd,$redirect) = @_;
    
    print localtime() . ": PreprocessSAMs.pl: $cmd\n";
    
    return if $dry_run;
    
    if ($redirect) { system ( "$cmd > $redirect" ) }
    else           { system ( $cmd ); }
}




################################
#                              #
#     CONTROL STARTS HERE      #
#                              #
################################


# Get the command-line arguments, or check syntax.
if ( @ARGV != 3 ) {
    print STDERR "\nPreprocessSAMs.pl: A script to prepare SAM or BAM files for use with Lachesis.\n\nSyntax: $0 <sam-or-bam-filename> <draft-assembly-fasta> enzyme(HINDIII/MBOI)\n\n";
    exit;
}


# Get the input filenames, and check that they actually exist.
my ( $SAM, $fasta) = @ARGV;
unless ( -e $SAM ) {
    print STDERR "$0: Can't find input SAM/BAM file `$SAM`\n";
    exit;
}
unless ( -e $fasta) {
    print STDERR "$0: Can't find draft assembly file `$fasta`\n";
    exit;
}

$ARGV[2] = uc $ARGV[2];
my $RE_site;
if($ARGV[2] eq "HINDIII"){
  $RE_site = 'AAGCTT';
  }elsif($ARGV[2] eq "MBOI"){
  $RE_site = 'GATC';
  }
# Find the input file's "head" and extension.
my ($head,$extension) = $SAM =~ /^(.*)\.(.*)$/;


# Examine the extension to determine whether this is a SAM or a BAM file.  If it's a SAM, convert it to BAM.  If it doesn't seem to be either, throw an error.
if    ( uc($extension) eq 'SAM' ) { run_cmd( "$samtools view -bS $SAM -o $head.bam" ); }
elsif ( uc($extension) eq 'BAM' ) {}
else {
    print STDERR "$0: Can't determine file type for input file `$SAM`.\nFilename should end in '.SAM' or '.BAM' (not case-sensitive.)\n";
    exit;
}


print "$0 @ARGV\n\n";



# COMMAND                                OUTPUT FILENAME                               WHAT THE COMMAND DOES
# make_bed_around_RE_site.pl             <fasta>.near_<RE>.<range>.bed                 Prepare the bed file for bedtools intersect (next command)
#
# Make the BED file for the restriction sites on the draft assembly.  This only needs to be done once.
my $BED_RE_file = "$fasta.near_$RE_site.500.bed";
run_cmd( "$make_bed_around_RE_site_pl $fasta $RE_site 500" ) unless -e $BED_RE_file;


# Do the pre-processing on this file.
#
# COMMAND                                OUTPUT FILENAME                               WHAT THE COMMAND DOES
# bedtools intersect                     <head>.REduced.bam                            Remove all reads that aren't within 500 bp of a restriction site
### picard SortSam.jar                     <head>.REduced.sort_coord.bam                 Sort the file in coordinate order so PCR duplicates can be removed
### picard MarkDuplicates.jar              <head>.REduced.sort_coord.nodups.bam          Remove PCR duplicates
### picard SortSam.jar                     <head>.REduced.nodups.bam                     Sort the file in query-name order so Lachesis can read it
# samtools view -F12                     <head>.REduced.paired_only.bam         Filter out all pairs in which both reads are not aligned
# samtools flagstat                      <head>.REduced.paired_only.flagstat    Make a flagstat file that describes the contents of the BAM file

my $opts = "VALIDATION_STRINGENCY=SILENT";
my $nodups = ""; # or ".nodups", if removing PCR duplicates

run_cmd( "$bedtools intersect -abam $head.bam -b $BED_RE_file > $head.REduced.bam" );
#run_cmd( "${picard_head}SortSam.jar $opts I=$head.REduced.bam O=$head.REduced.sort_coord.bam SO=coordinate" );
#run_cmd( "${picard_head}MarkDuplicates.jar $opts I=$head.REduced.sort_coord.bam O=$head.REduced.sort_coord.nodups.bam M=$head.REduced.sort_coord.dup_metrics AS=true REMOVE_DUPLICATES=true" );
#run_cmd( "${picard_head}SortSam.jar $opts I=$head.REduced.sort_coord.nodups.bam O=$head.REduced.nodups.bam SO=queryname" );
run_cmd( "$samtools view -F12 $head.REduced$nodups.bam -b -o $head.REduced$nodups.paired_only.bam" );
run_cmd( "$samtools flagstat $head.REduced$nodups.paired_only.bam > $head.REduced$nodups.paired_only.flagstat" );
