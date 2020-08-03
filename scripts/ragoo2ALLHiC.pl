#!/usr/bin/perl -w

use Getopt::Std;
getopts "l:r:b:e:";


if ((!defined $opt_l)|| (!defined $opt_r) ||(!defined $opt_b)) {
    die "************************************************************************
    Usage: perl ragoo2ALLHiC -l orderings.list -r draft.asm.fasta -b sample.clean.bam 
      -h : help and usage.
      -l : ordering.list contains a list of output files from ragoo
      -r : draft contig assembly
      -b : sample.clean.bam
      -e : restriction sites, optional, default GATC
           MboI: GATC; HindIII: AAGCTT
************************************************************************\n";
}else{
  print "************************************************************************\n";
  print "Version demo\n";
  print "Copyright to Tanger\n";
  print "RUNNING...\n";
  print "************************************************************************\n";
	
	}

$opt_e = (defined $opt_e)?$opt_e:"GATC";


if(!(-e "draft.asm.fasta")){
	system("ln -s $opt_r ./draft.asm.fasta");
}else{
	print "check draft.asm.fasta file, exist\n";
	}

if(!(-e "sample.clean.bam")){
	system("ln -s $opt_b ./sample.clean.bam");
}else{
	print "check sample.clean.bam file, exist\n";
	}


my $num_g = 0;
my %cntdb = ();
open(IN, $opt_l) or die"";
while(<IN>){
	chomp;
	$num_g++;
	my @linedb = split(/\n/,$_);
	foreach my $file (@linedb){
		$gid = (split/\//,$file)[-1];
		$gid =~ s/_orderings.txt//g;
		open(my $fh, $file) or die"";
		while(<$fh>){
			chomp;
			my $ctg = (split/\s+/,$_)[0];
			$cntdb{$gid}->{$ctg}++;
			}
		close $fh;
		}
	}
close IN;


open(OUT, ">clusters.txt") or die"";
print OUT "#Group	nContigs	Contigs\n";
foreach my $g (sort keys %cntdb){
	my $num = keys %{$cntdb{$g}};
	print OUT "$g	$num	";
	foreach my $c (keys %{$cntdb{$g}}){
		print OUT "$c	";
		}
	print OUT "\n";
	}
close OUT;

print "#### Counting restriction sites from draft assembly\n";
print "allhic extract sample.clean.bam draft.asm.fasta --RE $opt_e\n...\n\n";
system("allhic extract sample.clean.bam draft.asm.fasta --RE $opt_e");

print "### Rescue unanchored contigs\n";
my $countRE = "sample.clean.counts_".$opt_e.".txt";
print "ALLHiC_rescue -r draft.asm.fasta -b sample.clean.bam -c clusters.txt -i $countRE\n...\n\n";
system("ALLHiC_rescue -r draft.asm.fasta -b sample.clean.bam -c clusters.txt -i $countRE -m 1");

foreach my $i (1..$num_g){
	my $gn = "group".$i.".txt";
	print "### Scaffolding $gn\n";
	print "allhic optimize $gn sample.clean.clm\n...\n\n";
	system("allhic optimize $gn	sample.clean.clm");
	}

print "### Build ALLHiC assembly\n";
system("ALLHiC_build draft.asm.fasta");

system("Done ...\n");

