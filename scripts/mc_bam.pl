#!/usr/bin/perl -w

use Getopt::Std;
getopts "b:a:r:";


if ((!defined $opt_b)|| (!defined $opt_a) ||(!defined $opt_r)) {
    die "************************************************************************
    Usage: mc_bam.pl -b mapping.bam -r groups.asm.fasta -a agp
      -h : help and usage.
           This script is used for modification the coordinates 
           in bam based on agp file
      -b : mapping.bam
      -r : reference genome, fasta format
      -a : agp file
************************************************************************\n";

}

my %posidb = ();
open(IN, $opt_a) or die"";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	next if($data[-1] eq "map");
	my $ga   = $data[1]; 
	my $gb   = $data[2]; 
	my $ta   = $data[6];
	my $tb   = $data[7];
	my $ctg  = $data[5];
	if($data[8] eq "+"){
		my $gi = $ga;
		foreach my $ti($ta..$tb){
			$posidb{$ctg}->{$ti} = $data[0].",".$gi;
			$gi++;
			}
	 }elsif($data[8] eq "-"){
	  my $gi = $gb;
	  foreach my $ti ($ta..$tb){
	  		$posidb{$ctg}->{$ti} = $data[0].",".$gi;
	  		$gi--;
	  		}
	  	}

	}
close IN;

#open(OUT, "> posi.txt") or die"";
#foreach my $ctg (keys %posidb){
#  foreach my $i(sort {$a<=>$b} keys %{$posidb{$ctg}}){
#	  print OUT "$ctg	$i	$posidb{$ctg}->{$i}\n";
#	  }
#	}
#close OUT;

my $outsam = "mc.sam";
my $outbam = "mc.bam";
open(OUT, "> $outsam") or die"";
open(IN, "samtools view $opt_b |") or die"";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	next if($data[6] eq "=");
	my $ctgA = $data[2]; 
	my $tiA  = $data[3]; 
	my $ctgB = ($data[6] eq "=")?$data[2]:$data[6];
	my $tiB  = $data[7];
	my $gidA; my $giA; my $gidB; my $giB;
	next if(!exists($posidb{$ctgA}->{$tiA}));
	next if(!exists($posidb{$ctgB}->{$tiB}));
  ($gidA,$giA) = split(/,/,$posidb{$ctgA}->{$tiA});
	($gidB,$giB) = split(/,/,$posidb{$ctgB}->{$tiB});
	$data[2]       = $gidA;
	$data[3]       = $giA;
	$data[6]       = ($gidB eq $gidA)?"=":$gidB;
	$data[7]       = $giB;
	map {print OUT "$_	"} @data;
	print OUT "\n";
	}
close IN;
close OUT;

system("samtools faidx $opt_r");
my $fai = $opt_r.".fai";
system("samtools view -bt $fai $outsam > $outbam");
