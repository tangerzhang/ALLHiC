#!/usr/bin/perl -w
use Getopt::Std;
getopts "c:b:o:";


if ((!defined $opt_c)|| (!defined $opt_b)||(!defined $opt_o) ) {
    die "************************************************************************
    Usage: bam2net.pl -c draft.asm.fasta -b file.bam -o out.net
      -h : help and usage.
      -c : draft.asm.fasta
      -b : mapping.bam
      -o : output
************************************************************************\n";
}
my $bam    = $opt_b;
my $refSeq = $opt_c;

open(IN, $refSeq) or die"";
my $name;
while(<IN>){
	chomp;
	if(/>/){
		$name = $_;
		$name =~ s/>//g;
	}else{
		$refdb{$name} .= $_;
		}
	}
close IN;

foreach $name (keys %refdb){
  $refdb{$name} =~ s/\s+//g;
	}

my %infordb;
open(IN, "samtools view $bam |") or die"";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	next if($data[6] eq "=");
	next if($data[6] eq "*");
        my ($ctg1,$ctg2) = sort ($data[2],$data[6]);
	$infordb{$ctg1}->{$ctg2}++;
	}
close IN;

open(OUT, "> $opt_o") or die"";
print OUT "ctg1	ctg1_size	ctg2	ctg2_size	signalDensity\n";
foreach my $ctg1(keys %infordb){
	my $len1 = length $refdb{$ctg1};
	foreach my $ctg2(keys %{$infordb{$ctg1}}){
		my $len2  = length $refdb{$ctg2};
                my $normL = ($len1 + $len2)/100000;
                my $sigD  = $infordb{$ctg1}->{$ctg2}/$normL;
                   $sigD  = sprintf("%.2f",$sigD);
		print OUT "$ctg1	$len1	$ctg2	$len2	$sigD\n";
		}
	}
close OUT;
