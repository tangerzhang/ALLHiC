#!/usr/bin/perl -w

die "Usage: perl $0 ref.gff3\n" if(!defined ($ARGV[0]));
my $refGFF = $ARGV[0];
open(IN, "grep 'gene' gmap.gff3 |") or die"";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	my $gene = $1 if(/Name=([^;\n]*)/);
	$infordb{$gene} .= $data[0]."	";
	}
close IN;

open(OUT, "> Allele.ctg.table") or die"";
open(IN, "awk '\$3==\"gene\"' $refGFF |") or die"";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	my $gene = $1 if(/Name=(\S+)/);
	   $gene =~ s/;.*//g;
	next if(!exists($infordb{$gene}));
	my @tdb = split(/\s+/,$infordb{$gene});
	my %tmpdb = ();
	map {$tmpdb{$_}++} @tdb;
	print OUT "$data[0]	$data[3]	";
	map {print OUT "$_	"} keys %tmpdb;
	print OUT "\n";
	}
close IN;
close OUT;
