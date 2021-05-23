#!/usr/bin/perl -w

die "Usage: perl $0 chr.agp\n" if(!defined $ARGV[0]);
my %infordb;
my $cnt = 0;
open(IN, "grep -v contig $ARGV[0]|") or die"";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	my $chrn = $data[0];
	next if(!($chrn=~/Chr/) and !($chrn=~/group/));
	if(!exists($infordb{$chrn})){
		$cnt = 1;
		$infordb{$chrn}->{$cnt} .= $data[5]."".$data[8];
		}else{
			$cnt++;
			$infordb{$chrn}->{$cnt} .= $data[5]."".$data[8];
			}
	}
close IN;


foreach my $c (sort keys %infordb){
	my $outfile = $c.".tour";
	open(my $out, ">$outfile") or die"";
	foreach my $i (sort {$a<=>$b} keys %{$infordb{$c}}){
		print $out "$infordb{$c}->{$i} ";
		}
	close $out;
	}
