#!/usr/bin/perl -w
### Convert ALLHiC output AGP file to ALLMAPS input csv file
print "Convert ALLHiC output AGP file to ALLMAPS input csv file\n";
die "Usage: perl $0 groups.agp\n" if(!defined $ARGV[0]);

my $agp = $ARGV[0];
open(OUT, "> hic.csv") or die"";
print OUT "Scafffold ID,scaffold position,LG,genetic position\n";
open(IN, "grep -v 'contig' $agp|") or die"";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	if($data[8] eq "+"){
		$a = $data[6]; $b = $data[7];
	}elsif($data[8] eq "-"){
		$a = $data[7]; $b = $data[6];
		}
	print OUT "$data[5],$a,$data[0],$data[1]\n";
	print OUT "$data[5],$b,$data[0],$data[2]\n";
	}
close IN;
close OUT;

