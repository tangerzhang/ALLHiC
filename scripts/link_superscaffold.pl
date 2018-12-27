#!/usr/bin/perl -w

my %namedb;
my %removedb;
while(<DATA>){
	chomp;
	my ($id,$name) = (split/\s+/,$_)[0,1];
	$namedb{$name} = $id;
	my @data = split(/\s+/,$_);
	my $key  = "";
	foreach my $i (2..$#data){
		my ($sa,$sb) = sort ($data[1],$data[$i]);
		$key   = $sa."	".$sb;
		$removedb{$key}++;
		}
	}
	
my %infordb;
open(IN, "grep -v 'tig' all.clm|") or die"";
<IN>;
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	my $scf1 = $data[0];
	my $scf2 = $data[1];
	   $scf1 =~ s/[+|-]//g;
	   $scf2 =~ s/[+|-]//g;
#	my ($s1,$s2) = sort ($scf1,$scf2);
	my $key1  = $scf1."	".$scf2;
	my $key2  = $scf2."	".$scf1;
	next if(exists($removedb{$key1}));
	if(!exists($infordb{$key1})){
		$infordb{$key1} = $data[4];
	}elsif(exists($infordb{$key1}) and $data[4]>$infordb{$key1}){
		$infordb{$key1} = $data[4];
		}
	if(!exists($infordb{$key2})){
		$infordb{$key2} = $data[4];
	}elsif(exists($infordb{$key2}) and $data[4]>$infordb{$key2}){
		$infordb{$key2} = $data[4];
		}	
	}
close IN;

my %bestdb;
open(OUT, "> tmp.txt") or die"";
foreach my $key (keys %infordb){
	my ($sa,$sb) = split(/\s+/,$key);
	my $ida      = $namedb{$sa};
	my $idb      = $namedb{$sb};
	print OUT "$ida	$idb	$sa	$sb	$infordb{$key}\n";
	if(!exists($bestdb{$ida}->{$idb})){
		$bestdb{$ida}->{$idb} = $infordb{$key};
	}elsif($infordb{$key}>$bestdb{$ida}->{$idb}){
		$bestdb{$ida}->{$idb} = $infordb{$key};
		}
	if(!exists($bestdb{$idb}->{$ida})){
		$bestdb{$idb}->{$ida} = $infordb{$key};
	}elsif($infordb{$key}>$bestdb{$idb}->{$ida}){
		$bestdb{$idb}->{$ida} = $infordb{$key};
		}	
	}
close OUT;

open(OUT, "> best_link.txt") or die"";
my $ln = 0;
my %linkdb;
open(IN, "sort -k5,5nr -k1,1n tmp.txt|") or die"";
while(<IN>){
	chomp;
	$ln++;
  my @data = split(/\s+/,$_);
  my $ida  = $data[0];
  my $idb  = $data[1];
  my $key  = $ida."	".$idb;
  if($ln==1){
  	$linkdb{$key} = $_;
  	$tmpdb{$ida}++;
  	$tmpdb{$idb}++;
  }else{
  	next if(exists($tmpdb{$ida}) or exists($tmpdb{$idb}));
   	$linkdb{$key} = $_;
  	$tmpdb{$ida}++;
  	$tmpdb{$idb}++; 	
  	}
	}
close IN;

foreach my $key (keys %linkdb){
	print OUT "$linkdb{$key}\n";
	}

close OUT;


### Below are the information that listed allelic super-scaffolds for each target.
#Format:
#ID	target	allelic_superscaffold1	allelic_superscaffold2 ...
__DATA__
1	group1	group2	group4	group6	group8	group9	group11	group14	group15	group16
2	group2	group1	group3	group4	group5	group6	group7	group8	group9	group10	group11	group12	group13	group14	group15	group16
3	group3	group2	group5	group7	group9	group10	group11	group12	group13	group14
4	group4	group1	group3	group6	group8	group9	group11	group14	group15	group16
5	group5	group2	group3	group7	group9	group10	group11	group12	group13	group14
6	group6	group1	group2	group4	group8	group9	group11	group14	group15	group16
7	group7	group2	group3	group5	group9	group10	group11	group12	group13	group14
8	group8	group1	group3	group4	group6	group9	group11	group14	group15	group16
9	group9	group1	group3	group4	group5	group6	group7	group8	group2	group10	group11	group12	group13	group14	group15	group16
10	group10	group2	group3	group5	group7	group9	group11	group12	group13	group14
11	group11	group1	group3	group4	group5	group6	group7	group8	group9	group10	group2	group12	group13	group14	group15	group16
12	group12	group2	group3	group5	group7	group9	group10	group11	group13	group14
13	group13	group12	group2	group3	group5	group7	group9	group10	group11
14	group14	group12	group2	group3	group5	group7	group9	group10	group11	group16
15	group15	group1	group2	group4	group6	group8	group9	group14
16	group16	group1	group2	group4	group6	group8	group9	group11	group14	group15
