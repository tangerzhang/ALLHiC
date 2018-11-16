#!/usr/bin/perl -w

use Getopt::Std;
getopts "i:p:r:g:";


if ((!defined $opt_i)|| (!defined $opt_p) || (!defined $opt_r)|| (!defined $opt_g)) {
    die "************************************************************************
    Usage: perl $0 -i blast.out -p polyploid -r ref.gff3 -g target.gff3 
      -h : help and usage.
      -i : blast.out
      -p : number of alleles
      -r : reference.gff3, annotation from close relative species
      -g : target.gff3, annotation from target species

************************************************************************\n";
}

### Parameter reading
my $blast     = $opt_i;
my $polyn     = $opt_p;
my $rGFF      = $opt_r;
my $tGFF      = $opt_g;
my $geneTable = "Allele.gene.table";
my $ctgTable  = "Allele.ctg.table";

my %infordb;
my $count = 0;
open(IN, "sort -k2,2 -k12,12nr $blast|") or die"";
while(<IN>){
	chomp;
	my @data  = split(/\s+/,$_);
	my $tgene = $data[0];
	my $rgene = $data[1];
	my $bits  = $data[11];
	if(!exists($infordb{$rgene})){
		$count = 1;
		$infordb{$rgene}->{$count} = $tgene;
	}else{
		$count++;
		next if($count>$polyn);
		$infordb{$rgene}->{$count} = $tgene;
		}
	
	}
close IN;

my %tdb;  ### store target genome gff information, e.g het rice
open(IN, "awk '\$3==\"gene\"' $tGFF | ") or die"";
while(<IN>){
	chomp;
	my @data     = split(/\s+/,$_);
	my $tgene    = $1 if(/Name=(\S+)/);
	   $tgene    =~ s/;.*//g;
	$tdb{$tgene} = $data[0];
	}
close IN;


open(OUT, "> $geneTable") or die"";
open(IN, "awk '\$3==\"gene\"' $rGFF |sort -k1,1 -k4,4n |") or die"";
while(<IN>){
	chomp;
	my @data  = split(/\s+/,$_);
	my $rgene = $1 if(/Name=(\S+)/);
	$rgene    =~ s/;.*//g;
	next if(!exists($infordb{$rgene}));
	print OUT "$rgene	$data[0]	$data[3]	";
	foreach my $i(sort {$a<=>$b} keys %{$infordb{$rgene}}){
		my $tgene = $infordb{$rgene}->{$i};
		   $tctg  = $tdb{$tgene}; 
		print OUT "$tgene,$tctg	";     ###print out target gene order and contig name
		}
	print OUT "\n";
	}
close IN;

close OUT;


my %alleledb;
my $ln = 0; ###store line number
open(IN, "Allele.gene.table") or die"";
while(<IN>){
	chomp;
	$ln++;
	my @data = split(/\s+/,$_);
	my %tmpdb = ();
	foreach my $i(3..$#data){
		my $ctg = (split/,/,$data[$i])[1];
		$tmpdb{$ctg}++;
		}
	map {$alleledb{$ln}->{'ctg'} .= $_." "} keys %tmpdb;
	$alleledb{$ln}->{'chrn'}      = $data[1];
	$alleledb{$ln}->{'posi'}      = $data[2];
	}
close IN;

open(OUT, "> remove.log") or die"";
my %removedb = ();
for(my $i=2;$i<=$ln;$i++){
	my $chrI = $alleledb{$i}->{'chrn'};
	my $ctgI = $alleledb{$i}->{'ctg'};
	my $chrR; my $ctgR; my $R;
	for(my $j=1;$j<$i;$j++){
		next if(exists($removedb{$j}));
		my $chrJ = $alleledb{$j}->{'chrn'};
		next if($chrI ne $chrJ);
		my $ctgJ = $alleledb{$j}->{'ctg'};
		my $flag = & compare($ctgI,$ctgJ);
		print OUT "$i	$chrI	$ctgI	$j	$chrJ	$ctgJ	$flag\n" if($flag==1);
### flag=1, remove
		$removedb{$i}++ if($flag==1);
		}
	
	}
close OUT;


open(OUT, ">$ctgTable") or die"";
$ln = 0;
open(IN, $geneTable) or die"";
while(<IN>){
	chomp;
	$ln++;
	next if(exists($removedb{$ln}));
	my @data = split(/\s+/,$_);
	print OUT "$data[1]	$data[2]	";
	foreach my $i(3..$#data){
		my $ctg = (split/,/,$data[$i])[1];
		print OUT "$ctg	";
		}	
	print OUT "\n";
	}
close IN;
close OUT;



sub compare{
	my $ctgT   = shift;
	my $ctgR   = shift;
	my @ctgTdb = split(/\s+/,$ctgT);
	my @ctgRdb = split(/\s+/,$ctgR);
	my %tdb = ();
	my $num_T  = @ctgTdb;
  map {$tdb{$_}++} @ctgTdb;
  my $num_S  = 0;   ###Number of Same contigs
	my $num_D  = 0;   ###Number of Different contigs	  
  foreach my $ctg(@ctgRdb){
  	if(exists($tdb{$ctg})){
  		$num_S++;
  	}else{
  		$num_D++;
  		}
  	} 
 	if($num_S == $num_T){
 		return 1;
 	}else{
 		return 0;
 		}
	}




