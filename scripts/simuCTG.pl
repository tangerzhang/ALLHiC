#!/usr/bin/perl -w

use Getopt::Std;
getopts "i:m:s:";

if ((!defined $opt_i)|| (!defined $opt_m)  || (!defined $opt_s)) {
    die "************************************************************************
    Usage: perl $0 -i input.fasta -m mean -s SD
      -h : help and usage.
      -i : input.fasta, chromosome assembly
      -m : mean length
      -s : sd
************************************************************************\n";
}else{
  print "************************************************************************\n";
  print "Version 1.1\n";
  print "Copyright to Tanger\n";
  print "RUNNING...\n";
  print "************************************************************************\n";
        
        }

my $mean = lc $opt_m;
my $sd   = lc $opt_s;

if($mean =~ /m/){
	$mean =~ s/m//g;
	$mean = $mean * 1000000;
}elsif($mean =~ /k/){
	$mean =~ s/k//g;
	$mean = $mean * 1000;
	}


if($sd =~ /m/){
	$sd =~ s/m//g;
	$sd = $sd * 1000000;
}elsif($sd =~ /k/){
	$sd =~ s/k//g;
	$sd = $sd * 1000;
	}

print "1. generate a contig assembly with Average length = $mean bp ...\n";

my %chrdb;
open(CTG, "> chrUn.fasta") or die"";
open(OUT, "> new_genome.posi.bed") or die"";
open(IN, $opt_i) or die"";
$/='>';
<IN>;
while(<IN>){
	chomp;
	my ($gene,$seq) = split(/\n/,$_,2);
	$seq =~ s/\s+//g;
	if($gene=~/[C|c]hrUn/){
		print CTG ">$gene\n$seq\n";
		next;
		}
	$chrdb{$gene}   = $seq;
	my $total_len   = length $seq;
	my $num_seq     = int $total_len/$sd + 500;
  system("echo \"data<-rnorm($num_seq,mean=$mean,sd=$sd)\" >>Rscript.txt");
  system("echo \"write.table\(data,file\=\'x.txt\'\) \" >> Rscript.txt");
  system("chmod +x Rscript.txt");
  my $Rcmd = "R CMD BATCH --no-save ./Rscript.txt";
  system($Rcmd);
  my $start = 0; my $l  =  0; my $end = 0;
  open(F, "x.txt") or die"";
  my $content = <F>;
  my @linedb = split(/\n/,$content);
  foreach my $i(1..$#linedb){
  	my $line     = $linedb[$i];
  	$start       = $end+1;
  	$l           = (split/\s+/,$line)[1];
  	$l           = int $l;
#  	next if($l<=0);
    $l           = 0 - $l if($l<0);
  	if($end>$total_len){
  		$end       = $total_len;
  	}else{
  		$end       = $start + $l - 1;
  		}
  	next if($start>=$total_len);
  	print OUT "$gene	$start	$end\n";
  	}
  close F;
  system("rm x.txt");
  system("rm Rscript.*");
	}
close IN;
close OUT;
close CTG;

my $count = 0;
my %tdb;
open(OUT, "> ctg.tmp.fasta") or die"";
open(IN, "new_genome.posi.bed") or die"";
$content = <IN>;
@linedb  = split(/\n/,$content);
foreach $line(@linedb){
	my ($chrn,$a,$b) = split(/\s+/,$line);
	my $L            = $b - $a + 1;
	my $subseq       = substr($chrdb{$chrn},$a-1,$L);
	if(!exists($tdb{$chrn})){
		$count = 0;
		$tdb{$chrn}++;
		$count++;
		$outname       = $chrn.".ctg".$count;
	}else{
		$count++;
		$outname       = $chrn.".ctg".$count;	
		}
	print OUT ">$outname\n$subseq\n";	
	}
close OUT;

system("cat ctg.tmp.fasta chrUn.fasta > ctg.fasta");
system("rm ctg.tmp.fasta");

print "2. get statistics for the contig assembly ...\n";
system("perl ~/software/script/faSize.pl ctg.fasta");

$content = `perl ~/software/script/faSize.pl ctg.fasta`;
my $N50  = $1 if($content=~/N50:\s+(\d+)/);
my $ave  = $1 if($content=~/Average\s+length:\s+(\d+)/);
my $ctgname = "ctg."."n".$N50."_m".$ave.".fasta";
system("mv ctg.fasta ./$ctgname");
