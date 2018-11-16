#!/usr/bin/perl -w

###This script was used to parse blast+ result (outfmt 6)
###you can get best hit with parameter -b 1
###or -b 0 to get more results 
###The default coverage and identity are 60%, respectively

use Getopt::Std;
getopts "i:o:b:c:d:q:";


if ((!defined $opt_i)|| (!defined $opt_o)  || (!defined $opt_q)) {
    die "************************************************************************
    Usage: perl $0 -i input -o output -q query.fasta -b 0||1  
      -h : help and usage.
      -q : query file, fasta format
      -i : input file is the result of blast+
      -b : (optioanl, default 1)1 means only output best hit; 0 means get more results
      -d : identity (optional, default is 0.6)
      -c : coverage (optional, defalut is 0.6)
      -o : output
************************************************************************\n";
}

$input         = $opt_i;
$output        = $opt_o;
$BestHit_model = (defined $opt_b) ? $opt_b : 1;
$coverage 	   = (defined $opt_c) ? $opt_c : 0.6;
$identity      = (defined $opt_d) ? $opt_d : 0.6;

open(IN, $opt_q) or die"No query file: $opt_q\n";
while(<IN>){
	if(/>/){
		$gene = $_;
		$gene =~ s/>//g;
		$gene =~ s/\s+.*//g;
	}else{
		$infordb{$gene} .= $_;
		}
	}
close IN;

open(OUT, "> $output") or die"No output file: $output\n";
open(IN, $input) or die"No input file: $input\n";
while(<IN>){
	chomp;
	@data    = split(/\s+/,$_);
	$query   = $data[0];
	$countdb{$query} += 1;
	next if($countdb{$query}>1 and $BestHit_model==1);
	$q_len   = length $infordb{$query};
#	$subject = $data[1];
	$blst_i  = $data[2]/100;
	$blst_c  = ($data[7]-$data[6])/$q_len;
  if($blst_i>=$identity and $blst_c>=$coverage){
  	print OUT "$_\n";
  	}
	}
close IN;
close OUT;





