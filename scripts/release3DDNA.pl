#!/usr/bin/perl -w

die "Usage: perl $0 No_of_chr seq.FINAL.fasta\n" if(!defined ($ARGV[0]) or !defined($ARGV[1]));
my $Kchr = $ARGV[0];

open(IN, $ARGV[1]) or die"";
$/='>';
<IN>;
while(<IN>){
	chomp;
	my ($name,$seq) = split(/\n/,$_,2);
	$seq =~ s/\s+//g;
	my $len = length $seq;
	$infordb{$name}->{'seq'} = $seq;
	$infordb{$name}->{'len'} = $len;
	}
close IN;

open(OUT, "> chr.fasta") or die"";
my $count = 0;
foreach my $scaf (sort {$infordb{$b}->{'len'}<=>$infordb{$a}->{'len'}} keys %infordb){
	$count++;
	my $chrname = "";
	if($count<=$Kchr){
		$chrname = 'Chr'.$count;
	}else{
		$chrname = 'scaffold'.$count;
		}
	print OUT ">$chrname\n$infordb{$scaf}->{'seq'}\n";
	}

close OUT;


my $ctgn = 0;
open(OUT, ">tig.HiCcorrected.fasta") or die"";
open(IN, "chr.fasta") or die"";
$/='>';
<IN>;
while(<IN>){
	chomp;
	my ($chrn,$seq) = split(/\n/,$_,2);
	print "Process $chrn\n";	
	$seq            =~ s/N/\n/g;
	my $tour        = "";
	my $ctgname     = "";
	my $otour       = $chrn.".tour";
	my @seqdb = split(/\n/,$seq);
	foreach my $i (0..$#seqdb){
		next if ($seqdb[$i] eq "");
		$ctgn++;
		$ctgn = sprintf("%07d",$ctgn);
		$ctgname = "tig".$ctgn;
		$tour      .= $ctgname."+ ";
		print OUT ">$ctgname\n$seqdb[$i]\n";
		}
	next if($chrn =~ /scaffold/);
	open(my $out, ">$otour") or die"";
	print $out ">$chrn\n$tour\n";
	close $out;
	}
close IN;
close OUT;

system("ALLHiC_build tig.HiCcorrected.fasta");
