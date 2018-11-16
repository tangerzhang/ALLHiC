#!/usr/bin/perl -w

use Getopt::Std;
getopts "g:d:b:r:";


if ((!defined $opt_g)|| (!defined $opt_r)) {
    die "************************************************************************
    Usage: perl $0 -g Allele.gene.table -r draft.asm.fasta
      -h : help and usage.
      -g : Allele.gene.table 
      -b : optional,default prunning.bam
      -r : reference ctg assembly
      -d : optional, default wrk_dir
************************************************************************\n";
}

my $bam    = (defined $opt_b)?$opt_b:"prunning.bam";
my $table  = $opt_g;
my $wrkd   = (defined $opt_d)?$opt_d:"wrk_dir";
my $refSeq = $opt_r;

### Read referece ctg fasta
my %refdb = ();
my $ctgn;
open(IN, $refSeq) or die"";
while(<IN>){
	chomp;
	if(/>/){
		$ctgn = $_;
		$ctgn =~ s/>//g;
		$ctgn =~ s/\s+//g;
	}else{
		$refdb{$ctgn} .= $_;
		}
	}
close IN;

foreach $ctgn (keys %refdb){
	$refdb{$ctgn} =~ s/\s+//g;
	}

### Read prunning BAM file
my %bamdb = ();
my $count = 1;
my %rdb;
open(IN, "samtools view $bam |") or die"";
while(<IN>){
	chomp;
	my $rname = (split/\s+/,$_)[0];
	next if(exists($rdb{$rname}));    ### only retain single-end reads
	$rdb{$rname}++;        
	$bamdb{$count++} = $_;
	}
close IN;
### Assign ctgs to pre-defined clusters

my %ctgdb;
open(IN, $table) or die"";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	my $chrn = $data[1];
	foreach my $i(3..$#data){
		my $ctg = (split/,/,$data[$i])[1];
		$ctgdb{$ctg}->{$chrn}++;
		}
	}
close IN;

my %chrdb; ### pre-defined cluster based on chromosomes of close-releative species
foreach my $ctg (keys %ctgdb){
	my $count = 0;
	foreach my $chrn (sort {$ctgdb{$ctg}->{$b}<=>$ctgdb{$ctg}->{$a}} keys %{$ctgdb{$ctg}}){
		$count++;
		next if($count>1);
#		print "$ctg	$chrn	$ctgdb{$ctg}->{$chrn}\n";
		$chrdb{$chrn} .= $ctg.",";
		}
	}

system("rm -rf $wrkd");
system("mkdir $wrkd");
foreach my $chrn (keys %chrdb){
	next if($chrn=~/tig/);
	next if($chrn=~/ctg/);
	system("rm -rf $wrkd/$chrn");
	system("mkdir $wrkd/$chrn");
	my @ctgdb  = split(/,/,$chrdb{$chrn});
	my %tmpdb = (); $tmpdb{'='}++; ### need retain intra-contig links
### output ctg list	to each cluster
	open(my $out, ">$wrkd/$chrn/ctg.list") or die"";
	map {print $out "$_\n";$tmpdb{$_}++} @ctgdb;
	close $out;
### output ctg sequence to each cluster
	open(my $faout, ">$wrkd/$chrn/seq.fasta") or die"";
	map {chomp;print $faout ">$_\n$refdb{$_}\n" if(exists($refdb{$_}))} @ctgdb;
	close $faout;
### output bam file to each cluster
	open(my $bamout, "> $wrkd/$chrn/sample.clean.sam") or die"";
	foreach my $i(keys %bamdb){
		my ($c1,$c2) = (split/\s+/,$bamdb{$i})[2,6];
		next if(!exists($tmpdb{$c1}) or !exists($tmpdb{$c2}));
		print $bamout "$bamdb{$i}\n";
		}
	close $bamout;
	system("samtools faidx $wrkd/$chrn/seq.fasta");
	system("samtools view -bt $wrkd/$chrn/seq.fasta.fai $wrkd/$chrn/sample.clean.sam > $wrkd/$chrn/sample.clean.bam");
	system("rm $wrkd/$chrn/sample.clean.sam");
	
	}





