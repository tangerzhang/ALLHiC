#!/usr/bin/perl -w


use Getopt::Std;
getopts "i:b:r:";


if ((!defined $opt_i)|| (!defined $opt_b)|| (!defined $opt_r)) {
    die "************************************************************************
    Usage: perl $0 -i Allele.ctg.table -b bam.list -r draft.asm.fasta
      -h : help and usage.
      -i : Allele.ctg.table 
      -b : bam.list, a file contains input bam files
      -r : draft.sam.fasta
************************************************************************\n";
}

my $bamfile = $opt_b;
my $table   = $opt_i;
my $refSeq  = $opt_r;
### Read bam files

my %pairdb = ();
my %ctgdb  = ();
my %bamdb  = ();
open(IN, $bamfile) or die"";
while(<IN>){
	chomp;
	my $bam = $_;
	   $bam =~ s/\s+//g;
	next if(!($bam =~ /.bam/));
	$bamdb{$bam}++;
	open(my $fh, "samtools view $bam |") or die"";
	while(<$fh>){
		chomp;
		my @data = split(/\s+/,$_);
		my $ctg1 = $data[2];
		my $ctg2 = $data[6];
		next if($ctg2 eq "=");
		my ($sa,$sb) = sort ($ctg1,$ctg2);
		$pairdb{$sa}->{$sb}  .= $data[0].",";
		$ctgdb{$ctg1}++; $ctgdb{$ctg2}++;
		}
	close $fh;
	}
close IN;

### Read allele information
### Remove signal between alleles
open(OUT1, ">removedb_Allele.txt") or die"";
open(LOG, "> log.txt") or die"";
open(IN, $table) or die"";
while(<IN>){
	chomp;
	my @data     = split(/\s+/,$_);
	next if(@data<=3);
	my %tmpdb    = (); ### Record alelle contigs
	my $n        = $#data;
	for(my $i=2;$i<$n;$i++){
		my $ctg1 = $data[$i];
		for(my $j=$i+1;$j<=$n;$j++){
			my $ctg2 = $data[$j];
			my ($sa,$sb) = sort ($ctg1,$ctg2);
			my $key      = $sa.",".$sb;
			$tmpdb{$key}++;
			print OUT1 "$sa	$sb	$pairdb{$sa}->{$sb}\n" if(exists($pairdb{$sa}->{$sb}));
			}
		}
	print LOG ">$_\n";
	foreach my $i(2..$#data){
		my $ctg1    = $data[$i];
		foreach my $ctg2 (keys %ctgdb){
			my ($sa,$sb) = sort ($ctg1,$ctg2);
			my $key      = $sa.",".$sb;
			next if(exists($tmpdb{$key})); 
			next if(!exists($pairdb{$sa}->{$sb}));
			my @rnamedb = split(/,/,$pairdb{$sa}->{$sb});
			my $num_r   = @rnamedb;
			print LOG "$ctg2	$ctg1	$num_r	$pairdb{$sa}->{$sb}\n";
			}
		}
	}
close IN;
close OUT1;
close LOG;

### Remove signal which are not best match with listed alleles (ctgs)
open(OUT2, "> removedb_nonBest.txt") or die"";
open(IN, "log.txt") or die"";
$/='>';
<IN>;
while(<IN>){
	chomp;
	my %hashdb = ();
	my ($name,$info) = split(/\n/,$_,2);
	my @linedb   = split(/\n/,$info);
	foreach my $line(@linedb){
		my @data   = split(/\s+/,$line);
		if(!exists($hashdb{$data[0]})){
		  $hashdb{$data[0]}->{'retain'} = $data[1];
		  $hashdb{$data[0]}->{'num'}    = $data[2];				
		}elsif(exists($hashdb{$data[0]}) and $data[2]>$hashdb{$data[0]}->{'num'}){
		  $hashdb{$data[0]}->{'retain'} = $data[1];
		  $hashdb{$data[0]}->{'num'}    = $data[2];				
	 	}
	 }
	foreach $line (@linedb){
		@data = split(/\s+/,$line);
		if($hashdb{$data[0]}->{'retain'}  eq $data[1]){
#			print OUT2 "$data[0]	$data[1]	$data[2]	retain	$data[3]\n";
      next;
		}else{
			print OUT2 "$data[0]	$data[1]	$data[2]	remove	$data[3]\n";
			}
		}	
	}
close IN;
close OUT2;
system("remove_reads.pl");
### Reading removed reads
#my %removedb = ();
#open(IN, "removedb_Allele.txt") or die"";
#my $content = <IN>;
#my @linedb  = split(/\n/,$content);
#foreach my $line (@linedb){
#	my $info    = (split/\s+/,$line)[2];
#	my @rnamedb = split(/,/,$info);
#	map {$removedb{$_}++} @rnamedb;
#	
#	}
#close IN;
#
#open(IN, "removedb_nonBest.txt") or die"";
#$content = <IN>;
#@linedb  = split(/\n/,$content);
#foreach my $line (@linedb){
#	my $info  = (split/\s+/,$line)[4];
#	my @rnamedb = split(/,/,$info);
#	map {$removedb{$_}++} @rnamedb;	
#	}
#close IN;

#my $num_of_remove_reads = keys %removedb;
#print "Removing $num_of_remove_reads reads\n";

#open(OUT, "> prunning.sam") or die"";
#foreach my $bam (keys %bamdb){
#	open(my $fh, "samtools view $bam |") or die"";
#	$content = <$fh>;
#	@linedb  = split(/\n/,$content);
#	foreach my $line (@linedb){
#		my $rname = (split/\s+/,$line)[0];
#		next if(exists($removedb{$rname}));
#		print OUT "$line\n";
#		}
#	close $fh;
#	}
#close OUT;

system("samtools faidx $refSeq");
my $fai  = $refSeq.".fai";
system("samtools view -bt $fai prunning.sam > prunning.bam");




