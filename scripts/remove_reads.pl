#!/usr/bin/perl -w

my %bamdb  = ();
open(IN, "bam.list") or die"";
while(<IN>){
        chomp;
        my $bam = $_;
           $bam =~ s/\s+//g;
        next if(!($bam =~ /.bam/));
        $bamdb{$bam}++;

        }
close IN;


my %removedb = ();
open(IN, "removedb_Allele.txt") or die"";
while(<IN>){
	chomp;
	my $info    = (split/\s+/,$_)[2];
	my @rnamedb = split(/,/,$info);
	map {$removedb{$_}++} @rnamedb;
	}
close IN;

open(IN, "removedb_nonBest.txt") or die"";
while(<IN>){
	chomp;
	my $info  = (split/\s+/,$_)[4];
	my @rnamedb = split(/,/,$info);
	map {$removedb{$_}++} @rnamedb;
	}
close IN;

my $num_of_remove_reads = keys %removedb;
print "Removing $num_of_remove_reads reads\n";


open(OUT, "> prunning.sam") or die"";
foreach my $bam (keys %bamdb){
        open(my $fh, "samtools view $bam|") or die"";
        while(<$fh>){
        	chomp;
                my @data  = split(/\s+/,$_);
        	my $rname = (split/\s+/,$_)[0];
                my $ctg2  = $data[6];
                next if($ctg2 eq "*");
        	print OUT "$_\n" if(!exists($removedb{$rname}));
        	}
        close $fh;
        }
close OUT;



