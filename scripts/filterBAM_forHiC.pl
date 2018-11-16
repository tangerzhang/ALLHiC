#!/usr/bin/perl -w

die "Usage: perl $0 file.bam out.sam\n" if(!defined($ARGV[0]) or !defined($ARGV[1]));
open(OUT, "> $ARGV[1]") or die"";
my ($NM,$XM,$XO,$XG) = 0;
open(IN, "samtools view $ARGV[0] |") or die"";
while(<IN>){
	chomp;
	my $mapq = (split/\s+/,$_)[4];
	if(/NM:i:(\d)/){
		$NM = $1;
		}
	if(/XM:i:(\d)/){
		$XM = $1;
		}
	if(/XO:i:(\d)/){
		$XO = $1;
		}
	if(/XG:i:(\d)/){
		$XG = $1;
		}
        next if($mapq<30);
	next if(!(/XT:A:U/));	
	next if(!(defined $NM) or $NM>5);
	next if(!(defined $XM) or $XM>3);
	next if(!(defined $XO) or $XO>2);
	next if(!(defined $XG) or $XG>2);
	next if(/XA:/);
	print OUT "$_\n";
	}
close IN;
close OUT;


#Tag	Meaning
#NM	Edit distance
#MD	Mismatching positions/bases
#AS	Alignment score
#BC	Barcode sequence
#X0	Number of best hits
#X1	Number of suboptimal hits found by BWA
#XN	Number of ambiguous bases in the referenece
#XM	Number of mismatches in the alignment
#XO	Number of gap opens
#XG	Number of gap extentions
#XT	Type: Unique/Repeat/N/Mate-sw
#XA	Alternative hits; format: (chr,pos,CIGAR,NM;)*
#XS	Suboptimal alignment score
#XF	Support from forward/reverse alignment
#XE	Number of supporting seeds
