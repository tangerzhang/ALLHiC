#!/usr/bin/perl -w

die "Usage: perl $0 chr.agp\n" if(!defined $ARGV[0]);
my $agp = $ARGV[0];
my %uctgdb;
my %actgdb;
my %chrdb;
my $sumL = 0;
my $sumC = 0;
my $sumU = 0;
open(IN, "grep -v 'contig' $agp |grep -v '#'|") or die"";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	$sumC++;
	$sumL   += $data[7];
	if($data[0] eq $data[5]){
	  $uctgdb{$data[5]} = "Unanchor"; 
	  $sumU += $data[7];
	}else{
		$actgdb{$data[5]} = "Anchor";
		$chrdb{$data[0]}->{'ctg'}++;
		$chrdb{$data[0]}->{'len'} = $data[2];
		}
	}
close IN;

my $numU = keys %uctgdb;
my $numA = keys %actgdb;
my $sumA = 0;
print "ChrID	Anchored_ctg	Length\n";
foreach my $chrn (sort {$chrdb{$b}->{'len'}<=>$chrdb{$a}->{'len'} } keys %chrdb){
	$sumA += $chrdb{$chrn}->{'len'};
	print "$chrn	$chrdb{$chrn}->{'ctg'}	$chrdb{$chrn}->{'len'}\n";
	}

print "Total number of contigs (bp): $sumC\n";
print "Total length of contigs (bp): $sumL\n";
print "Total number of anchored contgis: $numA\n";
print "Total length of chromosome level assembly (bp): $sumA\n";
print "Number of unanchored contigs: $numU\n";
print "Length of unanchored contigs: $sumU\n";
my $arate = (1-$sumU/$sumL)*100;
$arate = sprintf("%.2f",$arate);
print "Anchor rate (%): $arate\n";
