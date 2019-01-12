#!/usr/bin/perl -w

die "Usage: perl $0 mapping.bam refSeq.fasta\n" if((!defined $ARGV[0]) or (!defined $ARGV[1]));

my %seqdb;
my $ctg;
open(IN, $ARGV[1]) or die"";
while(<IN>){
	chomp;
	if(/>/){
		$ctg = $_;
		$ctg =~ s/>//g;
		$ctg =~ s/\s+.*//g;
    $ctg =~ s/_pilon//g;
	}else{
		$seqdb{$ctg} .= $_;
		}
	}
close IN;

foreach $ctg(keys %seqdb){
	$seqdb{$ctg} =~ s/\s+//g;
	}
	
print "Reading and filtering $ARGV[0] file ...\n";
my %tmprdb = (); ###store reads name
my %infordb;     ###store contig pairs with directions: e.g. A+B+,A+B-,A-B+,A-B-
my $count = 0;   ###used for sorting
open(IN, "samtools view $ARGV[0] |awk \'\$7!=\"*\" && \$7!=\"=\"\' |") or die"";
while(<IN>){
	chomp;
	$_                =~ s/_pilon//g;
	my @data          = split(/\s+/,$_);
	next if(exists($tmprdb{$data[0]}));
	$tmprdb{$data[0]}++;
	my ($ctgA,$ctgB)  = sort ($data[2], $data[6]);
  my $ctgAL         = length $seqdb{$ctgA};
  my $ctgBL         = length $seqdb{$ctgB};
  my $RAP           = ($data[2] le $data[6])?$data[3]:$data[7];
  my $RBP           = ($data[2] le $data[6])?$data[7]:$data[3]; 
  my $A1            = $RAP;
  my $A2            = $ctgAL - $RAP;
  my $B1            = $RBP;
  my $B2            = $ctgBL - $RBP;
###calculate distance for contig pairs                    
  my $ApBp          = $A2 + $B1;
  my $ApBm          = $A2 + $B2;
  my $AmBp          = $A1 + $B1;
  my $AmBm          = $A1 + $B2;

  my $PApBp         = $ctgA."+ ".$ctgB."+"; #P means pair
  my $PApBm         = $ctgA."+ ".$ctgB."-"; 
  my $PAmBp         = $ctgA."- ".$ctgB."+";
  my $PAmBm         = $ctgA."- ".$ctgB."-";
  $infordb{$PApBp}->{'d'} .= $ApBp." ";    #d means distance
  $infordb{$PApBm}->{'d'} .= $ApBm." ";
  $infordb{$PAmBp}->{'d'} .= $AmBp." ";
  $infordb{$PAmBm}->{'d'} .= $AmBm." ";
  $infordb{$PApBp}->{'g'}  = $ctgAgid ;    #g means group id
  $infordb{$PApBm}->{'g'}  = $ctgAgid ; 
  $infordb{$PAmBp}->{'g'}  = $ctgAgid ; 
  $infordb{$PAmBm}->{'g'}  = $ctgAgid ;   
  $infordb{$PApBp}->{'c'}  = $count++ ;    #c means count
  $infordb{$PApBm}->{'c'}  = $count++ ; 
  $infordb{$PAmBp}->{'c'}  = $count++ ; 
  $infordb{$PAmBm}->{'c'}  = $count++ ;     
	}
close IN;

###Get CLM FILES####
print "Output CLM files ...\n";
open(ALLCLM, "> all.clm") or die"";
print ALLCLM "groupA	groupB	num_of_link	Average_distance	signalDensity	distance_list\n";
foreach my $key (sort {$infordb{$a}->{'c'}<=>$infordb{$b}->{'c'}} keys %infordb){
	my @t           = split(/\s+/,$infordb{$key}->{'d'});
	my $num_of_link = @t;
        my $sum = 0; my $ave = 0;
        map {$sum+=$_} @t;
        $ave   = $sum/$num_of_link;
        $ave   = sprintf("%.2f",$ave);
        my ($g1,$g2) = split(/\s+/,$key);
           $g1      =~ s/[+|-]//g;
           $g2      =~ s/[+|-]//g;
        my $l1      = length $seqdb{$g1};
        my $l2      = length $seqdb{$g2};
        my $len     = $l1 + $l2;
        my $signalD = $num_of_link/$len * 1000;
	print ALLCLM "$key	$num_of_link	$ave	$signalD	$infordb{$key}->{'d'}\n";
	}
close ALLCLM;
