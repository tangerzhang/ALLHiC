#!/usr/bin/perl -w

use Getopt::Std;
getopts "b:r:d:";


if ((!defined $opt_b)|| (!defined $opt_r) || (!defined $opt_d) ) {
    die "************************************************************************
    Usage: perl $0 -b mapping.bam -r refSeq.fasta -d main_results/
      -h : help and usage.
      -b : mapping.bam
      -r : reference genome, fasta format
      -d : LACHESIS main_results/
************************************************************************\n";

}



my %seqdb;
my $ctg;
open(IN, $opt_r) or die"";
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


######GET GROUP IDS######
print "a. Getting group ids ...\n";
print "Reading anchored contigs ...\n";
my %anchordb;
my %gidb;
while(my $file = glob "$opt_d/group*ordering"){
	my $gid = $1 if($file=~/group(\d+).ordering/);
        open(my $in, $file) or die"";
	while(<$in>){
		chomp;
		next if(/#/);
		my $ctg = (split/\s+/,$_)[1];
                $ctg    =~ s/_pilon//g;
		$anchordb{$ctg}->{'gid'}  = $gid;
		$anchordb{$ctg}->{'stat'} = "An";
		$gidb{$gid}->{$ctg}       = "An";
		}
	close $in;
	}

my $ufile = "unanchor.signal.txt";
if(!(-e $ufile)){
  system("touch $ufile");
  }

my $num_of_group  = keys %gidb;
print "Number of groups: $num_of_group\n";
print "Reading unanchored contigs ...\n";
open(IN, "unanchor.signal.txt") or die"";
<IN>;
while(<IN>){
	chomp;
	my $i          = $num_of_group + 1;
        my ($ctg,$gid) = (split/\s+/,$_)[0,$i];
            $ctg       =~ s/_pilon//g;
	$gid =~ s/group//g;
	$anchordb{$ctg}->{'gid'}  = $gid;
	$anchordb{$ctg}->{'stat'} = "Un";
	$gidb{$gid}->{$ctg}       = "Un";
	}
close IN;

print "Output group ids ...\n";
foreach my $gid(sort keys %gidb){
	my $outid = "group".$gid.".ids";
	open(my $out, ">$outid") or die"";
	foreach my $ctg (keys %{$gidb{$gid}}){
	  my $len = length $seqdb{$ctg};
		print $out "$ctg	$len\n" if($gidb{$gid}->{$ctg} eq "An");
		print $out "$ctg	$len	recover\n" if($gidb{$gid}->{$ctg} eq "Un");
		}
	close $out;
	}

print "b. Getting CLM files ...\n";

print "Reading and filtering $opt_b file ...\n";
my %tmprdb = (); ###store reads name
my %infordb;     ###store contig pairs with directions: e.g. A+B+,A+B-,A-B+,A-B-
my $count = 0;   ###used for sorting
open(IN, "samtools view $opt_b |awk \'\$7!=\"*\" && \$7!=\"=\"\' |") or die"";
while(<IN>){
	chomp;
	$_                =~ s/_pilon//g;
	my @data          = split(/\s+/,$_);
	next if(exists($tmprdb{$data[0]}));
	$tmprdb{$data[0]}++;
	my ($ctgA,$ctgB)  = sort ($data[2], $data[6]);
###determine gid for the contig pairs
  next if(!exists($anchordb{$ctgA}->{'gid'}));
  next if(!exists($anchordb{$ctgB}->{'gid'}));
  my $ctgAgid       = $anchordb{$ctgA}->{'gid'};
	my $ctgBgid       = $anchordb{$ctgB}->{'gid'};
	next if($ctgAgid ne $ctgBgid);
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
#  print ">$_\n";
#  print "$ctgA length=$ctgAL	and $ctgB length=$ctgBL\n";
#  print "$ctgA+ $ctgB+: $ApBp\n";
#  print "$ctgA+ $ctgB-: $ApBm\n";
#  print "$ctgA- $ctgB+: $AmBp\n";
#  print "$ctgA- $ctgB-: $AmBm\n";
#  print "\n";
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
foreach my $key (sort {$infordb{$a}->{'c'}<=>$infordb{$b}->{'c'}} keys %infordb){
	my @t           = split(/\s+/,$infordb{$key}->{'d'});
	my $num_of_link = @t;
	print ALLCLM "$key	$num_of_link	$infordb{$key}->{'d'}\n";
	}
close ALLCLM;

foreach my $gid (keys %gidb){
 	my $outfile        = "group".$gid.".clm";
	open(my $out, ">$outfile") or die"";
	foreach my $key (sort {$infordb{$a}->{'c'}<=>$infordb{$b}->{c}}  keys %infordb){
		my @t            = split(/\s+/,$infordb{$key}->{'d'});
		my $num_of_link  = @t;
		print $out "$key	$num_of_link	$infordb{$key}->{'d'}\n" if($infordb{$key}->{'g'} eq $gid);
		}
	
	close $out;
	}


