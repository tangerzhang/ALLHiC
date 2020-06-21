#!/usr/bin/perl -w

while(my $file=glob "*_orderings.txt"){
	my $name = $file; 
		 $name =~ s/_orderings.txt//g;
		 $name .= ".tour";
  open(my $out, "> $name") or die"";
	open(my $fh, $file) or die"";
	while(<$fh>){
		chomp;
		my ($ctg,$dir) = (split/\s+/,$_)[0,1];
		my $line = $ctg."".$dir;
		print $out "$line	";
		}
	close $fh;
	close $out;
	}

