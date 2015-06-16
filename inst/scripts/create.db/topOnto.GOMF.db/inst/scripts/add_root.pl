#!/usr/bin/perl
###################################
# Adapted from the DO.db package.
###################################
use strict;
use warnings;


my $infile='ont.obo';


open(IN,"child2parent.txt") or die $!;
open(OBO,$infile) or die $!;
open(OBOADDALL,">my.obo") or die $!;

my $header=<IN>;
my (%ps,%cs);
while(<IN>){
	chomp;
	my @item=split('\t',$_);
	$cs{$item[0]}=1;
	$ps{$item[1]}=1;
}
##find root
my @roots;

map {push(@roots,$_) if !defined($cs{$_})} keys %ps;
my %roots=map{$_=>1} @roots;
print "found root terms:";
print "@roots";
print "\n";

if(scalar(@roots)==1 and $roots[0] eq 'all'){
	#do nothing
	while(<OBO>){
		my $line=$_;
		print OBOADDALL $line;
	}
}else{
	##generate new obo file with 'all' as root
	my $flag=0;
	while(<OBO>){
		my $line=$_;
		if($line =~ /^\[Term\]/ and $flag==0){
			$flag=1;
			$line = "[Term]\nid: all\nname: all\n\n".$line;
		}
		my $id;
		if($line=~/id: (.*)\n/){
			$id=$1;
			if(defined($roots{$id})){
				$line .= "is_a: all ! all\n";
			}
		}
		print OBOADDALL $line;
	}
}



close IN;
close OBO;
close OBOADDALL;
