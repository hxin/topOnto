#!/usr/bin/perl
###################################
# Adapted from the DO.db package.
###################################
use strict;
use warnings;

open(II,"term.txt") or die $!;<II>;
my %id2_id;
while(<II>){
	s/\r|\n//g;
	next unless($_);
	my($_id,$id,$term) = split "\t";
	$id2_id{$id}=$_id;
}

open(OUT,">offspring.txt") or die $!;
print OUT join "\t",("_id","_offspring_id");
print OUT "\n";
open(IN,"parent2offspring.txt") or die $!;<IN>;
while(<IN>){
	s/\r|\n//g;
	next unless($_);
	my($c,$p) = split "\t";
	print OUT join  "\t",($id2_id{$c},$id2_id{$p});
	print OUT "\n";
}

