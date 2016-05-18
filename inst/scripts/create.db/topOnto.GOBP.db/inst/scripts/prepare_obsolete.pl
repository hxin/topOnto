#!/usr/bin/perl
###################################
# Adapted from the DO.db package.
###################################
use strict;
use warnings;

my $usage=<<USAGE;
perl $0 inputfile
		--inputfile is the obo format file you downloaded from Disease Ontology
		  for example perl $0 HumanDO.obo
USAGE

#check the parameter
if(@ARGV<1){
	print $usage;
	exit(1);
}

my $infile=$ARGV[0];

#read the inputed file, and construct a file which corresponds to the table
#do_obsolete in DO.db
open(IN,$infile) or die $!;
open(O1,">obsolete.txt") or die $!;


print O1 "id\tterm\n";

my $str="";
my $is_finish=0;
my %content;
while(<IN>){
	if($is_finish){last;}
	if(/^\[Term\]/ || /^\[Typedef\]/ || eof){
		if($str){
			#get the id and names
			unless($str=~/^\[Term\]/){$str=$_;next ;}
			my ($id,$name);
			if($str=~/id: (.*)\n/){
				$id=$1;
			}
			if($str=~/name: (.*?)\n/){
				$name=$1;
			}
			
			
			my $flag=0;
			#determine whether have string "is_obsolete: true"
			if($str=~/is_obsolete/){
				$flag=1;
			}
			
			if($flag){
				$name=~s/\t/ /g;
				print O1 join "\t",($id,$name);
				print O1 "\n";
			}
		}
		$str=$_;
		if( /^\[Typedef\]/){
			$is_finish=1;
		}
	}else{
		$str.=$_;
	}
}
close IN;
close O1;

