#!/usr/bin/perl
###################################
# Adapted from the DO.db package.
###################################
use strict;
use warnings;


=head sample data from Disease Ontology
format-version: 1.2
date: 06:04:2010 14:44
saved-by: laronhughes
auto-generated-by: OBO-Edit 2.1-beta3
default-namespace: disease_ontology
remark: This is an alpha version and is only for experimental implementation.

[Term]
id: DOID:0000109
name: maturation disease
is_obsolete: true

[Term]
id: DOID:0000634
name: body growth disease
is_obsolete: true

[Term]
id: DOID:0050012
name: chikungunya
synonym: "Chikungunya fever" RELATED []
synonym: "Chikungunya virus disease " RELATED []
xref: ICD10:A92.0
is_a: DOID:1329 ! arbovirus infectious disease
=cut

my $usage=<<USAGE;
perl $0 inputfile
		--inputfile is the obo format file you downloaded from Disease Ontology
		  for example perl $0 HumanDO.obo
USAGE

##check the parameter
#if(@ARGV<1){
#	print $usage;
#	exit(1);
#}
#
my $infile=$ARGV[0];
#my $infile='hp.obo';
open(my $fh, '<', 'config') or die "Unable to open file:$!\n"; 
my %config = map { chomp;split(/=/,$_,2); } <$fh>;

#read the inputed file, and construct a file with four columns:
#_id
#id
#term name
#definition

open(IN,$infile) or die $!;
open(O1,">term.txt") or die $!;


print O1 "_id\tid\tterm\tontology\tdefinition\n";

my $str="";
my $is_finish=0;
my %content;
while(<IN>){
	if($is_finish){last;}
	if(/^\[Term\]/ || /^\[Typedef\]/ || eof){
		if($str){
			#get the id and names
			unless($str=~/^\[Term\]/){$str=$_;next ;}
			my ($id,$name,$def);
			if($str=~/id: (.*)\n/){
				$id=$1;
			}
			if($str=~/name: (.*?)\n/){
				$name=$1;
			}
			if($str=~/def: (.*?)\n/){
				$def=$1;
			}else{
				$def='';
			}
			
			my $flag=0;
			#determine whether have string "is_obsolete: true"
			if($str=~/is_obsolete/){
				$flag=1;
			}
			
			unless($flag){
				$content{$id}->{'name'}=$name;
				$content{$id}->{'def'}=$def;
				$content{$id}->{'ontology'}=$config{'SOURCENAME'};
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

my $index=0;
foreach my $id(sort keys %content){
	$index++;
	print O1 join "\t",($index,$id,$content{$id}->{'name'},$content{$id}->{'ontology'},$content{$id}->{'def'});
	print O1 "\n";
}

close IN;
close O1;
