#!/usr/bin/perl
###################################
# Adapted from the DO.db package.
# Print out the sqlite queries. This is faster than execute one query at a time.
###################################
use strict;
use warnings;
#use Data::Dumper;

use DBI;

#my $dbh = DBI->connect("dbi:SQLite:dbname=DB.sqlite","","");

my($tablename,$inputdata) = @ARGV;
#$tablename='term';
#$inputdata='term.txt';

##########################
#first delete all
print "BEGIN TRANSACTION;\n";
my $sql="delete  from ".$tablename . ";";
print $sql."\n";
#$dbh->do($sql);


open(IN,$inputdata) or die $!;
my $head=<IN>;
my @a=split "\t",$head;
while(<IN>){
	s/\r|\n//g;
	next unless($_);
	@a=split("\t",$_,-1);
	@a=map {$_ =~ s/\'/\'\'/g; $_} @a;
	my $sth = "INSERT INTO ".$tablename." VALUES ('".join("','", @a)."');";
	print $sth."\n";
}
print "COMMIT;\n";

