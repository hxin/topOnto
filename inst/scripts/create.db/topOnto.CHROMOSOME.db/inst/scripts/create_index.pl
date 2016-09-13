#!/usr/bin/perl
###################################
# Create index for db tables
###################################
use strict;
use warnings;


print createIndex('term','id');
print createIndex('parents','_id');
print createIndex('parents','_parent_id');
print createIndex('offspring','_id');
print createIndex('offspring','_offspring_id');

sub createIndex{
	my ($table,$colume)=@_;
	my $n=int(rand(100));	
	my $sth = "CREATE INDEX ${table}_${n} ON $table ($colume);\n";
	return $sth;
}
