#!/usr/bin/perl
###################################
# Adapted from the DO.db package.
###################################
use strict;
use warnings;
#use PadWalker;

###########################
#use RSQLite to create the DO.sqlite

##readin config
#my($DBSCHEMA,$DBSCHEMAVERSION,$SOURCENAME,$SOURCURL,$SOURCEDATE,$VERSION);
open(my $fh, '<', 'config') or die "Unable to open file:$!\n"; 
my %config = map { chomp;split(/=/,$_,2); } <$fh>;


open(II,"term.txt") or die $!;<II>;
my %id2_id;
while(<II>){
	s/\r|\n//g;
	next unless($_);
	my($_id,$id,$term) = split "\t";
	$id2_id{$id}=$_id;
}

#my $currentdir=`pwd`;
my $rfile="create_DO.db.R";
open(R,">${rfile}") or die $!;
my $schema="";
while(<DATA>){
	$schema.=$_;
}


my $rstr=<<RDOC;
library(RSQLite)
drv<-dbDriver("SQLite")
dbfile="DB.sqlite"
db <- dbConnect(drv, dbname=dbfile)
schema.text<-'$schema'
create.sql <- strsplit(schema.text, "\\n")[[1]]
create.sql <- paste(collapse="\\n", create.sql)
create.sql <- strsplit(create.sql, ";")[[1]]
create.sql <- create.sql[-length(create.sql)] # nothing to run here
tmp <- sapply(create.sql, function(x) sqliteQuickSQL(db, x))

RDOC


###########################################
#add metadata
my $metadata=<<META;
metadata<-rbind(c("DBSCHEMA","$config{'DBSCHEMA'}"),
		c("DBSCHEMAVERSION","$config{'DBSCHEMAVERSION'}"),
		c("SOURCENAME","$config{'SOURCENAME'}"),
		c("SOURCURL","$config{'SOURCURL'}"),
		c("SOURCEDATE","$config{'SOURCEDATE'}"),
		c("VERSION","$config{'VERSION'}"));
q<-paste(sep="","INSERT INTO 'metadata' VALUES('",metadata[,1],"','",metadata[,2],"');")
tmp<-sapply(q,function(x) sqliteQuickSQL(db,x))		
META

$rstr.=$metadata;

	

##################################
#add map_metadata
my $map_metadata=<<MAP;
map_metadata<-rbind(c("TERM","$config{'SOURCENAME'}","$config{'SOURCURL'}","$config{'SOURCEDATE'}"),
		    c("OBSOLETE","$config{'SOURCENAME'}","$config{'SOURCURL'}","$config{'SOURCEDATE'}"),
		    c("CHILDREN","$config{'SOURCENAME'}","$config{'SOURCURL'}","$config{'SOURCEDATE'}"),
		    c("PARENTS","$config{'SOURCENAME'}","$config{'SOURCURL'}","$config{'SOURCEDATE'}"),
		    c("ANCESTOR","$config{'SOURCENAME'}","$config{'SOURCURL'}","$config{'SOURCEDATE'}"),
		    c("OFFSPRING","$config{'SOURCENAME'}","$config{'SOURCURL'}","$config{'SOURCEDATE'}")	
);
q<-paste(sep="","INSERT INTO 'map_metadata' VALUES('",map_metadata[,1],"','",map_metadata[,2],"','",map_metadata[,3],"','",map_metadata[,4],"');")
tmp<-sapply(q,function(x) sqliteQuickSQL(db,x))	

MAP

$rstr.=$map_metadata;

##############################
#data for map_counts
my $term_counts=getCounts("term.txt");
my $obsolete_counts=getCounts("obsolete.txt");
my $children_counts=`sed -n '2,\$p' parents.txt |cut -f2|sort|uniq|wc -l |cut -f 1|tr -d '\n'`;
my $parents_counts=`sed -n '2,\$p' parents.txt |cut -f1|sort|uniq|wc -l |cut -f 1|tr -d '\n'`;
my $ancestor_counts=`sed -n '2,\$p' offspring.txt |cut -f2|sort|uniq|wc -l  |cut -f1 |tr -d '\n'`;
my $offspring_counts=`sed -n '2,\$p' offspring.txt |cut -f1|sort|uniq|wc -l  |cut -f1 |tr -d '\n'`;

my $map_counts=<<MAP;
map_counts<-rbind(c("TERM","$term_counts"),
		c("OBSOLETE","$obsolete_counts"),
		c("CHILDREN","$children_counts"),
		c("PARENTS","$parents_counts"),
		c("ANCESTOR","$ancestor_counts"),
		c("OFFSPRING","$offspring_counts"));
q<-paste(sep="","INSERT INTO 'map_counts' VALUES('",map_counts[,1],"','",map_counts[,2],"');")
tmp<-sapply(q,function(x) sqliteQuickSQL(db,x))			
MAP

$rstr.=$map_counts;

$rstr.="\ndbDisconnect(db)\n";


print R $rstr;

close R;

 `R --no-save < $rfile`;


sub getCounts{
	my($in) = @_;
	my $l=`sed -n '2,\$p' $in |wc -l`;
	my $count=0;
	if($l=~/(\d+)/){
		$count=$1;
	}
	return $count;
}

__DATA__
--
-- DB schema
-- ====================
--
CREATE TABLE term (
  _id INTEGER PRIMARY KEY,
  id VARCHAR(12) NOT NULL UNIQUE,               -- DI ID
  term VARCHAR(255) NOT NULL,                 -- textual label for the DO term
  ontology VARCHAR(9) NOT NULL,
  definition TEXT default NULL
);

CREATE TABLE synonym (
  _id INTEGER NOT NULL,                     -- REFERENCES do_term
  synonym VARCHAR(255) NOT NULL,                -- label or DO ID
  secondary VARCHAR(12) NULL,                      -- DO ID
  like_term_id SMALLINT,                          -- boolean (1 or 0)
  FOREIGN KEY (_id) REFERENCES term (_id)
);


CREATE TABLE parents ( 
  _id INTEGER NOT NULL,                     -- REFERENCES do_term
  _parent_id INTEGER NOT NULL,                   -- REFERENCES do_term
  relationship_type VARCHAR(7) NOT NULL,                 -- type of DO child-parent relationship
  FOREIGN KEY (_id) REFERENCES term (_id),
  FOREIGN KEY (_parent_id) REFERENCES term (_id)
);

CREATE TABLE offspring (
  _id INTEGER NOT NULL,                     -- REFERENCES do_term
  _offspring_id INTEGER NOT NULL,                -- REFERENCES do_term
  FOREIGN KEY (_id) REFERENCES term (_id),
  FOREIGN KEY (_offspring_id) REFERENCES term (_id)
);


CREATE TABLE obsolete (
  id VARCHAR(12) PRIMARY KEY,                   -- DO ID
  term VARCHAR(255) NOT NULL                   -- textual label for the DO term
)
;

CREATE TABLE map_counts (
  map_name VARCHAR(80) PRIMARY KEY,
  count INTEGER NOT NULL
);

CREATE TABLE map_metadata (
  map_name VARCHAR(80) NOT NULL,
  source_name VARCHAR(80) NOT NULL,
  source_url VARCHAR(255) NOT NULL,
  source_date VARCHAR(20) NOT NULL
);

CREATE TABLE metadata (
  name VARCHAR(80) PRIMARY KEY,
  value VARCHAR(255)
);

-- Indexes

