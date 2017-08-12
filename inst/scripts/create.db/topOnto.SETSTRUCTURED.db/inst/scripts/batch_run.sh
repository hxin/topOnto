#!/bin/bash
###################################
# Adapted from the DO.db package.
###################################
#ls |grep -E 'txt$|R$|sqlite$|sql$' |grep -v 'readme' |xargs rm

[ ! -e ont.obo ] && echo "ont.obo does not exist" && exit;

#to run the add_root.pl
echo "###################################";
echo "To run the perl code add_root.pl";
echo "This will add a 'all' root to the ontology";
perl get_relationship.pl ont.obo
perl add_root.pl
echo "###################################";

[ ! -e my.obo ] && echo "my.obo does not exist" && exit;

###################################
#to run the get_relationship.pl
echo "To run the perl code get_relationship.pl";
echo "And it will generate two files:";
echo "1,child2parent.txt";
echo "2,parent2offspring.txt";
perl get_relationship.pl my.obo
echo "###################################";
echo "";


###################################
#to run the prepare_term.pl
echo "To run the perl code prepare_term.pl";
echo "And it will generate one file:";
echo "1,term.txt";
perl prepare_term.pl my.obo
echo "###################################";
echo "";


###################################
#to run the prepare_synonym.pl
echo "To run the perl code prepare_synonym.pl";
echo "And it will generate one file:";
echo "1,synonym.txt";
perl prepare_synonym.pl my.obo
echo "###################################";
echo "";


###################################
#to run the prepare_obsolete.pl
echo "To run the perl code prepare_obsolete.pl";
echo "And it will generate one file:";
echo "1,obsolete.txt";
perl prepare_obsolete.pl my.obo
echo "###################################";
echo "";



###################################
#to run the prepare_offspring.pl
echo "To run the perl code prepare_offspring.pl";
echo "And it will generate one file:";
echo "1,offspring.txt";
perl prepare_offspring.pl my.obo
echo "###################################";
echo "";



###################################
#to run the prepare_parents.pl
echo "To run the perl code prepare_parents.pl";
echo "And it will generate one file:";
echo "1,parents.txt";
perl prepare_parents.pl my.obo
echo "###################################";
echo "";



###################################
#to run the create.db.pl
echo "To run the perl code create.db.pl";
echo "And it will invoke R and generate the DB.sqlite database file";
echo "meanwhile it will insert the meta data into DB.sqlite"
perl create.db.pl
echo "###################################";
echo "";


####################################
#to create insert data 
echo "To run the perl code insert_data for each table file";
echo "And it will generate one file:";
echo "1,sql.sql";
echo "" >sql.sql
perl insert_data.pl term term.txt >>sql.sql
perl insert_data.pl obsolete obsolete.txt >>sql.sql
perl insert_data.pl offspring offspring.txt >>sql.sql
perl insert_data.pl parents parents.txt >>sql.sql
perl insert_data.pl synonym synonym.txt >>sql.sql
perl create_index.pl >>sql.sql
echo "inserting into db...";
sqlite3 DB.sqlite < sql.sql
echo "###################################";
echo "";




####################################
#to remove the template files
echo "cleaning up...";
ls |grep -E 'txt$|R$|sql$|obo$' |grep -v 'readme' |xargs rm


