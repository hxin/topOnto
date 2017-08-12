########Readme File##################
# Adapted from the DO.db package.
#####################################
This file will describe how we convert the OBO format file of neccessary files which are needed by topOnto.xx.db


1)put your obo file into this folder
2)edit the ./config file
2)cd to this folder
3)chmod +x batch_run.sh
./batch_run.sh

4)This will create a file 'DB.sqlite' in the current folder.

Note:You should have perl installed and perl packages "DBI" and "DBD::SQLite" are required, and 
you also need to have R installed and R package "RSQLite" is required.


