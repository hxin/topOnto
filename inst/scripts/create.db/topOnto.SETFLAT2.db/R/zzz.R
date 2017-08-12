datacache <- new.env(hash=TRUE, parent=emptyenv())

ONT <- function() showQCData("ONT", datacache)
ONT_dbconn <- function() dbconn(datacache)
ONT_dbfile <- function() dbfile(datacache)
##@todo
ONT_dbschema <- function(){writeLines(strwrap(readLines(system.file("DBschemas","schemas_1.0","DB.sql", package ="topOnto.HPO.db")),indent=2, exdent=4))}
ONT_dbInfo <- function() dbInfo(datacache)



.onLoad <- function(libname, pkgname)
{
  #require("methods", quietly=TRUE)
  
  ##detach other topOnto.xx.db package to avoid name conflict
   while(length(grep("topOnto.\\w+.db", search(), perl=TRUE, value=FALSE)) >0 ){
     detach(pos = grep("topOnto.\\w+.db", search(), perl=TRUE, value=FALSE)[1], unload=TRUE,force=TRUE)
   }
  
  setClass("ONTTermsAnnDbBimap", contains="AnnDbBimap")
  
  ## Connect to the SQLite DB
  dbfile <- system.file("extdata", "DB.sqlite", package=pkgname, lib.loc=libname)
  assign("dbfile", dbfile, envir=datacache)
  dbconn <- dbFileConnect(dbfile)
  assign("dbconn", dbconn, envir=datacache)
  
  ## Create the AnnObj instances
  ann_objs<-createAnnObjs.ONT_DB("ONT", "ONT", dbconn, datacache)
  mergeToNamespaceAndExport(ann_objs, pkgname)
  
  #packageStartupMessage(paste(pkgname,"loaded!",sep=' '))
}

.onUnload <- function(libpath)
{
  dbFileDisconnect(ONT_dbconn())
}
