###################################
# The db schema is not in AnnotationDbi package. In order to create the AnnoObjs, we need to tell AnnotationDbi
# how to create the AnnoObjs. This file defines the object.
#
# The script is created on top of the DO.db package.
###############Notes#######################
#This file  contains two main objects: 
#1)DO_DB_AnnDbBimap_seeds;
#2)createAnnObjs.DO_DB
#which is created by taking file 'createAnnObjs.GO_DB.R' in AnnotationDbi package as a reference
#We defined 7 Bimaps(5 are defined in DO_DB_AnnDbBimap_seeds and 2 are defined by revmap2 which defined in createAnnObjs.DO_DB)
#Bimap object TERM, OBSOLETE, SYNONYM are similiar to GO but can't be assign to GOTermsAnnDbBimap class
#and currently assign them to AnnDbBimap class, so when applying method such as 'as.list','toTable' etc. to 
#them, it can not fetch the items defined in Rattribnames. So we may need to define an DOTermsAnnoDbBimap class
#and create functions as.list, toTable etc. to DOTermsAnnoDbBimap class
#
###################################
###################################
### =========================================================================
### An SQLite-based ann data package (AnnDbPkg) provides a set of pre-defined
### AnnObj objects that are created at load-time. This set depends only on
### the underlying db schema i.e. all the SQLite-based ann data packages that
### share the same underlying db schema will provide the same set of AnnObj
### objects.
###
### This file describes the set of AnnObj objects provided by any
### DO_DB-based package i.e. any SQLite-based ann data package based
### on the DO_DB schema.
### The createAnnObjs.DO_DB() function is the main entry point for
### this file: it is called by any DO_DB-based package at load-time.
### -------------------------------------------------------------------------
### Mandatory fields: objName, Class and L2Rchain
ONT_DB_AnnDbBimap_seeds <- list(
  list(
    objName="PARENTS",
    Class="AnnDbBimap",
    L2Rchain=list(
      list(
        tablename="term",
        Lcolname="id",
        Rcolname="_id"
      ),
      list(
        tablename="parents",
        Lcolname="_id",
        tagname=c(RelationshipType="{relationship_type}"),
        Rcolname="_parent_id"
      ),
      list(
        tablename="term",
        Lcolname="_id",
        Rcolname="id"
      )
    )
  ),
  
  list(
    objName="ANCESTOR",
    Class="AnnDbBimap",
    L2Rchain=list(
      list(
        tablename="term",
        Lcolname="id",
        Rcolname="_id"
        #filter="{ontology}='CC'"
      ),
      list(
        tablename="offspring",
        Lcolname="_offspring_id",
        Rcolname="_id"
      ),
      list(
        tablename="term",
        Lcolname="_id",
        Rcolname="id"
      )
    )
  ),
  
  list(
    objName="TERM",
    Class="ONTTermsAnnDbBimap",
    L2Rchain=list(
      list(
        tablename="term",
        Lcolname="id",
        Rcolname="id",
        Rattribnames=c(
          Term="{term}",
          Ontology="{ontology}",
          Definition="{definition}",
          Synonym="synonym.synonym",
          Secondary="synonym.secondary"
        ),
        Rattrib_join="LEFT JOIN synonym ON {_id}=synonym._id"
      )
    )
  ),
  list(
    objName="OBSOLETE",
    Class="ONTTermsAnnDbBimap",
    L2Rchain=list(
      list(
        tablename="obsolete",
        Lcolname="id",
        Rcolname="id",
        Rattribnames=c(
          Term="{term}",
          Ontology="{ontology}",
          Definition="{definition}",
          ## The RSQLite driver crashes on queries like
          ##   SELECT NULL, ... FROM ...
          ## so a temporary workaround is to use
          ##   SELECT '', ... FROM ...
          #Synonym="NULL",
          #Secondary="NULL"
          Synonym="''",
          Secondary="''"
        )
      )
    )
  ),
  list(
    objName="SYNONYM",
    Class="ONTTermsAnnDbBimap",
    L2Rchain=list(
      list(
        tablename="synonym",
        Lcolname="synonym",
        Rcolname="_id",
        filter="{like_term_id}=1"
      ),
      list(
        tablename="term",
        Lcolname="_id",
        Rcolname="id",
        Rattribnames=c(
          Term="{term}",
          Ontology="{ontology}",
          Definition="{definition}",
          Synonym="synonym.synonym",
          Secondary="synonym.secondary"
        ),
        Rattrib_join="LEFT JOIN synonym ON {_id}=synonym._id"
      )
    )
  )
)

createAnnObjs.ONT_DB <- function(prefix, objTarget, dbconn, datacache)
{
  #Now skip here
  #checkDBSCHEMA(dbconn, "DO_DB") 
  
  ## AnnDbBimap objects
  seed0 <- list(
    objTarget=objTarget,
    datacache=datacache
  )
  #ann_objs <- createAnnDbBimaps(DO_DB_AnnDbBimap_seeds, seed0)
  ann_objs <- AnnotationDbi:::createAnnDbBimaps(ONT_DB_AnnDbBimap_seeds, seed0)
  
  ## Reverse maps
  #I am not sure whether it is suitable for diease ontology
  revmap2 <- function(from, to)
  {
    map <- revmap(ann_objs[[from]], objName=to)
    L2Rchain <- map@L2Rchain
    tmp <- L2Rchain[[1]]@filter
    L2Rchain[[1]]@filter <- L2Rchain[[length(L2Rchain)]]@filter
    L2Rchain[[length(L2Rchain)]]@filter <- tmp
    map@L2Rchain <- L2Rchain
    map
  }
  ann_objs$CHILDREN <- revmap2("PARENTS", "CHILDREN")
  ann_objs$OFFSPRING <- revmap2("ANCESTOR", "OFFSPRING")
  
  ## 1 special map that is not an AnnDbBimap object (just a named integer vector)
  #ann_objs$MAPCOUNTS <- createMAPCOUNTS(dbconn, prefix)
  ann_objs$MAPCOUNTS <- AnnotationDbi:::createMAPCOUNTS(dbconn, prefix)
  
  #prefixAnnObjNames(ann_objs, prefix)  
  AnnotationDbi:::prefixAnnObjNames(ann_objs, prefix)
}

# selectDB <- function(ont_file){
#   ## Connect to the SQLite DB
#   dbfile <- system.file("extdata", ont_file, package="topOnto.db", lib.loc="/home/xin/R/x86_64-pc-linux-gnu-library/3.1")
#   assign("dbfile", dbfile, envir=datacache)
#   dbconn <- dbFileConnect(dbfile)
#   assign("dbconn", dbconn, envir=datacache)
#   
#   #later when Bimap object defined in AnnotationDbi package use
#   #ann_objs <- createAnnObjs.SchemaChoice("DO_DB", "DO", "DO", dbconn, datacache)
#   #but currently we use following
#   #Get the objects.R file in inst/scripts directory which defined Bimap object
#   #bimapfile<-system.file("scripts","objects.R", package=pkgname, lib.loc=libname)
#   #source(bimapfile)
#   ann_objs<-createAnnObjs.ONT_DB(ont_prefix, ont_objTarget, dbconn, datacache)
#   mergeToNamespaceAndExport(ann_objs, "topOnto")   
# }

