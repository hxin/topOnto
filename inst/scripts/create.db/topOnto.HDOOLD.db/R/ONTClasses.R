setClass("ONTTermsAnnDbBimap",contains="AnnDbBimap");
setClass("ONTTerms",
         representation(
           ID="character",       # a single string (mono-valued)
           Term="character",       # a single string (mono-valued)
           Ontology="character",   # a single string (mono-valued)
           Definition="character", # a single string (mono-valued)
           Synonym="character",    # any length including 0 (multi-valued)
           Secondary="character"   # any length including 0 (multi-valued)
         )
)

### The mono-valued slots are also the mandatory slots.
.TERMNODE_MONOVALUED_SLOTS <- c("ID", "Term", "Ontology", "Definition")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization.
###

setMethod("initialize", "ONTTerms",
          function(.Object, ...)
          {
            args <- list(...)
            argnames <- names(args)
            if (is.null(argnames) || any(argnames == ""))
              stop("all arguments must be named")
            argnames <- match.arg(argnames, slotNames(.Object), several.ok=TRUE)
            if (!(all(.TERMNODE_MONOVALUED_SLOTS %in% argnames))) {
              s <- paste(.TERMNODE_MONOVALUED_SLOTS, collapse=", ")
              stop("arguments ", s, " are mandatory")
            }
            for (i in seq_len(length(args))) {
              argname <- argnames[i]
              value <- args[[i]]
              if ((argname %in% .TERMNODE_MONOVALUED_SLOTS)) {
                if (length(value) != 1)
                  stop("can't assign ", length(value),
                       " values to mono-valued slot ", argname)
              } else {
                value <- value[!(value %in% c(NA, ""))]
              }
              slot(.Object, argname) <- value
            }
            .Object
          }
)

ONTTerms <- function(Id, term, ontology, synonym = "", secondary = "",
                    definition = ""){
  return(new("ONTTerms", ID = Id, Term = term,
             Synonym = synonym, Secondary = secondary,
             Definition = definition, Ontology = ontology))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "ID", "Term", "Ontology", "Definition", "Synonym" and "Secondary" 
### generics (accessor methods).
###

setGeneric("ID", function(object) standardGeneric("ID")) 
setGeneric("Term", function(object) standardGeneric("Term"))
setGeneric("Ontology", function(object) standardGeneric("Ontology"))
setGeneric("Definition", function(object) standardGeneric("Definition"))
setGeneric("Synonym", function(object) standardGeneric("Synonym"))
setGeneric("Secondary", function(object) standardGeneric("Secondary"))

setMethod("ID", "ONTTerms", function(object) standardGeneric("ID")) 
setMethod("Term", "ONTTerms", function(object) standardGeneric("Term")) 
setMethod("Ontology", "ONTTerms", function(object) standardGeneric("Ontology")) 
setMethod("Definition", "ONTTerms", function(object) standardGeneric("Definition")) 
setMethod("Synonym", "ONTTerms", function(object) standardGeneric("Synonym")) 
setMethod("Secondary", "ONTTerms",  function(object) standardGeneric("Secondary")) 


##.id2termField() retrieves ids of type field from term
.id2termField <- function(ids, field){
  ##@todo
  #require("GO.db")
  ##     message(cat("Before SQL \n")) ##test
  sql <- sprintf("SELECT id, %s
                 FROM term
                 WHERE id IN ('%s')",
                 field,
                 paste(ids, collapse="','"))
  res <- dbGetQuery(ONT_dbconn(), sql)
  if(dim(res)[1]==0 && dim(res)[2]==0){
    stop("None of your IDs match IDs.  Are you sure you have valid IDs?")
  }else{
    ans <- res[[2]]
    names(ans) <- res[[1]]
    return(ans[ids]) ##This only works because each GO ID is unique (and therefore a decent index ID)
  }
}

setMethod("ID", "ONTTermsAnnDbBimap",function(object) .id2termField(keys(object),"id") )
setMethod("ID", "character",function(object) .id2termField(object,"id") )

setMethod("Term", "ONTTermsAnnDbBimap",function(object) .id2termField(keys(object),"term") )
setMethod("Term", "character",function(object) .id2termField(object,"term") )

setMethod("Ontology", "ONTTermsAnnDbBimap",function(object) .id2termField(keys(object),"ontology") )
setMethod("Ontology", "character",function(object) .id2termField(object,"ontology") )

setMethod("Definition", "ONTTermsAnnDbBimap",function(object) .id2termField(keys(object),"definition") )
setMethod("Definition", "character",function(object) .id2termField(object,"definition") )


##.id2synonymField() retrieves ids of type field from synonym
.id2synonymField <- function(ids, field){
  ##@todo
  #require("GO.db")
  sql <- paste0("SELECT t.id, s.",field,"
                FROM term AS t, synonym AS s
                WHERE t._id=s._id AND id IN ('",paste(ids, collapse="','"),"')")
  res <- dbGetQuery(ONT_dbconn(), sql)
  if(dim(res)[1]==0 && dim(res)[2]==0){
    stop("None of your IDs match IDs.  Are you sure you have valid IDs?")
  }else{
    ans = split(res[,2],res[,1])
    return(ans[ids])##once again (this time in list context), we are indexing with unique IDs.
  }
}

setMethod("Synonym", "ONTTermsAnnDbBimap",function(object) .id2synonymField(keys(object),"synonym") )
setMethod("Synonym", "character",function(object) .id2synonymField(object,"synonym") )

setMethod("Secondary", "ONTTermsAnnDbBimap",function(object) .id2synonymField(keys(object),"secondary") )
setMethod("Secondary", "character",function(object) .id2synonymField(object,"secondary") )




### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" methods.
###

setMethod("show", "ONTTerms",
          function(object)
          {
            s <- character(0)
            for (slotname in slotNames(object)) {
              x <- slot(object, slotname)
              if ((slotname %in% .TERMNODE_MONOVALUED_SLOTS) && length(x) != 1) {
                warning("mono-valued slot ", slotname,
                        " contains ", length(x), " values")
              } else {
                if (length(x) == 0)
                  next
              }
              s <- c(s, paste0(slotname, ": ", x))
            }
            cat(strwrap(s, exdent=4), sep="\n")
          }
)







#######################################################################
## Now add a convenience method to just represent the GO as a graph.
## This method should take arg for which ontology the user wants, and
## return a graph.
## users who want less can just use subgraph.
## weight of graph edges will always be 1 for each edge

#makeGOGraph <- function(){
#  #require(graph)
#  ##@todo
#  #require(GO.db)
#  df <- toTable(ONTPARENTS)
#  ftM2graphNEL(as.matrix(df[, 1:2]), W=rep(1,dim(df)[1])) 
#}

## f = makeGOGraph("bp")
