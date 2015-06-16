#' @name annFUN
#' @title Functions which map gene identifiers to GO terms
#' @aliases annFUN.db annFUN annFUN.GO2genes annFUN.gene2GO annFUN.file annFUN.org inverseList readMappings
#' @description 
#'   These functions are used to compile a list of GO terms such that each
#'   element in the list is a character vector containing all the gene
#'   identifiers that are mapped to the respective GO term.  
#' 
#' @usage
#' annFUN.db( feasibleGenes = NULL, affyLib)
#' annFUN.org( feasibleGenes = NULL, mapping, ID = "entrez") 
#' annFUN( feasibleGenes = NULL, affyLib)
#' annFUN.gene2GO( feasibleGenes = NULL, gene2GO)
#' annFUN.GO2genes( feasibleGenes = NULL, GO2genes)
#' annFUN.file( feasibleGenes = NULL, file, ...)
#' 
#' readMappings(file, sep = "\t", IDsep = ",")
#' inverseList(l)
#' 
#'  @param feasibleGenes character vector containing a subset of gene
#'     identifiers. Only these genes will be used to annotate GO
#'     terms. Default value is \code{NULL} which means that there are no
#'     genes filtered.
#' 
#'  @param affyLib character string containing the name of the
#'     Bioconductor annotaion package for a specific microarray chip.
#' 
#'  @param gene2GO named list of character vectors. The list names are
#'     genes identifiers. For each gene the character vector contains the
#'     GO identifiers it maps to. Only the most specific annotations are required.
#'   
#'  @param GO2genes named list of character vectors. The list names are
#'     GO identifiers. For each GO the character vector contains the
#'     genes identifiers which are mapped to it. Only the most specific
#'     annotations are required.
#'   
#'  @param mapping character string specifieng the name of the
#'     Bioconductor package containing the gene mappings for a
#'     specific organism. For example: \code{mapping = "org.Hs.eg.db"}.
#'   
#'  @param ID character string specifing the gene identifier to
#'     use. Currently only the following identifiers can be used:
#'     \code{c("entrez", "genbank", "alias", "ensembl", "symbol",
#'       "genename", "unigene")}
#'   
#'  @param file character string specifing the file containing the annotations.
#' 
#'  @param ... other parameters
#' 
#'  @param sep the character used to separate the columns in the CSV file
#'  
#'  @param IDsep the character used to separate the annotated entities
#' 
#'  @param l a list containing mappings
#' 
#' 
#' @details
#'   All these function restrict the GO terms to the ones belonging
#'   to the specified ontology and to the genes listed in the
#'   \code{feasibleGenes} attribute (if not empty).
#'   
#'   The function \code{annFUN.db} uses the mappings provided
#'   in the Bioconductor annotation data packages. For example, if the
#'   Affymetrix \code{hgu133a} chip it is used, then the user should set
#'   \code{affyLib = "hgu133a.db"}.
#' 
#'   The functions \code{annFUN.gene2GO} and \code{annFUN.GO2genes} are
#'   used when the user provide his own annotations either as a gene-to-GOs
#'   mapping, either as a GO-to-genes mapping.
#'   
#'   The \code{annFUN.org} function is using the mappings from the
#'   "org.XX.XX" annotation packages. The function supports different gene
#'   identifiers.
#' 
#'   The \code{annFUN.file} function will read the annotationsof the type
#'   gene2GO or GO2genes from a text file.          
#' 
#' 
#' 
#' @return
#'   A named(GO identifiers) list of character vectors. 
#' @author Adrian Alexa
#' @seealso
#'   \code{\link{topGOdata-class}}
#' 
#' @examples
#' 
#' library(hgu133a.db)
#' set.seed(111)
#' 
#' ## generate a gene list and the GO annotations
#' selGenes <- sample(ls(hgu133aGO), 50)
#' gene2GO <- lapply(mget(selGenes, envir = hgu133aGO), names)
#' gene2GO[sapply(gene2GO, is.null)] <- NA
#' 
#' ## the annotation for the first three genes
#' gene2GO[1:3]
#' 
#' ## inverting the annotations
#' G2g <- revmap(gene2GO)
#' 
#' ## inverting the annotations and selecting an ontology
#' #go2genes <- annFUN.gene2GO(gene2GO = gene2GO)
#' 
#' 
#' ## generate a GO list with the genes annotations
#' selGO <- sample(ls(hgu133aGO2PROBE), 30)
#' GO2gene <- lapply(mget(selGO, envir = hgu133aGO2PROBE), as.character)
#' 
#' GO2gene[1:3]
#' 
#' ## select only the GO terms for a specific ontology
#' #go2gene <- annFUN.GO2genes(GO2gene = GO2gene)
#' 
#' 
#' ##################################################
#' ## Using the org.XX.xx.db annotations
#' ##################################################
#' \dontrun{
#' ## GO to Symbol mappings (only the BP ontology is used)
#' xx <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
#' head(xx)
#' allGenes <- unique(unlist(xx))
#' myInterestedGenes <- sample(allGenes, 500)
#' geneList <- factor(as.integer(allGenes %in% myInterestedGenes))
#' names(geneList) <- allGenes
#' 
#' GOdata <- new("topGOdata",
#'               ontology = "BP",
#'               allGenes = geneList,
#'               nodeSize = 5,
#'               annot = annFUN.org, 
#'               mapping = "org.Hs.eg.db",
#'               ID = "symbol") 
#' }
#' 
#' @keywords misc
NULL


########################################################
##
## functions related to annotations ...
## 
########################################################

readMappings <- function(file, sep = "\t", IDsep = ",") {
  a <- read.delim(file = file, header = FALSE,
                  quote = "", sep = sep, colClasses = "character")
  
  ## a bit of preprocesing to get it in a nicer form
  map <- a[, 2]
  names(map) <- gsub(" ", "", a[, 1]) ## trim the spaces
  
  ## split the IDs
  return(lapply(map, function(x) gsub(" ", "", strsplit(x, split = IDsep)[[1]])))
}

##not implemented
## function annFUN.org() to work with the "org.XX.eg" annotations
annFUN.org <- function( feasibleGenes = NULL, mapping, ID = "entrez") {
# 
#   tableName <- c("genes", "accessions", "alias", "ensembl",
#                  "gene_info", "gene_info", "unigene")
#   keyName <- c("gene_id", "accessions", "alias_symbol", "ensembl_id",
#                "symbol", "gene_name", "unigene_id")
#   names(tableName) <- names(keyName) <- c("entrez", "genbank", "alias", "ensembl",
#                                           "symbol", "genename", "unigene")
#   
#   
#   ## we add the .db ending if needed 
#   mapping <- paste(sub(".db$", "", mapping), ".db", sep = "")
#   require(mapping, character.only = TRUE) || stop(paste("package", mapping, "is required", sep = " "))
#   mapping <- sub(".db$", "", mapping)
#   
#   geneID <- keyName[tolower(ID)]
#   .sql <- paste("SELECT DISTINCT ", geneID, ", go_id FROM ", tableName[tolower(ID)],
#                 " INNER JOIN ", paste("go", tolower(whichOnto), sep = "_"),
#                 " USING(_id)", sep = "")
#   retVal <- dbGetQuery(get(paste(mapping, "dbconn", sep = "_"))(), .sql)
#   
#   ## restric to the set of feasibleGenes
#   if(!is.null(feasibleGenes))
#     retVal <- retVal[retVal[[geneID]] %in% feasibleGenes, ]
#   
#   ## split the table into a named list of GOs
#   return(split(retVal[[geneID]], retVal[["go_id"]]))
}


##not implemented
annFUN.db <- function( feasibleGenes = NULL, affyLib) {
# 
#   ## we add the .db ending if needed 
#   affyLib <- paste(sub(".db$", "", affyLib), ".db", sep = "")
#   require(affyLib, character.only = TRUE) || stop(paste("package", affyLib, "is required", sep = " "))
#   affyLib <- sub(".db$", "", affyLib)
# 
#   orgFile <- get(paste(get(paste(affyLib, "ORGPKG", sep = "")), "_dbfile", sep = ""))
#   
#   try(dbGetQuery(get(paste(affyLib, "dbconn", sep = "_"))(),
#                  paste("ATTACH '", orgFile(), "' as org;", sep ="")),
#       silent = TRUE)
#   
#   .sql <- paste("SELECT DISTINCT probe_id, go_id FROM probes INNER JOIN ",
#                 "(SELECT * FROM org.genes INNER JOIN org.go_",
#                 tolower(whichOnto)," USING('_id')) USING('gene_id');", sep = "")
# 
#   retVal <- dbGetQuery(get(paste(affyLib, "dbconn", sep = "_"))(), .sql)
# 
#   ## restric to the set of feasibleGenes
#   if(!is.null(feasibleGenes))
#     retVal <- retVal[retVal[["probe_id"]] %in% feasibleGenes, ]
#   
#   ## split the table into a named list of GOs
#   return(split(retVal[["probe_id"]], retVal[["go_id"]]))
}


##not implemented
annFUN <- function( feasibleGenes = NULL, affyLib) {
    
#   require(affyLib, character.only = TRUE) || stop(paste('package', affyLib, 'is required', sep = " "))
#   mapping <- get(paste(affyLib, 'GO2PROBE', sep = ''))
# 
#   if(is.null(feasibleGenes))
#     feasibleGenes <- ls(get(paste(affyLib, 'ACCNUM', sep = '')))
#       
#   ontoGO <- get("Term",envir=.GlobalEnv)
#   goodGO <- intersect(ls(ontoGO), ls(mapping))
# 
#   GOtoAffy <- lapply(mget(goodGO, envir = mapping, ifnotfound = NA),
#                      intersect, feasibleGenes)
#   
#   emptyTerms <- sapply(GOtoAffy, length) == 0
# 
#   return(GOtoAffy[!emptyTerms])
}


## the annotation function
annFUN.gene2GO <- function( feasibleGenes = NULL, gene2GO) {

  ## GO terms annotated to the specified ontology 
  #ontoGO <- get(paste("GO", whichOnto, "Term", sep = ""))
  ontoGO <- get("Term",envir=.GlobalEnv)
  
  ## Restrict the mappings to the feasibleGenes set  
  if(!is.null(feasibleGenes))
    gene2GO <- gene2GO[intersect(names(gene2GO), feasibleGenes)]
  
  ## Throw-up the genes which have no annotation
  if(any(is.na(gene2GO)))
    gene2GO <- gene2GO[!is.na(gene2GO)]
  
  gene2GO <- gene2GO[sapply(gene2GO, length) > 0]
  
  ## Get all the GO and keep a one-to-one mapping with the genes
  allGO <- unlist(gene2GO, use.names = FALSE)
  geneID <- rep(names(gene2GO), sapply(gene2GO, length))

  goodGO <- allGO %in% ls(ontoGO)
  return(split(geneID[goodGO], allGO[goodGO]))
}


## the annotation function
annFUN.GO2genes <- function( feasibleGenes = NULL, GO2genes) {
  
  ## GO terms annotated to the specified ontology 
  #ontoGO <- get(paste("GO", whichOnto, "Term", sep = ""))
  ontoGO <- get("Term",envir=.GlobalEnv)

  ## Select only the GO's form the specified ontology
  GO2genes <- GO2genes[intersect(ls(ontoGO), names(GO2genes))]

  ## Restrict the mappings to the feasibleGenes set  
  if(!is.null(feasibleGenes))
    GO2genes <- lapply(GO2genes, intersect, feasibleGenes)

  return(GO2genes[sapply(GO2genes, length) > 0])
}


## annotation function to read the mappings from a file
annFUN.file <- function( feasibleGenes = NULL, file, ...) {

  ## read the mappings from the file
  mapping <- readMappings(file = file, ...)

  ## we check which direction the mapping is ...
  GOs <- ls(Term)
  if(any(GOs %in% names(mapping)))
    return(annFUN.GO2genes( feasibleGenes, mapping))
  
  return(annFUN.gene2GO( feasibleGenes, mapping))
}  


##not implemented
## function returning all genes that can be used for analysis 
feasibleGenes.db <- function(affyLib) {

#   affyLib <- sub(".db$", "", affyLib)
#   mapping <- get(paste(affyLib, 'GO2PROBE', sep = ''))
#   
#   #ontoGO <- get(paste('GO', whichOnto, "Term", sep = ''))
#   ontoGO <- get("Term",envir=.GlobalEnv)
#   goodGO <- intersect(ls(ontoGO), ls(mapping))
# 
#   return(unique(unlist(mget(goodGO, envir = mapping, ifnotfound = NA))))
}

