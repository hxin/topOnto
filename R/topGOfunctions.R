
#################### misc functions ####################
## 
##
##
##
##
########################################################

## Function that split GOTERM in different ontologies.
## Every new environment contain only the terms from one
## of the ontologies 'BP', 'CC', 'MF'
## @Not using!!
groupGOTerms <- function(where) {
#   if(missing(where))
#     where <- .GlobalEnv
#   where <- as.environment(where)
#  
# #   e <- new.env(hash = T, parent = emptyenv())
# #   assign(paste("GO", "CC", "Term", sep = ""), e, envir = where)
# #   
#   #require('topOnto.db') || stop('package topOnto.db is required  ')
# 
#   sql <- "SELECT term_id FROM term WHERE ontology IN"
#   for(onto in c("CC")) {
#     xx <- dbGetQuery(ONT_dbconn(), paste(sql, "('", onto, "');", sep = ""))$term_id
#     e <- new.env(hash = T, parent = emptyenv())
#     multiassign(xx, value = rep(TRUE, length(xx)), envir = e)
#     assign(paste("GO", onto, "Term", sep = ""), e, envir = where)
#   }
# 
#   cat("\ngroupTerms: \tGOCCTerm environments built.\n")
}

#' @title init methods
#' @description functions to load the ontology objects 
#' @aliases init initWHAT initONT
#' @usage 
#' initWHAT()
#' initONT(ontology)
#' @param ontology the ontology to be loaded. This can be check by running initWHAT()
#' @details initONT loads the ontology object from topOnto.<ontolog>.db packages.
#' @examples
#' initWHAT()
#' \dontrun{
#'  require(topOnto.HDO.db)
#'  initONT('HDO')
#' }
#' @author Xin He
init <- function(){
  cat("topOnto has been loaded. Now you need to specify which ontology you want to use.\n")
  cat("run topOnto::intiWHAT() to find out what ontologies are current supportted.\n")
  cat("run topOnto::intiONT(ontology_name) to choose ontology.\n")
  ##cat("Then follow topGO menu for all the analysis.\n")
  initWHAT()
  #cat("ontology_name can be one of the following:\n")
}

initWHAT <- function(){
  cat("These ontology package(s) are currently available in your libPath:\n")
  x=list.files(path=.libPaths(),pattern='topOnto.*.db')
  cat(sapply(x,function(y){ y=sub('topOnto.','',y,perl=TRUE,);y=sub('.db','',y,perl=TRUE,)},USE.NAMES = FALSE),"\n") 
  #x=list.files(system.file("extdata", package ="topOnto"),pattern='*.sqlite')
  #
}

initONT <- function(ontology='HDO'){
  #require("methods", quietly=TRUE)
  pkg <- paste('topOnto',ontology,'db',sep='.')
  
  ##detach other topOnto.xx.db package to avoid name conflict
  while(length(grep("topOnto.\\w+.db", search(), perl=TRUE, value=FALSE)) >0 ){
    detach(pos = grep("topOnto.\\w+.db", search(), perl=TRUE, value=FALSE)[1], unload=TRUE,force=TRUE)
  }
  
  require(pkg, character.only = TRUE) || stop(paste('package ',pkg,' is required',sep=''))
  #detach(pos = match(paste("package", pkg, sep = ":"), search()))

  
  sql <- "SELECT id FROM term;"
  xx <- dbGetQuery(ONT_dbconn(), sql)$id
  e <- new.env(hash = TRUE, parent = emptyenv())
  multiassign(xx, value = rep(TRUE, length(xx)), envir = e)
  assign("Term", e, envir = .GlobalEnv)
}



## given a list of vectors this function is returning the reverse list
## given a mapping form genes to GO terms as a list, compute which are
## the genes mapped to each GO.
inverseList <- function(l) {
  rId <- unlist(l, use.names = FALSE)
  lId <- rep(names(l), sapply(l, length))
  
  return(split(lId, rId))
}

## A test function for HPO
.samepleRunHDO<-function(){
  cat("Loading HDO objects from db...\n")
  topOnto::initONT('HPO')
  a<-system.file("extdata/annotation","human_g2d_omim", package ="topOnto")
  g<-system.file("extdata/genelist","testList", package ="topOnto")
  geneID2GO <- readMappings(file = a)
  geneNames=names(geneID2GO)
  myInterestingGenes=(read.csv(header = FALSE, file = g))$V1
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames
  GOdata <- new("topGOdata", ontology = "HDO", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
  resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  resultElimFis<- runTest(GOdata, algorithm = "elim", statistic = "fisher")
  allRes <- GenTable(GOdata,  elim = resultElimFis,elim = resultElimFis,topNodes = 30,useLevels=TRUE)
  print(allRes)
  cat("Demo done..Seems working!\n")
}



######################################################################
## function to print all the genes annotated to the specified GOs

if(!isGeneric("printGenes"))
  setGeneric("printGenes", function(object, whichTerms, file, ...) standardGeneric("printGenes"))


########
## TODO
##
## remove the chip argument and replace it with a function which,
## given a list of genes will provide information about these genes
##
## TODO
########


## if the file argument is missing the function will just return
## a list of data.frames, each data.frame containg the gene information for specified GO term
setMethod("printGenes",
          signature(object = "topGOdata", whichTerms = "character", file = "missing"),
          function(object, whichTerms, chip, numChar = 100, simplify = TRUE,
                   geneCutOff = 50, pvalCutOff) { 
            
            term.genes <- genesInTerm(object, whichTerms)
            all.genes <- character()
            all.genes <- unique(unlist(term.genes))
            #lapply(term.genes, function(x) all.genes <<- union(x, all.genes))

            chip <- sub(".db$", "", chip)
            
            LL.lib <- get(paste(chip, "ENTREZID", sep = ""))
            Sym.lib <- get(paste(chip, "SYMBOL", sep = ""))
            GNAME.lib <- get(paste(chip, "GENENAME", sep = ""))

            probeMapping <- data.frame(LL.id = as.integer(unlist(mget(all.genes,
                                         envir = LL.lib, ifnotfound = NA))),
                                       Symbol.id = unlist(mget(all.genes,
                                         envir = Sym.lib, ifnotfound = NA)))

            retList <- vector("list", length(whichTerms))
            names(retList) <- whichTerms
            
            for(gt in whichTerms) {
              affID <- term.genes[[gt]]

              pval <- sort(geneScore(object, affID, use.names = TRUE))
              if(missing(pvalCutOff))
                affID <- names(pval)
              else {
                pval <- pval[pval <= pvalCutOff]
                affID <- names(pval)
              }

              ## we restrict the output to the number of genes
              length(affID) <- min(length(affID), geneCutOff)
              
              ## if there are no genes, there is nothing to print
              if(length(affID) == 0) {
                cat("\n\t No genes over the cutOff for:", gt, "\n")
                next
              }
              
              genesNames <- sapply(mget(affID, envir = GNAME.lib),
                                   function(x) return(x[1]))
              genesNames <- paste(substr(genesNames, 1, numChar),
                                  ifelse(nchar(genesNames) > numChar, "...", ""), sep = "")
              
              retList[[gt]] <- cbind("Chip ID" = affID, probeMapping[affID, ], "Gene name" = genesNames,
                                     "raw p-value" = format.pval(pval[affID], dig = 3, eps = 1e-30))
            }            

            ## if we have only one GO term we return just the data.frame
            if(simplify && length(whichTerms) == 1)
              return(retList[[1]])

            return(retList)
          })


setMethod("printGenes",
          signature(object = "topGOdata", whichTerms = "character", file = "character"),
          ## numChar = "integer"),
          function(object, whichTerms, file, oneFile = FALSE, ...) { 

            infoList <- printGenes(object, whichTerms, ...)

            if(length(whichTerms) == 1) {
              write.table(infoList, quote = TRUE, sep = ",", append = FALSE,
                          file = paste(paste(file, sub(":", "_", whichTerms), sep = "_"), "csv", sep = "."),
                          col.names = colnames(infoList), row.names = 1:nrow(infoList))
              return()
            }
            
            if(oneFile) {
              fAppend <- FALSE
              for(gt in whichTerms) {
                write.table(cbind(GOID = rep(gt, nrow(infoList[[gt]])), infoList[[gt]]),
                            quote = TRUE, sep = ",", append = fAppend,
                            file = paste(file, "csv", sep = "."),
                            col.names = colnames(infoList[[gt]]), row.names = FALSE)
                fAppend <- TRUE
                return()
              }
            }
            
            for(gt in whichTerms) 
              write.table(infoList[[gt]], quote = TRUE, sep = ",", append = FALSE,
                          file = paste(paste(file, sub(":", "_", gt), sep = "_"), "csv", sep = "."),
                          col.names = colnames(infoList[[gt]]), row.names = 1:nrow(infoList[[gt]]))
          })

 ######################################################################


.getTermsDefinition <- function(whichTerms, ontology, numChar = 20, multipLines = FALSE) {
  
  termsNames=Term(ONTTERM)[whichTerms]
  
  if(!multipLines) 
    shortNames <- paste(substr(termsNames, 1, numChar),
                        ifelse(nchar(termsNames) > numChar, '...', ''), sep = '')
  else
    shortNames <- sapply(termsNames,
                         function(x) {
                           a <- strwrap(x, numChar)
                           return(paste(a, sep = "", collapse = "\\\n"))
                         })
  
  names(shortNames) <- names(termsNames)
  return(shortNames[whichTerms])
}


## methodsSig contains a named vector of p-values for each run method
.sigAllMethods <- function (methodsSig) 
{
  names.index <- names(methodsSig[[1]])
  retval <- as.data.frame(lapply(methodsSig, function(x) x[names.index]))
  names(retval) <- names(methodsSig)
  return(retval)
}



######################################################################
if(!isGeneric("GenTable"))
  setGeneric("GenTable", function(object, ...) standardGeneric("GenTable"))

#' @title Diagnostic functions for topGOdata and topGOresult objects.
#' @aliases printGenes-methods printGenes printGenes,topGOdata,character,character-method printGenes,topGOdata,character,missing-method GenTable GenTable,topGOdata-method showGroupDensity
#' @rdname diagnosticMethods
#' @name diagnosticMethods
#' @description The \code{GenTable} function generates a summary of the
#'   results of the enrichment analysis.
#'   
#'   The \code{showGroupDensity} function plots the distributions of the
#'   gene' scores/ranks inside a GO term.
#' 
#'   The \code{printGenes} function shows a short summary of the top genes annotated
#'   to the specified GO terms.
#' 
#' 
#' 
#' @usage
#' GenTable(object, ...)
#' 
#' showGroupDensity(object, whichGO, ranks = FALSE, rm.one = TRUE) 
#' 
#' printGenes(object, whichTerms, file, ...)
#' 
#'   @param object an object of class \code{topGOdata}.
#'   @param whichGO the GO terms for which the plot should be generated.
#'   @param ranks if ranks should be used instead of scores.
#'   @param rm.one the p-values which are 1 are removed. 
#'   @param whichTerms character vector listing the GO terms for which the summary should be printed.
#'   @param file character string specifying the file in which the results should be printed.
#'   @param \dots other
#'        \code{topNodes}  the number of top GO terms to be included in the table / the gene description is trimmed such that it has \code{numChar} characters.
#'     %%GenTable(object, ..., orderBy = 1, ranksOf = 2, topNodes = 10, numChar = 40)
#'     Extra arguments for \code{GenTable} can be:
#'       \code{orderBy} if more than one \code{topGOresult} object is given then \code{orderBy} gives the index of which scores will be used to order the resulting table. Can be an integer index  or a character vector given the name of the \code{topGOresult} object.
#'       \code{ranksOf} same as \code{orderBy} argument except that this parameter shows the relative ranks of the specified result.   
#'       \code{numChar} the GO term definition will be truncated such that only the first \code{numChar} characters are shown.
#'     %% printGenes(object, whichTerms, chip, numChar = 100, simplify = TRUE, geneCufOff = 50, pvalCutOff)
#'     %% printGenes(object, whichTerms, file, oneFile = FALSE, ...)
#'     Extra arguments for \code{printGenes} can be:
#'        \code{chip} character string containing the name of the Bioconductor annotation package for a microarray chip.
#'        \code{numChar} 
#'        \code{simplify} logical variable affecting how the results are returned.
#'        \code{geneCutOff} the maximal number of genes shown for each term.
#'        \code{pvalCutOff} only the genes with a p-value less than \code{pvalCutOff} are shown.
#'        \code{oneFile} if \code{TRUE} then a file for each GO term is generated.
#'   
#' @details
#' 
#'   \code{GenTable} is an easy to use function for summarising the most
#'   significant GO terms and the corresponding p-values. The function
#'   dispatches for \code{topGOdata} and \code{topGOresult} objects, and
#'   it can take an arbitrary number of the later, making comparison
#'   between various results easier.
#'   
#'   Note: One needs to type the complete attribute names (the exact name)
#'   of this function, like: \code{topNodes = 5}, \code{rankOf = "resultFis"}, etc. 
#'   This being the price paid for flexibility of specifying different
#'   number of \code{topGOdata} objects.
#'   
#' 
#'   
#'   The \code{showGroupDensity} function analyse the distribution of the
#'   gene-wise scores for a specified GO term.
#'   The function will show the distribution of the genes in a GO term
#'   compared with the complementary set, using a lattice plot.
#'   
#' 
#'   \code{printGenes} 
#'   The function will generate a table with all the probes annotated to
#'   the specified GO term. Various type of identifiers, the gene name and
#'   the gene-wise statistics are provided in the table. 
#'   
#'   One or more GO identifiers can be given to the function using the
#'   \code{whichTerms} argument. When more than one GO is specified, the
#'   function returns a list of \code{data.frames}, otherwise only one
#'   \code{data.frame} is returned.
#'   
#'   The function has a argument \code{file} which, when specified, will
#'   save the results into a file using the CSV format.
#' 
#'   For the moment the function will work only when the chip used has an
#'   annotation package available in Bioconductor. It will not work with
#'   other type of custom annotations.
#' 
#' 
#' @return A data.frame or a list of data.fames.
#' 
#' @author Adrian Alexa
#' 
#' @seealso
#'   \code{\link{groupStats-class}},
#'   \code{\link{getSigGroups-methods}}
#'

#' @examples
#' 
#' 
#' ########################################
#' ## GenTable
#' ########################################
#' 
#' ## load two topGOresult sample objects: resultFis and resultElimFis
#' data(ONTdata)
#' \dontrun{
#' ##These code needs the topOnto.xx.db package to run
#' ## generate the result of Fisher's exact test
#' sig.tab <- GenTable(GOdata, Fis = resultFis, topNodes = 20)
#' 
#' ## results of both test
#' sig.tab <- GenTable(GOdata, resultFis, resultElimFis, topNodes = 20)
#' 
#' ## results of both test with specified names
#' sig.tab <- GenTable(GOdata, Fis = resultFis, ElimFis = resultElimFis, topNodes = 20)
#' 
#' ## results of both test with specified names and specified ordering
#' sig.tab <- GenTable(GOdata, Fis = resultFis, ElimFis = resultElimFis, orderBy = "ElimFis", ranksOf = "Fis", topNodes = 20)
#' }
#' 
#' ########################################
#' ## showGroupDensity
#' ########################################
#' 
#' TERMID <- "DOID:10652"
#' print(showGroupDensity(GOdata, TERMID, ranks = FALSE, rm.one = FALSE))
#' 
#' 
#' ########################################
#' ## printGenes
#' ########################################
#' 
#' \dontrun{
#' library(hgu95av2.db)
#' goID <- "GO:0006629"
#' 
#' gt <- printGenes(GOdata, whichTerms = goID, chip = "hgu95av2.db", numChar = 40)
#' 
#' goIDs <- c("GO:0006629", "GO:0007076")
#' gt <- printGenes(GOdata, whichTerms = goIDs, chip = "hgu95av2.db", pvalCutOff = 0.01)
#' 
#' gt[goIDs[1]]
#' }
#' 
#' 
#' @keywords methods

setMethod("GenTable",
          signature(object = "topGOdata"),
          ## ... = list of topGOresult object
          ## orderBy = "ANY", ## integer or character (index/name)
          ## ranksOf = "ANY", ## which ranks to be computed (integer/character)
          ## topNodes = "integer",
          ## numChar = "integer",
          ## useLevels = "logical"),
          function(object, ..., orderBy = 1, ranksOf = 2,
                   topNodes = 10, numChar = 40,
                   format.FUN = format.pval, decreasing = FALSE,
                   useLevels = FALSE,cutoff=NULL,show.gene=FALSE) {
            
            resList <- list(...)
            
            ## first for the class of the elements in the list
            if(!all(sapply(resList, is, "topGOresult")))
              stop("Use: topGOdata, topGOresult_1, topGOresult_2, ..., \"parameters\".")
            
            ## if no names were provided we name them
            if(is.null(names(resList)))
              names(resList) <- paste("result", 1:length(resList), sep = "")

            ## obtain the score from the objects
            resList <- lapply(resList, score)
            
            ## order the scores and take care of the case in which only one result is provided
            ## in such case the orderBy and ranksOf parameters are ignored.
            if(length(resList) == 1) {
              orderBy <- ranksOf <- 1
              l <- data.frame(resList)
              names(l) <- ifelse(is.null(names(resList)), "", names(resList)) 
            } else {
              l <- .sigAllMethods(resList)
            }

            index <- order(l[, orderBy], decreasing = decreasing)
            l <- l[index, , drop = FALSE]

            if(decreasing)
              rr <- rank(-l[, ranksOf], ties = "first")
            else
              rr <- rank(l[, ranksOf], ties = "first")
            
            if(length(rownames(l)) < topNodes)
              topNodes = length(rownames(l))
                       
            whichTerms <- rownames(l)[1:topNodes]

            l <- l[whichTerms, , drop = FALSE]
            
            ##apply cutoff
            if(!is.null(cutoff)){
              cut<-sum(l[,orderBy]<=cutoff)
              if(cut<topNodes)
                topNodes=cut
              
              whichTerms <- rownames(l)[1:topNodes]           
              l <- l[whichTerms, , drop = FALSE]
            }
           
            
            rr <- as.integer(rr[1:topNodes])

            shortNames <- .getTermsDefinition(whichTerms, ontology(object), numChar = numChar)
            
            infoMat <- data.frame('TERM ID' = whichTerms, 'Term' = shortNames, stringsAsFactors = FALSE)
            
            ## put the levels of the GO
            if(useLevels) {
              nodeLevel <- buildLevels(graph(object), leafs2root = TRUE)
              nodeLevel <- unlist(mget(whichTerms, envir = nodeLevel$nodes2level))
              infoMat <- data.frame(infoMat, Level = as.integer(nodeLevel))
            }
            
            annoStat <- termStat(object, whichTerms)

            ## if orderBy == ranksOf then there is no need to put the ranks
            if(ranksOf != orderBy) {
              dim(rr) <- c(length(rr), 1)
              colnames(rr) <- paste("Rank in ", ifelse(is.character(ranksOf), ranksOf, colnames(l)[ranksOf]), sep = "")

              infoMat <- data.frame(infoMat, annoStat, rr,
                                    apply(l, 2, format.FUN, dig = 2, eps = 1e-30),
                                    check.names = FALSE, stringsAsFactors = FALSE)

            } else {
              infoMat <- data.frame(infoMat, annoStat,
                                    apply(l, 2, format.FUN, dig = 2, eps = 1e-30),
                                    check.names = FALSE, stringsAsFactors = FALSE)
            }
            
            ##rownames(infoMat) <- whichTerms
            rownames(infoMat) <- 1:length(whichTerms)
            
            if(show.gene){
              sig.genes<-.get.sig.gene.in.term(GOdata,myInterestingGenes)
              hits<-sapply(sig.genes[infoMat$TERM.ID],paste,collapse=',')
              infoMat<-data.frame(infoMat,sig.genes=hits)
            }
            
            return(infoMat)            
          })


## if(!isGeneric("genLatexTable"))
##   setGeneric("genLatexTable", function(object, resList, ...) standardGeneric("genLatexTable"))

## setMethod("genLatexTable",
##           signature(object = "topGOdata",
##                     resList = "list",
##                     orderBy = "ANY", ## integer or character (index/name)
##                     ranksOf = "ANY", ## which ranks to be computed (integer/character)
##                     topNodes = "integer",
##                     numChar = "integer"),
##           function(object, resList, orderBy, ranksOf, topNodes = 10, no.char = 40) {
            

## = infoMat, 
##                         pval.tab = xtable.matrix(infoMat, caption = 'GO terms p-value', label = 'tab:GOinfo')))
          



################################################################################
## function that computes the raw/adjusted p-values based on
## multtest package
#' @name getPvalues
#' @aliases getPvalues
#' @title Convenient function to compute p-values from a gene expression matrix.
#' @description Warping function of "mt.teststat", for computing p-values of a gene expression matrix.
#' @usage
#'    getPvalues(edata, classlabel, test = "t", alternative = c("greater", "two.sided", "less")[1],
#'    genesID = NULL, correction = c("none", "Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD",
#'    "BH", "BY")[8]) 
#'    
#'  @param edata Gene expression matrix.
#'  @param classlabel The phenotype of the data
#'  @param test Which test statistic to use
#'  @param alternative The alternative of the test statistic
#'  @param genesID if a subset of genes is provided
#'  @param correction Multiple testing correction procedure
#' 
#' @return
#'   An named numeric vector of p-values.
#' 
#' 
#' @author Adrian Alexa
#' @seealso
#'   \code{\link{GOKSTest}},  \code{\link{groupStats-class}},
#'   \code{\link{getSigGroups-methods}}
#' 
#' @examples
#' 
#' library(ALL)
#' data(ALL)
#' 
#' ## discriminate B-cell from T-cell
#' classLabel <- as.integer(sapply(ALL$BT, function(x) return(substr(x, 1, 1) == 'T')))
#' 
#' ## Differentially expressed genes
#' geneList <- getPvalues(exprs(ALL), classlabel = classLabel,
#'                        alternative = "greater", correction = "BY")
#' 
#' hist(geneList, 50)
#' 
#' @keywords graphs
getPvalues <- function(edata, classlabel, test = "t",
                       alternative = c("greater", "two.sided", "less")[1],
                       genesID = NULL,
                       correction = c("none", "Bonferroni", "Holm", "Hochberg",
                         "SidakSS", "SidakSD", "BH", "BY")[8]) {
  
  #require('multtest') || stop('package multtest is required')

  ## restrict the dataset
  if(!is.null(genesID))
  edata <- edata[intersect(genesID, rownames(edata)), ]
  genesID <- rownames(edata)
  
  t.stats <- multtest::mt.teststat(edata, classlabel = classlabel, test = test)

  if (alternative == "less") 
    p.values <- pt(t.stats, df = length(classlabel) - 2)
  else if (alternative == "greater") 
    p.values <- pt(t.stats, df = length(classlabel) - 2, lower.tail = FALSE)
  else 
    p.values <- 2 * pt(-abs(t.stats), df = length(classlabel) - 2)
  
  if(correction != "none") {
    p.values <- multtest::mt.rawp2adjp(p.values, correction)
    p.values <- p.values$adjp[order(p.values$index), correction]
  }
  names(p.values) <- genesID
  
  return(p.values)
}

################################################################################

## a <= b 
.sigRatio.ratio <- function(a, b, tolerance = 1e-50) {
  
  ## if a and b are almost equal we return 2
  if(identical(all.equal(a, b, tolerance = tolerance), TRUE))
    return(2)
  
  ## we want to compute b / a, thus we must take care for a = 0
  if(identical(all.equal(a, 0, tolerance = tolerance), TRUE))
    return(abs(b / (a + tolerance)))
  
  return(abs(b / a))
}

.sigRatio.log <- function(a, b, tolerance = 1e-50) {
  
  ## if a and b are almost equal we return 2
  if(identical(all.equal(a, b, tolerance = tolerance), TRUE))
    return(2)
  
  ## we want to compute log(a) / log(b), thus we must take care for b = 1
  if(identical(all.equal(log10(b), 0, tolerance = tolerance), TRUE))
    return(abs(log10(a) / tolerance))
  
  return(abs(log10(a) / log10(b)))
}

## a <= b 
.sigRatio.01 <- function(a, b, tolerance = 1e-50) {
  
  ## if a and b are almost equal we return 2
  if(identical(all.equal(a, b, tolerance = tolerance), TRUE))
    return(2)
  
  if(a < b)
    return(1e50)
  
  return(2)
}

################################################################################
################################################################################

## functions to get gene stats from the topGOdata object
.getGeneData <- function(object) {
  return(c(Annotated = numGenes(object),
           Significant = numSigGenes(object),
           NodeSize = object@nodeSize))
}


## functions to get gene stats from the topGOdata object
.printGeneData <- function(x) {
  cat("Annotation data:\n")
  if("Annotated" %in% names(x))
    cat("    Annotated genes:", x["Annotated"], "\n")
  if("Significant" %in% names(x))
    cat("    Significant genes:", x["Significant"], "\n")
  cat("    Min. no. of genes annotated to a term:", x["NodeSize"], "\n")
  if("SigTerms" %in% names(x))
    cat("    Nontrivial nodes:", x["SigTerms"], "\n")
}

################################################################################
################################################################################

## function to aggregate 2 or more topGOdata objects
combineResults <- function(..., method = c("gmean", "mean", "median", "min", "max")) {

  resList <- list(...)
  ## first for the class of the elements in the list
  if((length(resList) < 2) || !all(sapply(resList, is, "topGOresult")))
    stop("Use: topGOresult_1, topGOresult_2, ..., method = \"mean\"")
  
  combMethod <- match.arg(method)

  retVal <- resList[[1]] 

  ## update the infos
  description(retVal) <- paste(description(retVal),
                               paste(combMethod, "(",
                                     paste(sapply(resList, algorithm), collapse = ", "),
                                     ")", sep = ""),
                               sep = "\n")
  testName(retVal) <- paste(unique(sapply(resList, testName)), collapse = ", ")
  algorithm(retVal) <- "ensemble"
    
  ## get the list of GOs
  nn <- names(score(retVal))
  ## obtain the score from the objects
  resList <- do.call(cbind, lapply(resList, score, whichGO = nn))

  newRes <- switch(combMethod, 
                   mean = rowMeans(resList),
                   median = rowMedians(resList),
                   min = rowMin(resList),
                   max = rowMax(resList),
                   gmean = exp(rowMeans(log(resList))))
  names(newRes) <- rownames(resList) # just to make sure ... some of the functions we use might drop the names
  score(retVal) <- newRes

  return(retVal)
}


.get.sig.gene.in.term<-function(GOdata,my.gene.list){
  genes<-genesInTerm(GOdata)
  lapply(genes,intersect,my.gene.list)
}
