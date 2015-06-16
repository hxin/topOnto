#' ONTdata
#'
#' The \code{ONTdata} data contains just a small number of genes with the corespondent enrichment analysis result.
#' The information on where to find the GO annotations is stored in the \code{geneID2TERM} object.
#' @aliases ONTdata GOdata geneID2TERM geneList geneNames myInterestingGenes resultElimFis resultFis Term allRes
#' @title A toy example of a list of gene identifiers and their enrichment result with human disease ontology
#' @description  The \code{ONTdata} data contains just a small number of genes with the corespondent enrichment analysis result. The information on where to find the GO annotations is stored in the \code{geneID2TERM} object.
#' \itemize{
#'   \item GOdata. 
#'   \item geneID2TERM This object holds the gene annotation, i.e. link gene to ontology terms.
#'   \item myInterestingGenes A set of genes of interests
#'   \item geneList a factor indicate which gene can be found in the annotation data from your list
#'   \item geneNames gene entrez_ids
#'   \item resultElimFis the result from Fisher test when applied elim method
#'   \item resultFis the result from Fisher test
#'   \item Term obtology terms
#'   \item allRes everything together
#' }
#'
#' @name ONTdata
#' @keywords datasets
NULL