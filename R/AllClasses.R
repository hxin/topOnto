#################### topONTdata class ####################
## 
## The node atributes are environments containing the genes/probes
## annotated to the respective node
##
## If genes is a numeric vector than this should represent the gene's
## score. If it is factor it should discriminate the genes in
## interesting genes and the rest

## TODO: it will be a good idea to replace the allGenes and allScore with an
## ExpressionSet class. In this way we can use tests like global test, globalAncova....
##  -- ALL variables sarting with . are just for internal class usage (private)


#' @name topONTdata-class
#' 
#' @docType class
#' 
#' @aliases topONTdata-class allGenes attrInTerm countGenesInTerm description<- description feasible<- feasible geneSelectionFun<- geneSelectionFun genes genesInTerm graph<- graph numGenes ontology<- ontology sigGenes numSigGenes termStat updateGenes updateTerm<- usedGO expressionMatrix phenotype expressionMatrix,topONTdata-method phenotype,topONTdata-method geneScore geneScore,topONTdata,missing-method geneScore,topONTdata,character-method scoresInTerm scoresInTerm,topONTdata,missing-method scoresInTerm,topONTdata,character-method show,topONTdata-method allGenes,topONTdata-method attrInTerm,topONTdata,character,character-method attrInTerm,topONTdata,character,missing-method countGenesInTerm,topONTdata,character-method countGenesInTerm,topONTdata,missing-method description<-,topONTdata,ANY-method description,topONTdata-method feasible<-,topONTdata-method feasible,topONTdata-method geneScore,topONTdata-method geneSelectionFun<-,topONTdata-method geneSelectionFun,topONTdata-method genes,topONTdata-method genesInTerm,topONTdata,character-method genesInTerm,topONTdata,missing-method graph<-,topONTdata-method graph,topONTdata-method initialize,topONTdata-method numGenes,topONTdata-method ontology<-,topONTdata-method ontology,topONTdata-method print,topONTdata-method sigGenes,topONTdata-method numSigGenes,topONTdata-method termStat,topONTdata,character-method termStat,topONTdata,missing-method updateGenes,topONTdata,numeric,function-method updateGenes,topONTdata,factor,missing-method updateTerm<-,topONTdata,character-method usedGO,topONTdata-method
#' 
#' @title Class "topONTdata"
#' @description 
#'   TODO: The node attributes are environments containing the genes/probes
#'   annotated to the respective node
#' 
#'   If genes is a numeric vector than this should represent the gene's
#'   score. If it is factor it should discriminate the genes in
#'   interesting genes and the rest
#'   
#'   TODO: it will be a good idea to replace the allGenes and allScore with an
#'   ExpressionSet class. In this way we can use tests like global test, globalAncova....
#'   -- ALL variables starting with . are just for internal class usage (private)
#'   
#' 
#' 
#' @section Objects from the Class:
#'   Objects can be created by calls of the form \code{new("topONTdata", ontology, allGenes, geneSelectionFun, description, annotationFun, ...)}.
#'   ~~ describe objects here ~~ 
#'
#' 
#' @section Slots:
#'   \describe{
#'     \item{\code{description}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{ontology}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{allGenes}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{allScores}:}{Object of class \code{"ANY"} ~~ }
#'     \item{\code{geneSelectionFun}:}{Object of class \code{"function"} ~~ }
#'     \item{\code{feasible}:}{Object of class \code{"logical"} ~~ }
#'     \item{\code{nodeSize}:}{Object of class \code{"integer"} ~~ }
#'     \item{\code{graph}:}{Object of class \code{"graphNEL"} ~~ }
#'     \item{\code{expressionMatrix}:}{Object of class \code{"matrix"} ~~ }
#'     \item{\code{phenotype}:}{Object of class \code{"factor"} ~~ }
#'     
#'   }
#' 
#' @section Methods:
#'   \describe{
#'     \item{allGenes}{\code{signature(object = "topONTdata")}: ... }
#'     \item{attrInTerm}{\code{signature(object = "topONTdata", attr = "character", whichGO = "character")}: ... }
#'     \item{attrInTerm}{\code{signature(object = "topONTdata", attr = "character", whichGO = "missing")}: ... }
#'     \item{countGenesInTerm}{\code{signature(object = "topONTdata", whichGO = "character")}: ... }
#'     \item{countGenesInTerm}{\code{signature(object = "topONTdata", whichGO = "missing")}: ... }
#'     \item{description<-}{\code{signature(object = "topONTdata")}: ... }
#'     \item{description}{\code{signature(object = "topONTdata")}: ... }
#'     \item{feasible<-}{\code{signature(object = "topONTdata")}: ... }
#'     \item{feasible}{\code{signature(object = "topONTdata")}: ... }
#'     \item{geneScore}{\code{signature(object = "topONTdata")}: ... }
#'     \item{geneSelectionFun<-}{\code{signature(object = "topONTdata")}: ... }
#'     \item{geneSelectionFun}{\code{signature(object = "topONTdata")}: ... }
#' 
#'     \item{genes}{\code{signature(object = "topONTdata")}: A method for
#'       obtaining the list of genes, as a characther vector, which will be
#'       used in the further analysis.}
#'     \item{numGenes}{\code{signature(object = "topONTdata")}: A method for
#'       obtaining the number of genes, which will be used in the further
#'       analysis. It has the same effect as: \code{lenght(genes(object))}.}
#' 
#'     \item{sigGenes}{\code{signature(object = "topONTdata")}:  A method for
#'       obtaining the list of significant genes, as a charachter vector.}
#'         
#'     \item{genesInTerm}{\code{signature(object = "topONTdata", whichGO = "character")}: ... }
#'     \item{genesInTerm}{\code{signature(object = "topONTdata", whichGO = "missing")}: ... }
#'     \item{getSigGroups}{\code{signature(object = "topONTdata", test.stat = "classicCount")}: ... }
#'     \item{getSigGroups}{\code{signature(object = "topONTdata", test.stat = "classicScore")}: ... }
#'     \item{graph<-}{\code{signature(object = "topONTdata")}: ... }
#'     \item{graph}{\code{signature(object = "topONTdata")}: ... }
#'     \item{initialize}{\code{signature(.Object = "topONTdata")}: ... }
#' 
#'     \item{ontology<-}{\code{signature(object = "topONTdata")}: ... }
#'     \item{ontology}{\code{signature(object = "topONTdata")}: ... }
#' 
#'     \item{termStat}{\code{signature(object = "topONTdata", whichGO = "character")}: ... }
#'     \item{termStat}{\code{signature(object = "topONTdata", whichGO = "missing")}: ... }
#'     \item{updateGenes}{\code{signature(object = "topONTdata", geneList = "numeric", geneSelFun = "function")}: ... }
#'     \item{updateGenes}{\code{signature(object = "topONTdata", geneList = "factor", geneSelFun = "missing")}: ... }
#'     \item{updateTerm<-}{\code{signature(object = "topONTdata", attr = "character")}: ... }
#'     \item{usedGO}{\code{signature(object = "topONTdata")}: ... }
#'   }
#' 
#' @author Adrian Alexa
#' 
#' @seealso
#'   \code{\link{buildLevels}},
#'   \code{\link{annFUN}}
#' 
#' @examples
#' ## load the dataset 
#' data(ONTdata)
#' #library(package = affyLib, character.only = TRUE)
#' 
#' ## the distribution of the adjusted p-values
#' #hist(geneList, 100)
#' 
#' ## how many differentially expressed genes are:
#' #sum(topDiffGenes(geneList))
#' 
#' ## build the topONTdata class 
#' \dontrun{
#' require(topOnto.HDO.db)
#' GOdata <- new("topONTdata", ontology = "HDO", description="HDO example",allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2TERM)
#' ## display the GOdata object
#' GOdata
#' }
#' 
#' ##########################################################
#' ## Examples on how to use the methods
#' ##########################################################
#' 
#' ## description of the experiment
#' description(GOdata)
#' 
#' ## obtain the genes that will be used in the analysis
#' a <- genes(GOdata)
#' str(a)
#' numGenes(GOdata)
#' 
#' ## obtain the score (p-value) of the genes
#' selGenes <- names(geneList)[sample(1:length(geneList), 10)]
#' gs <- geneScore(GOdata, whichGenes = selGenes)
#' print(gs)
#' 
#' ## if we want an unnamed vector containing all the feasible genes
#' gs <- geneScore(GOdata, use.names = FALSE)
#' str(gs)
#' 
#' ## the list of significant genes
#' sg <- sigGenes(GOdata)
#' str(sg)
#' numSigGenes(GOdata)
#' 
#' ## to update the gene list 
#' .geneList <- geneScore(GOdata, use.names = TRUE)
#' GOdata ## more available genes
#' #GOdata <- updateGenes(GOdata, .geneList, topDiffGenes)
#' GOdata ## the available genes are now the feasible genes
#' 
#' ## the available GO terms (all the nodes in the graph)
#' go <- usedGO(GOdata)
#' length(go)
#' 
#' ## to list the genes annotated to a set of specified GO terms
#' sel.terms <- sample(go, 10)
#' ann.genes <- genesInTerm(GOdata, sel.terms)
#' str(ann.genes)
#' 
#' ## the score for these genes
#' ann.score <- scoresInTerm(GOdata, sel.terms)
#' str(ann.score)
#' 
#' ## to see the number of annotated genes
#' num.ann.genes <- countGenesInTerm(GOdata)
#' str(num.ann.genes)
#' 
#' ## to summarise the statistics
#' termStat(GOdata, sel.terms)
#' 
#' @keywords graphs,classes
setClass("topONTdata",
        representation = representation(
        ## some description of the data/experiment
        description = "character",
        ## which of the ontology to be used: BP, CC, MF, (BC, BM, CM, ...) 
        ontology = "character",
        ## the vector containing the genes/probes used in the experiment
        allGenes = "character",
        ## the vector containing the score or the  of the  genes/probes 
        allScores = "ANY",
        ## if each has a score, geneSelectionFun: function(x, ...)
        geneSelectionFun = "function",
        ## which genes are annotated to the specific ontology
        feasible = "logical",
        ## the GO ontology: graph structure and annotations: leaves to root
        graph = "graphNEL",
        ## nodes in the graph with fewer than "nodeSize" genes are removed
        nodeSize = "integer",
        ## expression matrix          #!
        expressionMatrix = "ANY",        #!
        ## phenotype information (groups that shall be compared) #!
        phenotype = "ANY",
        ## ontology id2term mapping
        termName="ANY",
        ##when using gsea, the gct file(or dataframe by read.gct)
        gct ="ANY",
        exp = "ANY",
        ##when using gsea, the cls file(or dataframe by read.cls)
        cls ="ANY",
        pty ="ANY",
        useScore = "ANY"
        ))        #!



######################## topONTresult class ######################
## probably will add more infos with time here
#' @name topONTresult-class
#' @docType class
#' @aliases topONTresult-class initialize,topONTresult-method description,topONTresult-method description<-,topONTresult,ANY-method score score,topONTresult-method score<-,topONTresult-method testName testName<- testName,topONTresult-method testName<-,topONTresult-method algorithm algorithm<- algorithm<-,topONTresult-method algorithm,topONTresult-method show,topONTresult-method print,topONTresult-method geneData<- geneData geneData<-,topONTresult-method geneData,topONTresult-method combineResults
#' @title Class "topONTresult"
#' 
#' @description Class instance created by \code{\link{getSigGroups-methods}} or by \code{runTest}
#' 
#' @section Objects from the Class:
#'   Objects can be created by calls of the form \code{new("topONTresult",
#'     description, score, testName, algorithm, geneData)}.
#' 
#' 
#' @section Slots:
#'   \describe{
#'     \item{\code{description}:}{character string containing a
#'       short description on how the object was build.}
#'     \item{\code{score}:}{named numerical vector containing the
#'       p-values or the scores of the tested GO terms.}
#'     \item{\code{testName}:}{character string containing the name
#'       of the test statistic used.}
#'     \item{\code{algorithm}:}{character string containing the name
#'       of the algorithm used.}
#'     \item{\code{geneData}:}{list containing summary statistics
#'       on the genes/gene universe/annotations.}
#'   }
#' 
#' 
#' @section Methods:
#'   \describe{
#'     \item{\code{score}:}{method to access the \code{score} slot.}
#'     
#'     \item{\code{testName}:}{method to access the \code{testName} slot.}
#' 
#'     \item{\code{algorithm}:}{method to access the \code{algorithm} slot.}
#' 
#'     \item{\code{geneData}:}{method to access the \code{geneData} slot.}
#' 
#'     \item{\code{show}:}{method to print the object.}
#' 
#'     \item{\code{combineResults}:}{method to aggregate two or more
#'       topONTresult objects. \code{method = c("gmean", "mean", "median",
#'       "min", "max")} provides the way the object scores (which most of
#'       the time are p-values) are combined.}.
#' }
#' 
#' 
#' @author Adrian Alexa
#' 
#' @seealso
#'   \code{\link{groupStats-class}},
#'   \code{\link{getSigGroups-methods}}
#' 
#' 
#' @examples
#' 
#' data(ONTdata)
#' 
#' s <- score(resultFis)
#' 
#' go <- sort(names(s))
#' go.sub<- sample(go, 100)
#' go.mixed <- c(sample(go, 50), sample(ls(Term), 20))
#' go.others <- sample(ls(Term), 100)
#' 
#' 
#' str(go)
#' str(go.sub)
#' str(go.mixed)
#' str(go.others)
#' 
#' str(score(resultFis, whichGO = go))
#' str(score(resultFis, whichGO = go.sub))
#' str(score(resultFis, whichGO = go.mixed))
#' str(score(resultFis, whichGO = go.others))
#' 
#' avgResult <- combineResults(resultFis, resultElimFis)
#' avgResult
#' combineResults(resultFis, resultElimFis, method = "min")
#' 
#' 
#' @keywords classes
setClass("topONTresult",
         representation = representation(
        ## some description of the data/experiment
        description = "character",
        ## the p-values of the GO terms (named vector)
        score = "numeric",
        ## which test statistic was used 
        testName = "character",
        ## which algorithm was used 
        algorithm = "character",
        ## stats about the genes ...
        geneData = "ANY"))



######################## groupInfo class ########################
## A virtual class containing basic group (GO term) data:
##     gene names, genes scores, etc...
#' @name groupStats-class
#' @docType class
#' @aliases groupStats-class allMembers<- allMembers members<- members Name<- Name numAllMembers numMembers testStatistic testStatPar updateGroup allMembers<-,groupStats-method allMembers,groupStats-method initialize,groupStats-method members<-,groupStats-method members,groupStats,missing-method Name<-,groupStats-method Name,groupStats-method numAllMembers,groupStats-method numMembers,groupStats-method runTest,groupStats-method runTest,groupStats,missing,missing-method testStatistic,groupStats-method testStatPar,groupStats-method updateGroup,groupStats,character,character-method
#' @title  Class "groupStats"
#' @description A virtual class containing basic gene set information:
#'   the gene universe, the member of the current group, the test statistic
#'   defined for this group, etc. 
#' 
#' 
#' @section Objects from the Class:A virtual Class: No objects may be created from it.
#' @section Slots:
#'   \describe{
#'     \item{\code{name}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{allMembers}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{members}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{testStatistic}:}{Object of class \code{"function"} ~~ }
#'     \item{\code{testStatPar}:}{Object of class \code{"ANY"} ~~ }      
#'   }
#' 
#' 
#' @section Methods:
#'   \describe{
#'     \item{allMembers<-}{\code{signature(object = "groupStats")}: ... }
#'     \item{allMembers}{\code{signature(object = "groupStats")}: ... }
#'     \item{initialize}{\code{signature(.Object = "groupStats")}: ... }
#'     \item{members<-}{\code{signature(object = "groupStats")}: ... }
#'     \item{members}{\code{signature(object = "groupStats")}: ... }
#'     \item{Name<-}{\code{signature(object = "groupStats")}: ... }
#'     \item{Name}{\code{signature(object = "groupStats")}: ... }
#'     \item{numAllMembers}{\code{signature(object = "groupStats")}: ... }
#'     \item{numMembers}{\code{signature(object = "groupStats")}: ... }
#'     \item{runTest}{\code{signature(object = "groupStats")}: ... }
#'     \item{testStatistic}{\code{signature(object = "groupStats")}: ... }
#'    }
#' 
#' 
#' @author Adrian Alexa
#' 
#' @examples
#' test.stat=new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
#' 
#' @seealso
#'   \code{\link{classicCount-class}},
#'   \code{\link{getSigGroups-methods}}
#' 
#' 
#' @keywords classes
setClass("groupStats",
        representation = representation(
        ## the name of the group (GO ID for example)
        name = "character",
        ## the names of all the mebers (gene ID)
        allMembers = "character",
        ## the names of the mebers in the group (gene ID)
        members = "character",
        ## function containg a statistical test takeing the
        ## parameter the object itself 
        testStatistic = "function",
        testStatPar = "list",
        "VIRTUAL"))

#################### classicCount class #########################
## A class that extends the virtual class groupStats by adding 
## a slot representing the significant members.
#'@name classicCount-class
#'@docType class
#'@aliases classicCount-class contTable numSigAll numSigMembers sigAllMembers sigMembers<- sigMembers contTable,classicCount-method initialize,classicCount-method numSigAll,classicCount-method numSigMembers,classicCount-method sigAllMembers,classicCount-method sigMembers<-,classicCount-method sigMembers,classicCount-method GOFisherTest,classicCount-method
#'@title Class "classicCount"
#'@description This class that extends the virtual class "groupStats" 
#'  by adding a slot representing the significant members.
#'@section Objects from the Class:
#'  Objects can be created by calls of the form
#'  \code{new("classicCount",
#'          testStatistic = "function",
#'          name = "character",
#'          allMembers = "character",
#'          groupMembers = "character",
#'          sigMembers = "character")}.
#'          
#' @section Slots:
#' \describe{
#'  \item{significant}{Object of class \code{"integer"}}
#'  \item{name}{Object of class \code{"character"}}
#'  \item{allMembers}{Object of class \code{"character"}}
#'  \item{members}{Object of class \code{"character"}}
#'  \item{testStatistic}{Object of class \code{"function"}}
#'   }
#'
#' @section Extends:
#'  Class \code{"groupStats"}, directly.
#' @section Methods:
#'  \describe{
#'    \item{contTable}{\code{signature(object = "classicCount")}: ... }
#'    \item{initialize}{\code{signature(.Object = "classicCount")}: ... }
#'    \item{numSigAll}{\code{signature(object = "classicCount")}: ... }
#'    \item{numSigMembers}{\code{signature(object = "classicCount")}: ... }
#'    \item{sigAllMembers}{\code{signature(object = "classicCount")}: ... }
#'    \item{sigMembers<-}{\code{signature(object = "classicCount")}: ... }
#'    \item{sigMembers}{\code{signature(object = "classicCount")}: ... }
#'  }
#' @details
#'  This class is used for test statistic based on counts, like Fisher's
#'  exact test
#' @author Adrian Alexa
#' @seealso
#'    \code{\link{classicScore-class}},
#'    \code{\link{groupStats-class}},
#'    \code{\link{getSigGroups-methods}}
#'    
#' @examples
#' data(ONTdata)
#' gene.universe <- genes(GOdata)
#' go.genes <- genesInTerm(GOdata, 'DOID:0050709')[[1]]
#' sig.genes <- sigGenes(GOdata)
#' my.group <- new("classicCount", testStatistic = GOFisherTest, name = "fisher",allMembers = gene.universe, groupMembers = go.genes,sigMembers = sig.genes)
#' contTable(my.group)
#' runTest(my.group)
#' 
#' @keywords classes
setClass("classicCount", contains = "groupStats", 
        representation = representation(
        ## the index of which members are significant  
        significant = "integer"))


##################### classicScore class ########################
## A class that extends the virtual class groupStats by adding 
## a slot representing the score of each gene. (used for KS test)
#' @name classicScore-class
#' @docType class
#' @aliases classicScore-class scoreOrder membersScore rankMembers score<- allScore allScore,classicScore,missing-method allScore,classicScore,logical-method scoreOrder,classicScore-method initialize,classicScore-method membersScore,classicScore-method rankMembers,classicScore-method score<-,classicScore-method GOKSTest,classicScore-method GOtTest,classicScore-method GOSumTest,classicScore-method GOKSTiesTest,classicScore-method
#' @title Class "classicScore"
#' @description A class that extends the virtual class "groupStats" by adding 
#'   a slot representing the score of each gene. It is used for tests like
#'   Kolmogorov-Smirnov test.
#' 
#' @section Objects from the Class:
#'   Objects can be created by calls of the form \code{new("classicScore", testStatistic, name, allMembers, groupMembers, score, decreasing)}.
#' 
#' @section Slots:
#'    \describe{
#'     \item{\code{score}:}{Object of class \code{"numeric"} ~~ }
#'     \item{\code{name}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{allMembers}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{members}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{testStatistic}:}{Object of class \code{"function"} ~~ }
#'     \item{\code{scoreOrder}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{testStatPar}:}{Object of class \code{"ANY"} ~~ }
#'   }
#'   
#' @section Extends:
#' Class \code{"groupStats"}, directly.
#' 
#' @section Methods:
#'   \describe{
#'     \item{allScore}{Method to obtain the score of all members.}
#'     \item{scoreOrder}{Returns TRUE if the score should be ordered
#'       increasing, FALSE otherwise.}
#'     \item{membersScore}{\code{signature(object = "classicScore")}: ... }
#'     \item{rankMembers}{\code{signature(object = "classicScore")}: ... }
#'     \item{score<-}{\code{signature(object = "classicScore")}: ... }
#'  }
#' 
#' @author Adrian Alexa
#' 
#' @seealso
#'   \code{\link{classicCount-class}},
#'   \code{\link{groupStats-class}},
#'   \code{\link{getSigGroups-methods}}
#' 
#' @examples
#' ## define the type of test you want to use
#' test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
#' 
#' @keywords classes
setClass("classicScore", contains = "groupStats",
        representation = representation(
        ## the score for each member, the most important
        ## member has the highest score
        score = "numeric",
        ## scoreOrder = TRUE
        ##(decreasing) the max(score) is considering the best score
        ## scoreOrder = FALSE
        ##(increasing) the min(score) is considre the best score
        scoreOrder = "logical"))


##################### classicExpr class #########################
## A class that extends the virtual class groupStats. Contains:
##  - a slot representing the gene expression matrix (for all genes).
##  - a slot representing the phenotipic data (can be a matrix or a vector!).
## Used for GlobalTest.

## the expression matrix will not have row names. The row names will be stored in
## allMembers slot of the parent class
#' @name classicExpr-class
#' @docType class
#' @aliases classicExpr-class pType<- pType membersExpr allMembers<-,classicExpr-method emptyExpr,classicExpr-method GOglobalTest,classicExpr-method initialize,classicExpr-method membersExpr,classicExpr-method pType<-,classicExpr-method pType,classicExpr-method
#' @title Class "classicExpr"
#' @description This class that extends the virtual class "groupStats" by adding two slots for accomodating gene expression data.
#' 
#' @section Objects from the Class:
#'   Objects can be created by calls of the form \code{new("classicExpr", testStatistic, name, groupMembers, exprDat, pType, ...)}.
#' 
#' @section Slots:
#'    \describe{
#'     \item{\code{eData}:}{Object of class \code{"environment"} ~~ }
#'     \item{\code{pType}:}{Object of class \code{"factor"} ~~ }
#'     \item{\code{name}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{allMembers}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{members}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{testStatistic}:}{Object of class \code{"function"} ~~ }
#'     \item{\code{testStatPar}:}{Object of class \code{"list"} ~~ }
#'   }
#'   
#' @section Extends:
#' Class \code{"\linkS4class{groupStats}"}, directly.
#' @section Methods:
#'   \describe{
#'     \item{allMembers<-}{\code{signature(object = "classicExpr")}: ... }
#'     \item{emptyExpr}{\code{signature(object = "classicExpr")}: ... }
#'     \item{getSigGroups}{\code{signature(object = "topONTdata", test.stat = "classicExpr")}: ... }
#'     \item{GOglobalTest}{\code{signature(object = "classicExpr")}: ... }
#'     \item{initialize}{\code{signature(.Object = "classicExpr")}: ... }
#'     \item{membersExpr}{\code{signature(object = "classicExpr")}: ... }
#'     \item{pType<-}{\code{signature(object = "classicExpr")}: ... }
#'     \item{pType}{\code{signature(object = "classicExpr")}: ... }
#'   }
#' 
#' @author Adrian Alexa
#' 
#' @seealso
#'   \code{\link{classicScore-class}},
#'   \code{\link{groupStats-class}},
#'   \code{\link{getSigGroups-methods}}
#' 
#' @examples
#' showClass("classicExpr")
#' @keywords classes
setClass("classicExpr", contains = "groupStats",
        representation = representation(
        ## the expression matrix
        eData = "environment", ## we store the expr matrix in this environment for fast class copying 
        ## scoreOrder = TRUE (decreasing) the max(score) is considering the best score
        pType = "factor"))




##################### weight01Count class #########################
## used for elim2 algorithm
setClass("weight01Count", contains = "classicCount", 
        representation = representation(
        ## the index of which of the group members should be removed
        elim = "integer"))


##################### weight01Score class #########################
setClass("weight01Score", contains = "classicScore",
        representation = representation(
        ## the score for each member, the most important
        ## member has the highest score
        elim = "integer"))


##################### weight01Expr class #########################
setClass("weight01Expr", contains = "classicExpr",
        representation = representation(
        elim = "integer"))



####################### elimCount class #########################
#' @name elimCount-class
#' @docType class
#' @aliases leaCount-class elimCount-class weight01Count-class cutOff cutOff,elimCount-method cutOff<- cutOff<-,elimCount-method elim elim,elimCount-method elim,weight01Count-method elim<- elim<-,elimCount-method elim<-,weight01Count-method depth depth,leaCount-method depth<- depth<-,leaCount-method initialize,leaCount-method initialize,elimCount-method initialize,weight01Count-method contTable,elimCount-method numSigAll,elimCount-method numSigMembers,elimCount-method numMembers,elimCount-method numAllMembers,elimCount-method sigAllMembers,elimCount-method sigMembers<-,elimCount-method sigMembers,elimCount-method numAllMembers,weight01Count-method numMembers,weight01Count-method numSigAll,weight01Count-method numSigMembers,weight01Count-method sigAllMembers,weight01Count-method sigMembers,weight01Count-method GOFisherTest,elimCount-method
#' @title Classes "elimCount" and "weight01Count"
#' @description Classes that extend the "classicCount" class by adding a slot representing the members that need to be removed.
#' 
#' @section Objects from the Class:
#'   Objects can be created by calls of the form \code{new("elimCount", testStatistic, name, allMembers, groupMembers, sigMembers, elim, cutOff, ...)}.
#' 
#' @section Slots:
#'   \describe{
#'     \item{\code{elim}:}{Object of class \code{"integer"} ~~ }
#'     \item{\code{cutOff}:}{Object of class \code{"numeric"} ~~ }
#'     \item{\code{significant}:}{Object of class \code{"integer"} ~~ }
#'     \item{\code{name}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{allMembers}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{members}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{testStatistic}:}{Object of class \code{"function"} ~~ }
#'     \item{\code{testStatPar}:}{Object of class \code{"list"} ~~ }
#'   }

#' @section Extends:
#'   Class \code{"\linkS4class{classicCount}"}, directly.
#' Class \code{"\linkS4class{groupStats}"}, by class "classicCount", distance 2.
#' 
#' @section Methods:
#'   No methods defined with class "elimCount" in the signature.
#' 
#' @author Adrian Alexa
#' 
#' @examples
#' data(ONTdata)
#' gene.universe <- genes(GOdata)
#' go.genes <- genesInTerm(GOdata, 'DOID:0050709')[[1]]
#' elim.genes <- sample(go.genes, length(go.genes) / 4)
#' sig.genes <- sigGenes(GOdata)
#' elim.group <- new("elimCount", testStatistic = GOFisherTest, name = "fisher",allMembers = gene.universe, groupMembers = go.genes,sigMembers = sig.genes, elim = elim.genes)
#' contTable(elim.group)
#' runTest(elim.group)
#' 
#' @seealso
#'   \code{\link{classicScore-class}},
#'   \code{\link{groupStats-class}},
#'   \code{\link{getSigGroups-methods}}
#' 
#' @keywords classes
setClass("elimCount", contains = "weight01Count", 
        representation = representation(
        cutOff = "numeric"))


####################### elimScore class #########################
#' @name elimScore-class
#' @docType class
#' @aliases leaScore-class elimScore-class weight01Score-class allMembers,elimScore-method allScore,elimScore,logical-method allScore,elimScore,missing-method cutOff,elimScore-method cutOff<-,elimScore-method elim,elimScore-method elim<-,elimScore-method depth,leaScore-method depth<-,leaScore-method alternative,elimScore-method membersScore,elimScore-method members,elimScore-method numMembers,elimScore-method numAllMembers,elimScore-method rankMembers,elimScore-method score<-,elimScore-method initialize,leaScore-method initialize,elimScore-method initialize,weight01Score-method allMembers,weight01Score-method allScore,weight01Score,missing-method allScore,weight01Score,logical-method elim<-,weight01Score-method elim,weight01Score-method members,weight01Score,missing-method membersScore,weight01Score-method numAllMembers,weight01Score-method numMembers,weight01Score-method rankMembers,weight01Score-method
#' @title Classes "elimScore" and "weight01Score" 
#' @description Classes that extend the "classicScore" class by adding a slot representing the members that need to be removed.
#' 
#' @section Objects from the Class:
#' Objects can be created by calls of the form \code{new("elimScore", testStatistic, name, allMembers, groupMembers, score, alternative, elim, cutOff, ...)}.
#'    ~~ describe objects here ~~ 
#' 
#' @section Slots:
#'  \describe{
#'     \item{\code{elim}:}{Object of class \code{"integer"} ~~ }
#'     \item{\code{cutOff}:}{Object of class \code{"numeric"} ~~ }
#'     \item{\code{score}:}{Object of class \code{"numeric"} ~~ }
#'     \item{\code{.alternative}:}{Object of class \code{"logical"} ~~ }
#'     \item{\code{name}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{allMembers}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{members}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{testStatistic}:}{Object of class \code{"function"} ~~ }
#'     \item{\code{testStatPar}:}{Object of class \code{"list"} ~~ }
#'   }
#'   
#' @section Extends:
#' Class \code{"\linkS4class{classicScore}"}, directly.
#' Class \code{"\linkS4class{groupStats}"}, by class "classicScore", distance 2.
#' 
#' @section Methods:
#' No methods defined with class "elimScore" in the signature.
#' 
#' 
#' @author Adrian Alexa
#' 
#' @seealso
#'   \code{\link{classicScore-class}},
#'   \code{\link{groupStats-class}},
#'   \code{\link{getSigGroups-methods}}
#' 
#' 
#' @examples
#' ##---- Should be DIRECTLY executable !! ----
#' 
#' @keywords classes
setClass("elimScore", contains = "weight01Score",
        representation = representation(
        cutOff = "numeric")) 


####################### elimExpr class ##########################
#' @name elimExpr-class
#' @docType class
#' @aliases leaExpr-class elimExpr-class weight01Expr-class depth,leaExpr-method depth<-,leaExpr-method cutOff<-,elimExpr-method cutOff,elimExpr-method initialize,leaExpr-method initialize,elimExpr-method initialize,weight01Expr-method allMembers,weight01Expr-method elim<-,weight01Expr-method elim,weight01Expr-method members,weight01Expr,missing-method numAllMembers,weight01Expr-method numMembers,weight01Expr-method
#' @title Class "elimExpr"
#' @description Classes that extend the "classicExpr" class by adding a slot representing the members that need to be removed.
#' 
#' @section Objects from the Class:
#' Objects can be created by calls of the form \code{new("elimExpr", testStatistic, name, groupMembers, exprDat, pType, elim, cutOff, ...)}.
#'    ~~ describe objects here ~~ 
#' 
#' @section Slots:
#'  \describe{
#'     \item{\code{cutOff}:}{Object of class \code{"numeric"} ~~ }
#'     \item{\code{elim}:}{Object of class \code{"integer"} ~~ }
#'     \item{\code{eData}:}{Object of class \code{"environment"} ~~ }
#'     \item{\code{pType}:}{Object of class \code{"factor"} ~~ }
#'     \item{\code{name}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{allMembers}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{members}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{testStatistic}:}{Object of class \code{"function"} ~~ }
#'     \item{\code{testStatPar}:}{Object of class \code{"list"} ~~ }
#'   }
#' 
#' @section Extends:
#' Class \code{"\linkS4class{weight01Expr}"}, directly.
#' Class \code{"\linkS4class{classicExpr}"}, by class "weight01Expr", distance 2.
#' Class \code{"\linkS4class{groupStats}"}, by class "weight01Expr", distance 3.
#' 
#' @section Methods:
#'   \describe{
#'     \item{cutOff<-}{\code{signature(object = "elimExpr")}: ... }
#'     \item{cutOff}{\code{signature(object = "elimExpr")}: ... }
#'     \item{getSigGroups}{\code{signature(object = "topONTdata", test.stat = "elimExpr")}: ... }
#'     \item{initialize}{\code{signature(.Object = "elimExpr")}: ... }
#'  }
#' 
#' 
#' @author Adrian Alexa
#' 
#' @seealso
#'   \code{\link{classicScore-class}},
#'   \code{\link{groupStats-class}},
#'   \code{\link{getSigGroups-methods}}
#' 
#' 
#' @examples
#' showClass("elimExpr")
#' 
#' @keywords classes
setClass("elimExpr", contains = "weight01Expr",
        representation = representation(
        cutOff = "numeric"))



###################### weightCount class ########################
#' @name weightCount-class
#' @docType class
#' 
#' @aliases  weightCount-class getSigRatio sigRatio penalise roundFun significant sigRatio<- Weights Weights<- Weights,weightCount-method Weights,weightCount,missing-method Weights,weightCount,logical-method Weights<-,weightCount-method sigRatio,weightCount-method sigRatio<-,weightCount-method penalise,weightCount,numeric,numeric-method roundFun,weightCount-method significant,weightCount-method Name,weightCount-method allMembers,weightCount-method members,weightCount-method testStatistic,weightCount-method testStatPar,weightCount-method updateGroup,weightCount,character,character-method getSigRatio,weightCount-method initialize,weightCount-method numSigAll,weightCount-method numMembers,weightCount-method numSigMembers,weightCount-method numAllMembers,weightCount-method
#' 
#' @title Class "weightCount"
#' 
#' @description    ~~ A concise (1-5 lines) description of what the class is.
#' 
#' @section Objects from the Class:
#'   Objects can be created by calls of the form \code{new("weightCount", testStatistic, name, allMembers, groupMembers, sigMembers, weights, sigRatio, penalise, ...)}.
#'    
#' @section Slots:
#'  \describe{
#'     \item{\code{weights}:}{Object of class \code{"numeric"} ~~ }
#'     \item{\code{sigRatio}:}{Object of class \code{"function"} ~~ }
#'     \item{\code{penalise}:}{Object of class \code{"function"} ~~ }
#'     \item{\code{roundFun}:}{Object of class \code{"function"} ~~ }
#'     \item{\code{significant}:}{Object of class \code{"integer"} ~~ }
#'     \item{\code{name}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{allMembers}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{members}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{testStatistic}:}{Object of class \code{"function"} ~~ }
#'     \item{\code{testStatPar}:}{Object of class \code{"list"} ~~ }
#'   }
#' 
#' @section Extends:
#'   Class \code{"\linkS4class{classicCount}"}, directly.
#'   Class \code{"\linkS4class{groupStats}"}, by class "classicCount", distance 2.
#' 
#' @section Methods:
#'   No methods defined with class "weightCount" in the signature.
#' 
#' 
#' @author Adrian Alexa
#' 
#' @examples
#' data(ONTdata)
#' test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
#' resultWeight <- getSigGroups(GOdata, test.stat)
#' resultWeight
#' 
#' @seealso
#'   \code{\link{groupStats-class}},
#'   \code{\link{getSigGroups-methods}}
#' 
#' 
#' @keywords classes
setClass("weightCount", contains = "classicCount", 
        representation = representation(
        weights = "numeric",
        sigRatio = "function",
        ## if more sig nodes ar penalized more (should be part of the sigRatio) 
        penalise = "function",  ## penalise = c("more", "less")
        roundFun = "function"))


#################### parentChild class #########################
##
## the parents information is stored in the class' slots
##
## allMembers - character vector containing all genes annotated in the parent(s).
##        If the group has two or more parents, the genes are
##        concatenated in this vector. The separation between
##        parents is given by the splitIndex.
## splitIndex - the possition in the allMembers where a new parent starts
## joinFun - the way the genes in the parents are joinded: union or intersection
#' @name parentChild-class
#' @docType class
#' @aliases  parentChild-class pC-class allParents joinFun allMembers<-,parentChild-method allMembers,parentChild-method allParents,parentChild-method initialize,parentChild-method joinFun,parentChild-method numAllMembers,parentChild-method numSigAll,parentChild-method sigAllMembers,parentChild-method sigMembers<-,parentChild-method updateGroup,parentChild,missing,character-method allMembers<-,pC-method initialize,pC-method sigMembers<-,pC-method updateGroup,pC,missing,character-method updateGroup,pC,missing,missing-method
#' @title Classes "parentChild" and "pC"
#' @description Classes that extend the "classicCount" class by adding support for the parent-child test.
#' 
#' @section Objects from the Class:
#'   Objects can be created by calls of the form \code{new("parentChild", testStatistic, name, groupMembers, parents, sigMembers, joinFun, ...)}.
#' 
#' @section Slots:
#'    \describe{
#'     \item{\code{splitIndex}:}{Object of class \code{"integer"} ~~ }
#'     \item{\code{joinFun}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{significant}:}{Object of class \code{"integer"} ~~ }
#'     \item{\code{name}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{allMembers}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{members}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{testStatistic}:}{Object of class \code{"function"} ~~ }
#'     \item{\code{testStatPar}:}{Object of class \code{"list"} ~~ }
#'   }
#'   
#' @section Extends:
#' Class \code{"\linkS4class{classicCount}"}, directly.
#' Class \code{"\linkS4class{groupStats}"}, by class "classicCount", distance 2.
#' 
#' @section Methods:
#'   \describe{
#'     \item{allMembers<-}{\code{signature(object = "parentChild")}: ... }
#'     \item{allMembers}{\code{signature(object = "parentChild")}: ... }
#'     \item{allParents}{\code{signature(object = "parentChild")}: ... }
#'     \item{getSigGroups}{\code{signature(object = "topONTdata", test.stat = "parentChild")}: ... }
#'     \item{initialize}{\code{signature(.Object = "parentChild")}: ... }
#'     \item{joinFun}{\code{signature(object = "parentChild")}: ... }
#'     \item{numAllMembers}{\code{signature(object = "parentChild")}: ... }
#'     \item{numSigAll}{\code{signature(object = "parentChild")}: ... }
#'     \item{sigAllMembers}{\code{signature(object = "parentChild")}: ... }
#'     \item{sigMembers<-}{\code{signature(object = "parentChild")}: ... }
#'     \item{updateGroup}{\code{signature(object = "parentChild", name = "missing", members = "character")}: ... }
#'  }
#' 
#' @author Adrian Alexa
#' 
#' @seealso
#'   \code{\link{classicCount-class}},
#'   \code{\link{groupStats-class}},
#'   \code{\link{getSigGroups-methods}}
#' 
#' @examples
#' showClass("parentChild")
#' showClass("pC")
#' 
#' @keywords classes
setClass("parentChild", contains = "classicCount", 
        representation = representation(splitIndex = "integer",
        joinFun = "character"))


## implements same functionality as parentChild but the join operation
## for the parents is done when the class is build (or updated)
setClass("pC", contains = "classicCount", 
        representation = representation(
        joinFun = "function"))


##################### LEA classes #########################
## used for LEA algorithm
## basically the same class as weight01xxx but 
## we need it for the dispatch mechanism 
setClass("leaCount", contains = "weight01Count",
        representation = representation(
        depth = "integer"))

setClass("leaScore", contains = "weight01Score",
        representation = representation(
        depth = "integer"))

setClass("leaExpr", contains = "weight01Expr",
        representation = representation(
        depth = "integer"))

##################### GSEA classes #########################
setClass("classicGsea", contains = "elimScore",
         representation = representation(
           ## the score for each member, the most important
           ## member has the highest score
           score = "numeric",
           ## scoreOrder = TRUE
           ##(decreasing) the max(score) is considering the best score
           ## scoreOrder = FALSE
           ##(increasing) the min(score) is considre the best score
           scoreOrder = "logical",
           annotation.weight ="numeric",
           min.size ="numeric",
           max.size ="numeric"))

# setClass("weight01Gsea", contains = "classicGsea",
#          representation = representation(
#            ## the score for each member, the most important
#            ## member has the highest score
#            elim = "integer"))

setClass("elimGsea", contains = "classicGsea",
         representation = representation(
           cutOff = "numeric",
           annotation.weight ="numeric",
           elim.type = "character",
           elim.gene.type = "character"))

setClass("topONTresultGSEA",contains = "topONTresult",
         representation = representation(
           global.report = "ANY",
           gs.report = "ANY",
           plots = "ANY",
           cutOff = "numeric"))
