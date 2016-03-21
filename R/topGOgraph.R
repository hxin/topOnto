#' @name Determines the levels of a Directed Acyclic Graph (DAG)
#' @rdname dagFunctions
#' @aliases buildLevels getNoOfLevels getGraphRoot reverseArch
#' @title Utility functions to work with Directed Acyclic Graphs (DAG)
#' @description Basic functions to work witg DAGs
#' 
#' @usage
#' buildLevels(dag, root = NULL, leafs2root = TRUE)
#' getNoOfLevels(graphLevels)
#' getGraphRoot(dag, leafs2root = TRUE)
#' reverseArch(dirGraph, useAlgo = "sparse", useWeights = TRUE)
#' 
#' 
#' @param dag A \code{graphNEL} object.
#'   
#' @param root A character vector specifing the root(s) of the DAG. If
#'     not specified the root node is autmatically computed.
#' 
#' @param leafs2root The leafs2root parameter tell if the graph has edges
#'     directed from the leaves to the root, or vice-versa
#' 
#' @param graphLevels An object of type list, returned by the \code{buildLevels} function.
#'   
#' @param dirGraph A \code{graphNEL} object containing a directed graph.
#'  
#' @param useAlgo A character string specifing one of the following options
#'     \code{c("sparse", "normal")}. By default, \code{useAlgo = "sparse"}, a sparce matrix object is
#'     used to transpose the adjacency matrix. Otherwise a standard R
#'     martix is used.
#'   
#' @param useWeights If weights should be used (if \code{useAlgo =
#'       "normal"} then the weigths are used anyway)
#' 
#' 
#' @details
#'   \code{buildLevels} function determines the levels of a Directed
#'   Acyclic Graph (DAG). The level of a node is defined as the longest
#'   path from the node to the root. The function take constructs a named
#'   list containg varios information about each nodes level. The root has
#'   level 1.
#' 
#'   \code{getNoOfLevels} - a convenient function to extract the number of
#'   levels from the object returned by \code{buildLevels}
#'   
#'   \code{getGraphRoot} finds the root(s) of the DAG
#' 
#'   \code{reverseArch} - simple function to invert the direction of edges
#'   in a DAG. The returned graph is of class graphNEL.  It can use either
#'   simple matrices or sparse matrices (SparseM library)
#' 
#' 
#' @return
#'   \code{buildLevels} returns a list containing:
#'   \item{level2nodes}{Environment where the key is the level number
#'     with the value being the nodes on that level.}
#'   \item{nodes2level}{Environment where the key is the node label (the GO
#'     ID) and the value is the level on which that node lies.}
#'   \item{noOfLevels}{The number of levels}
#'   \item{noOfNodes}{The number of nodes}
#'   An object of class \code{\link[graph]{graphNEL-class}} is returned.
#'       
#' 
#' @author Adrian Alexa
#' 
#' @seealso
#'   \code{\link{topONTdata-class}},
#'   \code{\link{inducedGraph}}
#' 
#' 
#' @examples
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--  or do  help(data=index)  for the standard data sets.
#' NULL
#' @keywords graphs
NULL


############################   reverseArch   ############################
## Simple function to invert the direction of edges in an directed graph.
## The returned graph is of class graphNEL.
## It can use either simple matrices or sparse matrices (SparseM library)
##
##    dirGraph   --  the graph to be transormed
##    useAlgo    --  'sparse' or 'normal'
##    useWeights --  if weights should be used (if useAlgo = 'normal' that
##                   the weigths are used anyway)                

reverseArch <- function(dirGraph,
                        useAlgo = 'sparse',
                        useWeights = TRUE) {

  if(edgemode(dirGraph) != 'directed')
    stop('The graph is not directed. Nothing to do')

  if(numEdges(dirGraph) == 0) {
    cat('Nothing to do:\n')
    return(dirGraph)
  }
     
  if(useAlgo == "sparse") {
    nodNames <- nodes(dirGraph)
    return(sparseM2Graph(t(graph2SparseM(dirGraph, useWeights)),
                         nodNames, edgemode = "directed"))
  }

  ## To slow for large graphs.  Probably is better to use only the SparseM package
  return(as(t(as(dirGraph, 'matrix')), 'graphNEL'))
}



############################   getGraphRoot   ############################
## find the root(roots) of the DAG

getGraphRoot <- function(dag,
                         leafs2root = TRUE) {

  if(!leafs2root)
    dag <- reverseArch(dag)
  
  aux <- sapply(adj(dag, nodes(dag)), length)
  return(names(which(aux == 0)))
}



########################   buildGOgraph.topology   ########################
## The structure of the GO graph is build recursivly.

buildGOgraph.topology <- function(knownNodes, whichOnto = "BP") {

  ## first build the lookUp table for the GO terms
  nodeLookUp <- new.env(hash = TRUE, parent = emptyenv())
  #GOOTerm <- get(paste('GO', whichOnto, 'Term', sep = ''), mode = 'environment')
  GOOTerm <- get('Term', mode = 'environment')
  
  ## warping functions for a easier acces to the lookUp table
  isNodeInDAG <- function(node) {
    return(exists(node, envir = nodeLookUp, mode = 'logical', inherits = FALSE))
  }
  setNodeInDAG <- function(node) {
    assign(node, TRUE, envir = nodeLookUp)
  }

  ## get the parents mapping
  GOParents <- get('ONTPARENTS')

  ## get the root of the ontology
  GENE.ONTO.ROOT <- as.character(names(as.list(GOParents)[match('NA',as.list(GOParents))]))

  ## we read all the database once and the access the list
  adjLookUP <- as.list(GOParents)
  
  ## we use an environment of environments to store edges: (this way is faster)
  ## in the end we will coerce it to a list of list and build a graphNEL obj. 
  edgeEnv <- new.env(hash = TRUE, parent = emptyenv())
  
  ## add the arc (u --> v) to edgeEnv of type :
  ##    0 for a is_a relation
  ##    1 for a part_of relation
  envAddEdge <- function(u, v, type) {
    assign(v, switch(type, isa = 0, partof = 1, -1), envir = get(u, envir = edgeEnv))
  }
  
  ## recursivly build the induced graph starting from one node
  buildInducedGraph <- function(node) {
    ## if we have visited the node, there is nothing to do
    if(isNodeInDAG(node))
      return(1)
    
    ## we put the node in the graph and we get his parents
    setNodeInDAG(node)    # we visit the node
    assign(node, new.env(hash = TRUE, parent = emptyenv()), envir = edgeEnv) # adj list
  
    if(node == GENE.ONTO.ROOT) 
      return(2)

    adjNodes <- adjLookUP[[node]]

    ## debuging point! should not happen!
    if(length(adjNodes) == 0)
      cat('\n There are no adj nodes for node: ', node, '\n')
          
    for(i in 1:length(adjNodes)) {
      x <- as.character(adjNodes[i])
      envAddEdge(node, x, names(adjNodes[i]))
      buildInducedGraph(x)
    }

    return(0)
  }

  ## we start from the most specific nodes
  lapply(knownNodes, buildInducedGraph)
  
  ## now we must transform our env into a Graph structure
  ## for now we use lapply, later we can do it with eapply
  .graphNodes <- ls(edgeEnv)
  .edgeList <- eapply(edgeEnv,
                      function(adjEnv) {
                        aux <- as.list(adjEnv)
                        return(list(edges = match(names(aux), .graphNodes),
                                    weights = as.numeric(aux)))
                      })
  
  ## now we can build the graphNEL object
  GOgraph.topo <- new('graphNEL',
                      nodes = .graphNodes,
                      edgeL = .edgeList,
                      edgemode = 'directed')
  
  return(GOgraph.topo)
}



############################   induceGraph   ############################
## Given a GO term (or a list of GO terms) this function is returning
## the subgraph induced by node.
#' @name inducedGraph
#' @title The subgraph induced by a set of nodes.
#' @aliases inducedGraph nodesInInducedGraph
#' @description 
#'   Given a set of nodes (GO terms) this function is returning
#'   the subgraph containing these nodes and their ancestors.
#' 
#' 
#' @usage
#' inducedGraph(dag, startNodes)
#' nodesInInducedGraph(dag, startNodes)
#' 
#'  @param dag An object of class \code{graphNEL} containing a directed graph.
#'  @param startNodes A character vector giving the starting nodes.
#' 
#' 
#' @return
#'   An object of class \code{\link[graph]{graphNEL-class}} is returned.
#' 
#' @author Adrian Alexa#' 
#' @seealso
#'   \code{\link{topONTdata-class}},
#'   \code{\link{reverseArch}},
#' 
#' @examples
#' data(ONTdata)
#' 
#' ## the GO graph
#' g <- graph(GOdata)
#' g
#' 
#' ## select 10 random nodes
#' sn <- sample(nodes(g), 10)
#' 
#' 
#' ## the subgraph induced by these nodes
#' sg <- inducedGraph(g, sn)
#' sg
#' 
#' @keywords graphs
inducedGraph <- function(dag,
                         startNodes) {

  return(subGraph(nodesInInducedGraph(dag, startNodes), dag))
}


############################   nodesInInducedGraph   ############################
## Given a GO term (or a list of GO terms) this function is returning
## the nodes in the subgraph induced by node.

nodesInInducedGraph <- function(dag,
                                 startNodes) {
  
  ## build a lookUp table with the nodes in the graph
  nodeLookUp <- new.env(hash = TRUE, parent = emptyenv())
  
  nodesDAG <- dag@nodes

  ## recursivly build the list of induced nodes
  buildInducedGraph <- function(node) {
    ## if we have visited the node, there is nothing to do
    if(exists(node, envir = nodeLookUp, mode = 'logical', inherits = FALSE))
      return(1)
    
    ## we put the node in the graph and we get his parents
    assign(node, TRUE, envir = nodeLookUp)
                                            
    adjNodes <- nodesDAG[dag@edgeL[[node]]$edges]
    
    if(length(adjNodes) == 0)
      return(2)
    
    for(i in 1:length(adjNodes))
      buildInducedGraph(adjNodes[i])
    
    return(0)
  }

  ## we start from the specified nodes
  lapply(startNodes, buildInducedGraph)

  return(ls(nodeLookUp))
}


nodesInInducedGraph2 <- function(dag,
                                startNodes) {
  
  ## build a lookUp table with the nodes in the graph
  nodeLookUp <- new.env(hash = TRUE, parent = emptyenv())

  edgesDAG <- edges(dag)

  ## warping functions for a easier acces to the lookUp table
  isNodeInDAG <- function(node) {
    return(exists(node, envir = nodeLookUp, mode = 'logical', inherits = FALSE))
  }
  setNodeInDAG <- function(node) {
    assign(node, TRUE, envir = nodeLookUp)
  }

  ## recursivly build the list of induced nodes
  buildInducedGraph <- function(node) {
    ## if we have visited the node, there is nothing to do
    if(isNodeInDAG(node))
      return(1)
    
    ## we put the node in the graph and we get his parents
    setNodeInDAG(node)  
                                        
    adjNodes <- edgesDAG[[node]]

    if(length(adjNodes) == 0)
      return(2)
    
    for(i in 1:length(adjNodes))
      buildInducedGraph(adjNodes[i])
    
    return(0)
  }

  ## we start from the specified nodes
  lapply(startNodes, buildInducedGraph)

  return(ls(nodeLookUp))
}


############################    buildLevels    ############################
## This function take the GOgraph and constructs a named vector which contain
## the level on which a node is. The root has level 1.
## The leafs2root parameter tell if the graph has edges directed from
## the leafs to the root, or vice-versa

## we return two environments:
## level2nodes   --   the key is the level no. and as value we have the
##                    nodes on that level
## nodes2level   --   as key we have the GO term and as value we have the
##                    level on which that node is

buildLevels <- function(dag,
                        root = NULL,
                        leafs2root = TRUE) {
  
  ## we need the root of the tree
  if(is.null(root))
    root <- getGraphRoot(dag, leafs2root)

  if(length(root) > 1)
    warning(paste('The graph has:', length(root), 'roots', sep = ' ')) 
  
  ## we need the the graph with the edges from the root to the leafs
  if(leafs2root == TRUE)
    dag <- reverseArch(dag)
    
  nodes2level <- new.env(hash = TRUE, parent = emptyenv())

  queue <- as.character(root)
  level <- 1
  adjList <- adj(dag, nodes(dag))
  
  while(length(queue)) {
    ## assign the curent level
    multiassign(queue, rep(level, length(queue)), envir = nodes2level)

    ## move to the next level
    level <- level + 1
    queue <- unique(unlist(adjList[queue], use.names = FALSE))
  }

  ## revert the node2level
  nl <- unlist(as.list(nodes2level))
  f.index <- rep(sort(unique(nl)), table(nl))
  level2nodes <- split(names(sort(nl)), f.index)

  return(list(nodes2level = nodes2level,
              level2nodes = list2env(level2nodes),
              noOfLevels = level - 1,
              noOfNodes = length(nodes2level)))
}
  

############################    buildLevels    ############################
## This function return the number of levels of a DAG

getNoOfLevels <- function(graphLevels) {
  return(max(as.integer(ls(graphLevels$level2nodes))))
}



###########################   mapGenes2GOgraph   ###########################
## This function builds for each node a vector containing all the genes/probes
## that can be annotated to that node.
## It starts with the nodes on the lowest level, and then pushes their genes
## to the parents/ancestors

## it returns the graph for which the attribute of each node contains a mapping
## of the genes/probes

mapGenes2GOgraph2 <- function(dag,
                             mostSpecificGOs,
                             nodeLevel = buildLevels(dag, leafs2root = TRUE)) {
  #browser()
  #init gene weight if not given
  if(names(mostSpecificGOs[[1]][1])=='NA' | is.null(names(mostSpecificGOs[[1]][1]))){
    mostSpecificGOs<-lapply(mostSpecificGOs,function(x){
      names(x)=rep(1,length(x))
      x})
  }
  
  allNodes <- nodes(dag)
  ## just in case .....
  if((ln.allNodes <- length(allNodes)) != nodeLevel$noOfNodes)
    stop('nodeLevel is corrupt')

  geneTerms <- new.env(hash = TRUE, parent = emptyenv())
  nn <- names(mostSpecificGOs)
  lapply(allNodes,
         function(x) {
           e <- new.env(hash = TRUE, parent = emptyenv())

           if(x %in% nn)
             multiassign(mostSpecificGOs[[x]], names(mostSpecificGOs[[x]]), envir = e)
             #multiassign(mostSpecificGOs[[x]], rep(TRUE, length(mostSpecificGOs[[x]])), envir = e)
             
           assign(x, e, envir = geneTerms)
         })
  
  
  ## get the levels list
  levelsLookUp <- nodeLevel$level2nodes
  noOfLevels <- nodeLevel$noOfLevels
  
  .aggregateScore<-function(currentE,destE){
    g1<-unlist(as.list.environment(currentE))
    g2<-unlist(as.list.environment(destE))
    ##
    for(i in names(g1)){
      if(i %in% names(g2)){
        g2[i]<-sum(as.numeric(g2[i]),as.numeric(g1[i]))
      }else{
        g2[i]<-as.numeric(g1[i])
      }
    }
    multiassign(names(g2), unname(g2), envir = destE)
  }
  
  for(i in noOfLevels:1) {
    currentNodes <- get(as.character(i), envir = levelsLookUp, mode = 'character')

    ## get all the adjacent nodes (teoreticaly nodes from level i - 1)
    adjList <- adj(dag, currentNodes)

    ## push the genes from level i to level i - 1
    lapply(currentNodes,
           function(node) {
             ## get the genes from this node
             genesID <- ls(get(node, envir = geneTerms, mode = 'environment'))

             ## debug option, just in case something goes wrong
             if(length(genesID) == 0)
               print(i)
            
             ## for each adiacent node mapp the genesID to them
             lapply(adjList[[node]],
                    function(destNode) {
                      .aggregateScore(get(node, envir = geneTerms, mode = 'environment'),get(destNode, envir = geneTerms, mode = 'environment'))
                      return(NULL)
                    })
             return(NULL)
           })
  }

  ## Assign for each node in the graph the coresponding environment
  nodeDataDefaults(dag, attr = "genes") <- emptyenv()
  nodeData(dag, allNodes, attr = "genes") <- as.list(geneTerms)[allNodes]
  
  return(dag)
}

mapGenes2GOgraph <- function(dag,
                             mostSpecificGOs,
                             nodeLevel = buildLevels(dag, leafs2root = TRUE)) {
  
  allNodes <- nodes(dag)
  ## just in case .....
  if((ln.allNodes <- length(allNodes)) != nodeLevel$noOfNodes)
    stop('nodeLevel is corrupt')
  #browser()
  geneTerms <- new.env(hash = TRUE, parent = emptyenv())
  nn <- names(mostSpecificGOs)
  lapply(allNodes,
         function(x) {
           e <- new.env(hash = TRUE, parent = emptyenv())
           
           if(x %in% nn)
            multiassign(mostSpecificGOs[[x]], rep(TRUE, length(mostSpecificGOs[[x]])), envir = e)
           
           assign(x, e, envir = geneTerms)
         })
  
  
  ## get the levels list
  levelsLookUp <- nodeLevel$level2nodes
  noOfLevels <- nodeLevel$noOfLevels
  
  for(i in noOfLevels:1) {
    currentNodes <- get(as.character(i), envir = levelsLookUp, mode = 'character')
    
    ## get all the adjacent nodes (teoreticaly nodes from level i - 1)
    adjList <- adj(dag, currentNodes)
    
    ## push the genes from level i to level i - 1
    lapply(currentNodes,
           function(node) {
             ## get the genes from this node
             genesID <- ls(get(node, envir = geneTerms, mode = 'environment'))
             
             ## debug option, just in case something goes wrong
             if(length(genesID) == 0)
               print(i)
             
             ## for each adiacent node mapp the genesID to them
             lapply(adjList[[node]],
                    function(destNode) {
                      destEnv <- get(destNode, envir = geneTerms, mode = 'environment')
                      multiassign(genesID, rep(FALSE, length(genesID)), envir = destEnv)
                      return(NULL)
                    })
             return(NULL)
           })
  }
  
  ## Assign for each node in the graph the coresponding environment
  nodeDataDefaults(dag, attr = "genes") <- emptyenv()
  nodeData(dag, allNodes, attr = "genes") <- as.list(geneTerms)[allNodes]
  
  return(dag)
}




######################################################################
########################### misc functions ###########################

## return for each specified node the information stored to attr
.getFromNode <- function(g, attr, nNames) {

  x <- nNames %in% nodes(g)
  if(!all(x)) {
    warning("Nodes not present in the graph:", nNames[!x])
    nNames <- nNames[x]
  }
  
  retValue <- nodeData(g, nNames, attr = attr)
  
  if(length(retValue) == 1)
    return(retValue[[1]])
  
  return(retValue[nNames])
}

.genesInNode <- function(g, nNames,score=F) {

  x <- nNames %in% nodes(g)
  if(!all(x)) {
    warning("Nodes not present in the graph:", nNames[!x])
    nNames <- nNames[x]
  }

  if(score)
    retValue <- lapply(nodeData(g, nNames, attr = "genes"), function(x){
      tmp<-unlist(as.list.environment(x))
      genes<-names(tmp)
      names(genes)<-unname(tmp)
      genes
    })
  else
    retValue <- lapply(nodeData(g, nNames, attr = "genes"), ls)

  ##  if(length(retValue) == 1)
  ##   return(retValue[[1]])

  return(retValue[nNames])
}

## return for each specified node the number of mapped genes
.countsInNode <- function(g, nNames) {

  x <- nNames %in% nodes(g)
  if(!all(x)) {
    warning("Nodes not present in the graph:", nNames[!x])
    nNames <- nNames[x]
  }

  retValue <- sapply(nodeData(g, nNames, attr = "genes"), length)
  
  return(retValue[nNames])
}

.writeToNodes <- function(g, attr, nodeList) {
  
  ## Assign for each node in the graph the coresponding value from nodeList
  nodeDataDefaults(g, attr = attr) <- vector(class(nodeList))
  nodeData(g, names(nodeList), attr = attr) <- nodeList
  
  return(g)
}


.get.all.node.annotation<-function(ONTdata,score=F){
  .genesInNode(ONTdata@graph,nodes(ONTdata@graph),score = score)
}
