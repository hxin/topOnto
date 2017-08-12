\name{ONTANCESTOR}
\alias{ONTANCESTOR}
\title{Annotation of ontology Identifiers to their Ancestors} 
\description{
  This data set describes associations between ontology
  terms and their ancestor  terms, based on the directed acyclic
  graph (DAG). The format is an R object mapping the terms to all ancestor terms, where an
  ancestor term is a more general term that precedes
  the given term in the DAG (in other words, the parents, and all
  their parents, etc.).
}
\details{
  Each term is mapped to a vector of ancestor terms.


}

\examples{
  # Convert the object to a list
  xx <- as.list(ONTANCESTOR)
  # Remove Terms that do not have any ancestor
  xx <- xx[!is.na(xx)]
  if(length(xx) > 0){
    # Get the ancestor IDs for the first two elents of xx
    ids <- xx[1:2]
  }
  
}
\keyword{datasets}

