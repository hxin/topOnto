\name{ONTPARENTS}
\alias{ONTPARENTS}
\title{Annotation of ontology Identifiers to their Parents} 
\description{
  This data set describes associations between ontology
  terms and their direct parent  terms, based on the directed acyclic
  graph (DAG). The format is an R
  object mapping the ontology  terms to all direct parent terms, where a
  direct parent term is a more general ontology term that immediately precedes
  the given ontology term in the DAG.
}
\details{
  Each ontology  term is mapped to a named vector of ontology  terms. The name
  associated with the parent term will be either \emph{isa}, \emph{partof},
  where \emph{isa} indicates that the child term is a more specific version
  of the parent, and  \emph{partof} indicate that the
  child term is a part of the parent. 



}

\examples{
  # Convert the object to a list
  xx <- as.list(ONTPARENTS)
  # Remove ontology IDs that do not have any parent
  xx <- xx[!is.na(xx)]
  if(length(xx) > 0){
     ids <- xx[[1]]
     # Find out the  terms for the first parent ID
     #ID(ONTTERM[[ids[1]]])
     Term(ONTTERM[[ids[1]]])
     Synonym(ONTTERM[[ids[1]]])
     Secondary(ONTTERM[[ids[1]]])
  }
}
\keyword{datasets}

