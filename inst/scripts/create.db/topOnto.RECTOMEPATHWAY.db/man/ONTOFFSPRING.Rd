\name{ONTOFFSPRING}
\alias{ONTOFFSPRING}
\title{Annotation of ontology Identifiers to their Offspring} 
\description{
  This data set describes associations between ontology 
  terms and their offspring  terms, based on the directed acyclic
  graph (DAG). The format is an R  object mapping the ontology  terms to all offspring terms, where an
  ancestor term is a more specific ontology term that is preceded
  by the given ontology term in the DAG (in other words, the children and all
  their children, etc.).
}
\details{
  Each ontology  term is mapped to a vector of offspring ontology  terms.


}

\examples{
  # Convert the object to a list
  xx <- as.list(ONTOFFSPRING)
  # Remove ontology IDs that do not have any offspring
  xx <- xx[!is.na(xx)]
   if(length(xx) > 0){
    # Get the offspring ontology identifiers for the first two elents of xx
    identifiers <- xx[1:2]
  }
}
\keyword{datasets}

