\name{ONTCHILDREN}
\alias{ONTCHILDREN}
\title{Annotation of ontology Identifiers to their Children} 
\description{
This data set describes associations between ontology 
  terms and their direct children  terms, based on the directed acyclic
  graph (DAG) defined by the Disease Ontology Consortium. The format is an R
  object mapping the ontology  terms to all direct children terms, where a
  direct child term is a more specific ontology term that is immediately preceded
  by the given ontology term in the DAG.
}
\details{
   Each ontologyterm is mapped to a vector of children ontology  terms.


}


\examples{
  # Convert the object to a list
  xx <- as.list(ONTCHILDREN)
  # Remove ontology IDs that do not have any children
  xx <- xx[!is.na(xx)]
  
  if(length(xx) > 0){
     # Get the parent ontology IDs for the first elents of xx
        ids <- xx[[1]]
        # Find out the terms for the first parent 
        #ID(ONTTERM[[ids[1]]])
        Term(ONTTERM[[ids[1]]])
        Synonym(ONTTERM[[ids[1]]])
        Secondary(ONTTERM[[ids[1]]])
  }
}
\keyword{datasets}

