\name{ONTOBSOLETE}
\alias{ONTOBSOLETE}
\title{Annotation of ontology identifiers by terms defined and their status are obsolete}
\description{
This is an R object mapping ontology identifiers to the specific
terms in defined and their definition are obsolete
}
\details{
  All the obsolete ontology terms that are collected in this index will no longer exist 
  in other mapping objects.

}

\examples{
    # Convert the object to a list
    xx <- as.list(ONTOBSOLETE)
    if(length(xx) > 0){
        # Get the TERMS for the first elent of xx
        #ID(xx[[1]])
        Term(xx[[1]])
    }
}
\keyword{datasets}

