\name{ONTTERM}
\alias{ONTTERM}
\title{Annotation of ontology Identifiers to ontology Terms}
\description{
  This data set gives mappings between ontology identifiers and their respective terms.
}
\details{
  Each ontology identifier is mapped to a \code{ONTTERMs} object that has 4 slots:
   ID: ontology Identifier;
   Term: The term for that ontology id;
   Secondary: Secondary terms that have been merged into this term;
   Synonym: other  ontology terms that are considered to be synonymous to the primary
      term attached to the ontology id
 
}

\examples{
    # Convert the object to a list
    FirstTenBimap <- ONTTERM[1:10] ##grab the 1st ten
    xx <- as.list(FirstTenBimap)
     if(length(xx) > 0){
        # Get the TERMS for the 2nd element of xx
        ID(xx[[2]])
        Term(xx[[2]])
        Synonym(xx[[2]])
        Secondary(xx[[2]])
    }
}
\keyword{datasets}

