\name{ONTSYNONYM}
\alias{ONTSYNONYM}
\title{Map from ontology synonyms to ontology terms}
\description{
  ONTSYNONYM is an R object that provides mapping
  from ontology synonyms to ontology terms
}
\details{
  Mappings were based on data provided by: Disease Ontology
  With a date stamp from the source of: 20140417 (sub_version 2653)


}

\examples{
    x <- ONTSYNONYM
    sample(x, 3)
    # ontology ID "DOID:8757" has a lot of synonyms
    ONTTERM[["DOID:8757"]]
    #  ID "DOID:9134" is a synonym of  ID "DOID:8757"
    #ID(ONTSYNONYM[["DOID:9134"]])
}
\keyword{datasets}

