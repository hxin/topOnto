\name{ONTTerms-class}
\docType{class}

\alias{class:ONTTerms}
\alias{ONTTerms-class}

\alias{ONTTerms}

\alias{initialize,ONTTerms-method}

\alias{ID}
\alias{ID,ONTTerms-method}
\alias{ID,ONTTermsAnnDbBimap-method}
\alias{ID,character-method}
\alias{Definition}
\alias{Definition,ONTTerms-method}
\alias{Definition,ONTTermsAnnDbBimap-method}
\alias{Definition,character-method}
\alias{Ontology}
\alias{Ontology,ONTTerms-method}
\alias{Ontology,ONTTermsAnnDbBimap-method}
\alias{Ontology,character-method}
\alias{Term}
\alias{Term,ONTTerms-method}
\alias{Term,ONTTermsAnnDbBimap-method}
\alias{Term,character-method}
\alias{Synonym}
\alias{Synonym,ONTTerms-method}
\alias{Synonym,ONTTermsAnnDbBimap-method}
\alias{Synonym,character-method}
\alias{Secondary}
\alias{Secondary,ONTTerms-method}
\alias{Secondary,ONTTermsAnnDbBimap-method}
\alias{Secondary,character-method}

\alias{show,ONTTerms-method}



\title{Class "ONTTerms"}

\description{A class to represent Disease Ontology nodes}

\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{ONTTerms(ID, term, synonym, secondary)}.
  ID, term are required.
}

\section{Slots}{
  \describe{
    \item{\code{ID}:}{Object of class \code{"character"} A character
      string for the ontology id of a primary node.}
    \item{\code{Definition}:}{Object of class \code{"character"} A character
      string for the ontology Definition of a primary node.}
    \item{\code{Ontology}:}{Object of class \code{"character"} A character
      string for the ontology name of a primary node.}
    \item{\code{Term}:}{Object of class \code{"character"} A
      character string that defines the role of gene product
      corresponding to the primary ontology id.}
    \item{\code{Synonym}:}{Object of class \code{"character"} other
      ontology terms that are considered to be synonymous to the primary
      term attached to the ontology id. Synonymous here can mean that the
      synonym is an exact synonym of the primary term, is related to the
      primary term, is broader than the primary term, is more precise
      than the primary term, or name is related to the term, but is not
      exact, broader or narrower.}
    \item{\code{Secondary}:}{Object of class \code{"character"} ontology ids
      that are secondary to the primary ontology id as results of merging ontology
      terms so that One ontology id becomes the primary ontology id and the rest
      become the secondary.}
  }
}

\section{Methods}{
  \describe{
    \item{ID}{\code{signature(object = "ONTTerms")}:
      The get method for slot ID.}
    \item{Term}{\code{signature(object = "ONTTerms")}:
      The get method for slot Term.}
    \item{Definition}{\code{signature(object = "ONTTerms")}:
      The get method for slot Term.}
    \item{Ontology}{\code{signature(object = "ONTTerms")}:
      The get method for slot Term.}
    \item{Synonym}{\code{signature(object = "ONTTerms")}:
      The get method for slot Synonym.}
    \item{Secondary}{\code{signature(object = "ONTTerms")}:
      The get method for slot Secondary.}
    \item{show}{\code{signature(x = "ONTTerms")}:
      The method for pretty print.}
  }
}


\note{ONTTerms objects are used to represent primary ontology nodes in the
   SQLite-based annotation data package topOnto.xxx.db}

\examples{
ONTnode <- new("ONTTerms", ID="DOID:1234567", Term="Test", Ontology="Test", Definition="Test")
\dontrun{##I want to show an ex., but don't want to require topOnto.xxx.db
  ID(ONTnode)
  Term(ONTnode)

  ##Or you can just use these methods on a ONTTermsAnnDbBimap

  #require(DO.db)
  FirstTenBimap <- ONTTERM[1:10] ##grab the 1st ten
  Term(FirstTenBimap)

  ##Or you can just use IDs directly
  ids = keys(FirstTenBimap)
  Term(ids)
}
}

\keyword{methods}
\keyword{classes}

