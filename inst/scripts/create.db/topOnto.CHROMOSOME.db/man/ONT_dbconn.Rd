
\name{ONT_dbconn}

\alias{ONT_dbconn}
\alias{ONT_dbfile}
\alias{ONT_dbschema}
\alias{ONT_dbInfo}


\title{Collect information about the package annotation DB}

\description{
  Some convenience functions for getting a connection object to (or collecting
  information about) the package annotation DB.
}

\usage{
  ONT_dbconn()
  ONT_dbfile()
%  ONT_dbschema(file="", show.indices=FALSE)
  ONT_dbschema()
  ONT_dbInfo()
}

%\arguments{
%  \item{file}{
%    A connection, or a character string naming the file to print to (see
%    the \code{file} argument of the \code{\link[base]{cat}} function for
%    the details).
%  }
%  \item{show.indices}{
%    The CREATE INDEX statements are not shown by default.
%    Use \code{show.indices=TRUE} to get them.
%  }
%}

\details{
  \code{ONT_dbconn} returns a connection object to the
  package annotation DB.  IMPORTANT: Don't call
  \code{\link[DBI:dbConnect]{dbDisconnect}} on the connection object
  returned by \code{ONT_dbconn} or you will break all the
  \code{\link[AnnotationDbi:AnnDbObj-class]{AnnDbObj}} objects defined
  in this package!

  \code{ONT_dbfile} returns the path (character string) to the
  package annotation DB (this is an SQLite file).

  \code{ONT_dbschema} prints the schema definition of the
  package annotation DB.

  \code{ONT_dbInfo} prints other information about the package
  annotation DB.
}

\value{
  \code{ONT_dbconn}: a DBIConnection object representing an
  open connection to the package annotation DB.

  \code{ONT_dbfile}: a character string with the path to the
  package annotation DB.

  \code{ONT_dbschema}: none (invisible \code{NULL}).

  \code{ONT_dbInfo}: none (invisible \code{NULL}).
}

\seealso{
  \code{\link[DBI:dbSendQuery]{dbGetQuery}},
  \code{\link[DBI:dbConnect]{dbConnect}},
  \code{\link[AnnotationDbi:AnnDbObj-class]{dbconn}},
  \code{\link[AnnotationDbi:AnnDbObj-class]{dbfile}},
  \code{\link[AnnotationDbi:AnnDbObj-class]{dbschema}},
  \code{\link[AnnotationDbi:AnnDbObj-class]{dbInfo}}
}

\examples{
  ## Count the number of rows in the "term" table:
  dbGetQuery(ONT_dbconn(), "SELECT COUNT(*) FROM term")

  ## The connection object returned by ONT_dbconn() was
  ## created with:
  dbConnect(SQLite(), dbname=ONT_dbfile(), cache_size=64000,
  synchronous='off')

  #ONT_dbschema()

#  ONT_dbInfo()
}

\keyword{utilities}
\keyword{datasets}

