\name{whop.go.connect}
\alias{whop.go.connect}
\title{
Establish a connection to the AmiGO database servers
}
\description{
Establish a connection to the AmiGO database servers or an arbitrary one with
the same database schema as the AmiGO DB.
}
\usage{
whop.go.connect(althost = NA, altport = NA, altuser = NA, altpass = NA,
					altdb = NA, altdbdrivername=NA, dbdrvpkgnam=NA)
}
\arguments{
  \item{althost}{Optional override for the hostname of the database server; default "mysql.ebi.ac.uk"}
  \item{altport}{Optional override for the port to connect to on the database server; default 4085}
  \item{altuser}{Optional override for the username to authenticate with; default "go_select" }
  \item{altpass}{Optional override for the password to authenticate with; default "amigo"}
  \item{altdb}{Optional override for the database name to connect to;default "go_latest"}
  \item{altdbdrivername}{Optional override for the DBMS driver to use;default "MySQL"}
  \item{dbdrvpkgnam}{Optional hint which R package provides the DBMS driver (e.g. " RMySQL" for the MySQL DBMS driver) }
}
%\details{}
\value{
Success status
}
\references{
AmiGO database
}
\author{
Ulrich Wittelsbuerger
}
