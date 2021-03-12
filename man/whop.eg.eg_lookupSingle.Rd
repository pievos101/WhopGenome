\name{whop.eg.eg_lookupSingle}
\alias{whop.eg.eg_lookupSingle}
\title{
Return the first entry in an EG organism's data table for a given identifier
}
\description{
Return the first entry in an EG organism's data table for a given identifier
}
\usage{
whop.eg.eg_lookupSingle(id, subdbname, db)
}
\arguments{
  \item{id}{Identifiers to look for in subtable}
  \item{subdbname}{Organism annotation table name}
  \item{db}{Optional, organism database if not using default active one}
}
\value{
First entry with any of the given id(s) in the table
}
\author{
Ulrich Wittelsbuerger
}
