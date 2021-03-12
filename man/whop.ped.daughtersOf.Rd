\name{whop.ped.daughtersOf}
\alias{whop.ped.daughtersOf}
\title{
Return all daughters of a given individual from a pedigree dataset
}
\description{
All individuals which are female and have at least one of the given IDs as either mother or father
}
\usage{
whop.ped.daughtersOf(p, lis)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{The pedigree dataset}
  \item{lis}{One or more individual IDs}
}
\value{
Table of rows from the pedigree
}
\author{
Ulrich Wittelsbuerger
}
