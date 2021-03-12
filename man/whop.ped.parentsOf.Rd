\name{whop.ped.parentsOf}
\alias{whop.ped.parentsOf}
\title{
Return the parents of individuals
}
\description{
Looks for all individuals which are listed as parents of certain other individuals.
}
\usage{
whop.ped.parentsOf(p, invids)
}
\arguments{
  \item{p}{The pedigree dataset}
  \item{invids}{One or more individuals' identifiers from the dataset}
}
\details{
All individuals which appear in the Maternal.ID and Paternal.ID columns of the given individuals.
}
\value{
Table of rows from the pedigree
}
\author{
Ulrich Wittelsbuerger
}
