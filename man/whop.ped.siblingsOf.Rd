\name{whop.ped.siblingsOf}
\alias{whop.ped.siblingsOf}
\title{
Return list of siblings
}
\description{
From the dataset 'p', all individuals which list at least one entry in 'lis' as mother or father are returned
}
\usage{
whop.ped.siblingsOf(p, lis)
}
\arguments{
  \item{p}{The pedigree dataset}
  \item{lis}{One or more individual identifiers from the dataset}
}
\details{
All entries which list one of the individuals in parameter 'lis' as either mother or father are returned.
}
\value{
Table of rows from the pedigree
}
\author{
Ulrich Wittelsbuerger
}
