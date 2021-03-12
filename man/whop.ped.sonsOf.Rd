\name{whop.ped.sonsOf}
\alias{whop.ped.sonsOf}
\title{
Returns all sons of the given individuals
}
\description{
All individuals in a pedigree data
}
\usage{
whop.ped.sonsOf(p, lis)
}
\arguments{
  \item{p}{The pedigree dataset to work on}
  \item{lis}{One or more individuals' identifiers from the dataset}
}
\details{
For each element in lis, finds all male individuals who refer to these elements as parent.
Essentially combines a whop.ped.males() with a whop.ped.siblingsOf() call.
}
\value{
Table of rows from the pedigree
}
\author{
Ulrich Wittelsbuerger
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
whop.ped.daughtersOf
}
