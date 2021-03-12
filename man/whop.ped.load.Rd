\name{whop.ped.load}
\alias{whop.ped.load}
\title{
Load a pedigree dataset from a .PED file
}
\description{
Returns a table with the pedigree data contained in the file
}
\usage{
whop.ped.load(filename)
}
\arguments{
  \item{filename}{Name of the file containing the pedigree data}
}
\details{
Expects the given file to be of the PLINK .PED format, i.e. a file with tab-separated columns of which
the first few are required to be of a certain order.
}
\value{
Table with pedigree data
}
\references{
PLINK .PED
}
\author{
Ulrich Wittelsbuerger
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
whop.ped.save
}
