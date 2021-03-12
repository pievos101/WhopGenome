\name{tabix_open}
\alias{tabix_open}
\title{
Open Tabix-indexed file for subsequent access with other tabix_ methods
}
\description{
Open Tabix-indexed file for subsequent access with other tabix_ methods
}
\usage{
tabix_open(filename)
}
\arguments{
  \item{filename}{String, name of tabix-indexed file to open}
}
\details{
As filename, specify the data file, not the index file ending in .tbi!
}
\value{
Tabix file handle
}
%XXXXXXXXXXXXXXX
\examples{

##
##	Example :
##
gffgzfile  <- system.file("extdata", "ex.gff3.gz", package = "WhopGenome" )
gffh <- tabix_open( gffgzfile )
gffh
tabix_close( gffh )
gffh

}
\author{
Ulrich Wittelsbuerger
}
\seealso{
\code{\link{tabix_open}}
\code{\link{tabix_read}}
}
