\name{tabix_close}
\alias{tabix_close}
\title{
Close Tabix-indexed file
}
\description{
Close Tabix-indexed file
}
\usage{
tabix_close(tabfh)
}
\arguments{
  \item{tabfh}{Tabix file handle}
}
\value{
None.
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
