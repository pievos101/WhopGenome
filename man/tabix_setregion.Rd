\name{tabix_setregion}
\alias{tabix_setregion}
\title{
Reopen a Tabix-indexed file if the filehandle became invalid.
}
\description{
Reopen a Tabix-indexed file if the filehandle became invalid.
}
\usage{
tabix_setregion( tabfh, tid, beginpos, endpos )
}
\arguments{
  \item{tabfh}{Tabix handle, once returned by tabix_open}
  \item{tid}{A string naming one of the contig/chromosome identifiers stored in the Tabix indexed file}
  \item{beginpos}{Earliest position from which subsequent tabix_read/tabix_readraw calls return lines}
  \item{endpos}{Last position to return lines from with tabix_read/tabix_readraw}
}
\details{
Use .Call("tabix_setRegion", tabfh, tid, beginpos, endpos ) to eliminate the slight overhead of using the R wrapper function.
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
tabix_setregion( gffh, "ex.1", 1, 400 )
tabix_close( gffh )
gffh


}
\author{
Ulrich Wittelsbuerger
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
tabix_open
}
