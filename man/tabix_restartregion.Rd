\name{tabix_restartregion}
\alias{tabix_restartregion}
\title{
Reset the currently selected region to the beginning.
}
\description{
Reset the currently selected region so that the next read call will return the
first line inside that region.
}
\usage{
tabix_restartregion( tabfh )
}
\arguments{
  \item{tabfh}{Tabix handle, once returned by tabix_open}
}
\details{
Use .Call("tabix_restartRegion", tabfh ) to eliminate the slight overhead of using the R wrapper function.
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
##
##
##
tabix_setregion( gffh, "ex.1", 1, 400 )
tabix_read( gffh )
tabix_read( gffh )
tabix_restartregion( gffh )
tabix_read( gffh )
tabix_read( gffh )
tabix_close( gffh )
gffh

}
\author{
Ulrich Wittelsbuerger
}
\seealso{
tabix_open
}
