\name{tabix_getregion}
\alias{tabix_getregion}
\title{
Return the currently selected region of the given tabix file.
}
\description{
Return the currently selected region of the given tabix file. The resulting value does
not reflect the current read position inside that region, i.e. you cannot infer whether
there are any lines left for reading from that region.
}
\usage{
tabix_getregion( tabfh )
}
\arguments{
  \item{tabfh}{Tabix handle, once returned by tabix_open}
}
\details{
Use .Call("tabix_getRegion", tabfh ) to eliminate the slight overhead of using the R wrapper function.
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
tabix_getregion( gffh )
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
