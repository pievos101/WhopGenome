\name{fai_close}
\alias{fai_close}
\title{
Closes a file previously opened with fai_open
}
\description{
Closes a file previously opened with fai_open
}
\usage{
fai_close( faifh )
}
\arguments{
  \item{faifh}{A FAIhandle as returned by fai_open}
}
\details{
Use .Call("FAI_close", faifh ) to eliminate the slight overhead of using the R wrapper function.
}
\value{
TRUE if call succeeds, FALSE if it fails.
}
%XXXXXXXXXXXXXXX
\examples{

##
##	Example :
##
faifile  <- system.file("extdata", "ex.fasta", package = "WhopGenome")
faifh <- fai_open( faifile )
stopifnot( !is.null(faifh) )
fai_close( faifh )

}
%XXXXXXXXXXXXXXX
\author{
Ulrich Wittelsbuerger
}
\seealso{
fai_open
}
