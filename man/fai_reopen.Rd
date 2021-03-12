\name{fai_reopen}
\alias{fai_reopen}
\title{
Reopen a FAIhandle that has become stale.
}
\description{
Reopen a FAIhandle that has become stale, e.g. by restarting R or loading a workspace containing a FAIhandle variable.
}
\usage{
fai_reopen( faifh )
}
\arguments{
  \item{faifh}{A FAIhandle to a .fai-indexed FASTA file}
}
\details{
Use .Call("FAI_reopen", faifh ) to eliminate the slight overhead of using the R wrapper function.
}
\value{
TRUE if call succeeds, FALSE if it fails.
}
\author{
Ulrich Wittelsbuerger
}
\examples{
	##
	##	Example :
	##
	faifile  <- system.file("extdata", "ex.fasta", package = "WhopGenome")
	faifh <- fai_open( faifile )
	stopifnot( !is.null(faifh) )
	result <- fai_query4( faifh , "1", 100 , 200 )
	print( result )
	fai_close( faifh )
	fai_reopen( faifh )
	result <- fai_query4( faifh , "1", 100 , 110 )
	print( result )
}
\seealso{
fai_open
}
