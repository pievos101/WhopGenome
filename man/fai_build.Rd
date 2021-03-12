\name{fai_build}
\alias{fai_build}
\title{
Build a .fai-index for the given FASTA file.
}
\description{
Build a .fai-index for the given FASTA file.
}
\usage{
fai_build( filename )
}
\arguments{
\item{filename}{Name of the FASTA file for which an index file should be built.}
}
\details{
Use .Call("FAI_build", filename ) to eliminate the overhead of using the R wrapper function.
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
faiindexfile <- paste( sep="", faifile, ".fai" )	# construct name of index file
file.remove( faiindexfile )				# remove existing index
fai_build( faifile )					# re-create index
stopifnot( file.exists(faiindexfile) )	# check whether index file exists
print( "All OK" )
}
%XXXXXXXXXXXXXXX
\author{
Ulrich Wittelsbuerger
}
\seealso{
fai_open
}
