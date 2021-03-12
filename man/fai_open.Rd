\name{fai_open}
\alias{fai_open}
\title{
Open a faidx-indexed FASTA file
}
\description{
Opens a FASTA file that has an associated .fai index file
}
\usage{
fai_open( filename )
}
\arguments{
  \item{filename}{File name of the FASTA file. A file filename.fai is expected to reside in the same path.}
}
\details{
Use .Call("FAI_open", filename ) to eliminate the slight overhead of using the R wrapper function.
}
\value{
Returns a FAIhandle that is required for fai_query3, fai_query5, fai_close
}
%XXXXXXXXXXXXXXX
\examples{

##
##	Example :
##
faifile  <- system.file("extdata", "ex.fasta", package = "WhopGenome")
faifh <- fai_open( faifile )
stopifnot( !is.null(faifh) )
}
%XXXXXXXXXXXXXXX
\author{
Ulrich Wittelsbuerger
}
\seealso{
fai_reopen, fai_query3, fai_query5
}
