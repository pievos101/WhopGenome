\name{fai_query4}
\alias{fai_query4}
\alias{fai_query5}
\title{
Extract a part of a FASTA sequence.
}
\description{
Return a part of the a FASTA sequence.
}
\usage{
fai_query4( faifh, sequencename, beginpos, endpos )
fai_query5( faifh, sequencename, beginpos, endpos, resultstring )
}
\arguments{
\item{faifh}{FAIhandle as returned by fai_open}
\item{sequencename}{Identifier of a sequence in the fasta file}
\item{beginpos}{Start position of the subsequence to extract}
\item{endpos}{End position of the subsequence to extract}
\item{resultstring}{Variable to store the results into}
}
\details{
Note:
	the fai_query3 and fai_query5 methods are DEPRECATED : to be as fast as possible, they modified a given variable's contents (resultstring) which
	will cause issues in R's internals!
Use .Call("FAI_query4", faifh, sequencename, beginpos, endpos ) to eliminate the overhead of using the R wrapper function.
Use this function in combination with a
	while( ( seq = fai_query4(F,region) ) != FALSE )
if you need to loop.
(This exploits the fact that only the string "FALSE" has a boolean value of FALSE, all others have a boolean value of TRUE.)
}
\value{
A string containing the (sub-)sequence, FALSE if it fails.
}
%XXXXXXXXXXXXXXX
\examples{

##
##	Example :
##
faifile  <- system.file("extdata", "ex.fasta", package = "WhopGenome")
faifh <- fai_open( faifile )
stopifnot( !is.null(faifh) )
result = fai_query4( faifh , "1", 9 , 20 )
if( result != FALSE )
{
	print( result )
}
fai_close( faifh )

}
\author{
Ulrich Wittelsbuerger
}
\seealso{
fai_open
}
