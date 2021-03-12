\name{fai_query2}
\alias{fai_query2}
\alias{fai_query3}
\title{
Extract a part of a FASTA sequence.
}
\description{
Return a part of a FASTA sequence.
}
\usage{
fai_query2( faifh, regionstring )
fai_query3( faifh, regionstring , resultstring )
}
\arguments{
\item{faifh}{FAIhandle as returned by fai_open}
\item{regionstring}{String of the form sequencename:beginpos-endpos e.g. "MTRR1mouse:20-40" specifying the sequence and region}
\item{resultstring}{String variable to store results into}
}
\details{
Note:
	the fai_query3 and fai_query5 methods are DEPRECATED : to be as fast as possible, they modified a given variable's contents (resultstring) which
	will cause issues in R's internals!

Use .Call("FAI_query2", faifh, regionstring ) to eliminate the overhead of using the R wrapper function.
Use this function in combination with a
	while( ( seq = fai_query2(F,region) ) != FALSE )
if you need to loop. Only the string "FALSE" has a boolean value of FALSE, all others have a boolean value of TRUE.
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
result = fai_query2( faifh , "1:100-200" )
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
