\name{tabix_read}
\alias{tabix_read}
\alias{tabix_readraw}
\title{
Read a line from a tabix_open()'ed file
}
\description{
Read a line from a tabix_open()'ed file
}
\usage{
tabix_read( tabfh )
tabix_readraw( tabfh )
}
\arguments{
  \item{tabfh}{Tabix file handle as returned by tabix_open}
}
\value{
A line of data from the indexed data file. tabix_read splits the line up into its
fields and returns a vector. tabix_readraw returns the line as stored in the file.
}
\details{
Instead of tabix_readraw() you can use .Call("tabix_readLine", tabfh ) to eliminate the slight overhead of using the R wrapper function.
}
\author{
Ulrich Wittelsbuerger
}
%XXXXXXXXXXXXXXX
\examples{

##
##	Example : (NOT RUN)
##

print("Opening and reading")
gffgzfile  <- system.file("extdata", "ex.gff3.gz", package = "WhopGenome" )
if( file.exists(gffgzfile) )
{
 gffgzfile
 gffh <- tabix_open( gffgzfile )
 gffh
 stopifnot( !is.null(gffh) )
 tabix_read( gffh )
 tabix_close( gffh )
 gffh
}
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
tabix_open
}
