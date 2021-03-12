\name{bgzf_compress}
\alias{bgzf_compress}
\title{
Compress file with bgzip
}
\description{
Write contents of file <infilename> bgzip-compressed to file named <outfilename>.
}
\usage{
bgzf_compress( infilename, outfilename )
}
\arguments{
\item{infilename}{Name of file to read data from for compression}
\item{outfilename}{Name of file to write compressed data to}
}
\details{
Compresses the file specified by <infilename> with the bgzip compression scheme,
as developed by Bob Handsaker and modified by Heng Li, and creates a compressed
file under the name given by <outfilename>.
}
\value{
TRUE if call succeeds, FALSE if it fails).
}
%XXXXXXXXXXXXXXX
\examples{

##
##	Example :
##
gfffile  <- system.file("data", "ex.gff3", package = "WhopGenome" )
gffgzfile <- paste( sep="", gfffile, ".gz" )
file.remove( gffgzfile )
bgzf_compress( gfffile , gffgzfile )
file.exists( gffgzfile )

}
%XXXXXXXXXXXXXXX
\author{
Ulrich Wittelsbuerger
}
