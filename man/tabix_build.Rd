\name{tabix_build}
\alias{tabix_build}
\title{
Build a tabix index file for fast access to tab-separated-value formatted files.
}
\description{
Given a pre-sorted and compressed file in a compatible tab-separated-columns format,
create a Tabix index file to perform fast queries on regions of data.
}
\usage{
tabix_build( filename, sc, bc, ec, meta, lineskip )
}
\arguments{
\item{filename}{Name of file to create index for}
\item{sc}{Number of sequence column}
\item{bc}{Number of start column}
\item{ec}{Number of end column}
\item{meta}{Symbol used to begin comment/meta-information lines}
\item{lineskip}{Number of lines to skip from the top}
}
\details{
Tabix is a tool that has been developed to quickly retrieve data on an arbitrary chromosomal region from files that store their data in tab-separated columns,
such as VCF, BED, GFF and SAM. As long as there is a column for named groups (e.g. chromosomes) and another column giving a numerical order (e.g. chromosomal position),
it can be used for other data as well.
As a required preprocessing step, it creates an index file for a file which has been sorted by group names (e.g. chromosome) and location as well as gzip/bgzf-compressed.
After sorting, compressing and indexing, specific portions of such a file can be very efficiently retrieved, e.g. using the other tabix_XXX functions.
}
\value{
TRUE or FALSE.
}
%XXXXXXXXXXXXXXX
\examples{

##
##	Example :
##

gfffile  <- system.file("extdata", "ex.gff3", package = "WhopGenome" )
gfffile

gffbasename <- tempfile()
file.copy( from=gfffile, to=gffbasename )
gffgzfile <- paste( sep="", gffbasename, ".gz" )
gffgzfile

##
##
gffindexfile <- paste( sep="", gffgzfile, ".tbi" )
gffindexfile
stopifnot( ! file.exists( gffindexfile ) )
print( "Index file does not exist yet!" )

###
###	compress GFF file
###
bgzf_compress( gffbasename , gffgzfile )
stopifnot( file.exists( gffgzfile ) )
###
###	build index
###
tabix_build( filename = gffgzfile,
			 sc = as.integer(1),
			 bc = as.integer(2),
			 ec = as.integer(3),
			 meta = "#",
			 lineskip = as.integer(0)
			)
# [1] TRUE
stopifnot( file.exists( gffindexfile ) )
print( "Index file has been built" )
#
gffh <- tabix_open( gffgzfile )
gffh
}
\author{
Ulrich Wittelsbuerger
}
\seealso{
tabix_open, tabix_setregion, tabix_read
}
