\name{vcf_setregion}
\alias{vcf_setregion}
\title{
Set region from which to return genome variation data.
}
\description{
Set region from which to return genome variation data.
}
\usage{
vcf_setregion( vcffh, tid, from=NA, to=NA )
}
\arguments{
  \item{vcffh}{VCF file handle}
  \item{tid}{Either a chromosome identifier (from and to MUST be specified) or a region string (rendering from and to unnecessary)}
  \item{from}{Start position of the region from which to return data, if str is a chromosome identifier}
  \item{to}{End position of the region from which to return data, if str is a chromosome identifier}
}
\details{
Parameter 'regionstr' is of the form "chr:beg-end", e.g. "1:102910-210030" for chromosome 1, positions >= 102910 and <= 210030.
Use .Call("VCF_setRegion", vcffh, chromosomeid, from, to ) to eliminate the overhead of using the R wrapper function.
}
\value{
TRUE or FALSE , whether the call succeeded or not.
}
\author{
Ulrich Wittelsbuerger
}
\examples{
	##
	##	Example:
	##
	vcffile <- vcf_open( system.file( "extdata" , "ex.vcf.gz" , package="WhopGenome" ) )
	
	vcf_setregion(vcffile, "Y", 1, 100000 )
	vcf_readLineVec( vcffile )
}
\seealso{
\code{\link{vcf_open}}
\code{\link{vcf_getregion}}
}
