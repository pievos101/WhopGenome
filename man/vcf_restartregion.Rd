\name{vcf_restartregion}
\alias{vcf_restartregion}
\title{
Let subsequent read calls return from the start of the currently set region.
}
\description{
Once the read-functions reached the end of the previously set region, no more results are returned.
If, for example for a two-pass algorithm, the same region should be scanned again from the start,
this function is the key.
}
\usage{
vcf_restartregion(vcffh)
}
\arguments{
  \item{vcffh}{Handle of a VCF file, as returned by vcf_open()}
}
\details{
Alternative to calling vcf_setregion() with the same parameters again.
Use .Call("VCF_restartRegion", vcffh ) to eliminate the overhead of using the R wrapper function.
}
\value{
TRUE if the region could be rewound, FALSE if not.
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
	vcf_readLineVec( vcffile )

	vcf_restartregion( vcffile )

	vcf_readLineVec( vcffile )
	vcf_readLineVec( vcffile )
}
\seealso{
vcf_setregion, vcf_open
}
