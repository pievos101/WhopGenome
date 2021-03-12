\name{vcf_isINDEL}
\alias{vcf_isINDEL}
\title{
Determines whether the last vcf_parse-call returned a InDel (instead of SNP)
}
\description{
Returns TRUE if the last call to vcf_parse/VCF_parse returned an InDel.
}
\usage{
vcf_isINDEL(vcff)
}
\arguments{
  \item{vcff}{VCF file handle}
}
\details{
Use .Call("VCF_isInDel", vcff ) to eliminate the overhead of using the R wrapper function.
}
\value{
TRUE or FALSE.
}
\author{
Ulrich Wittelsbuerger
}
\examples{
	vcffile <- vcf_open( system.file( "extdata", "ex.vcf.gz" , package="WhopGenome") )
	vcf_parseNextSNP( vcffile )
	vcf_getPos( vcffile )
	vcf_isINDEL( vcffile )
}
\seealso{
vcf_isSNP
}
