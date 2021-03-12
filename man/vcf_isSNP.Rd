\name{vcf_isSNP}
\alias{vcf_isSNP}
\title{
Determines whether the last vcf_parse-call returned a SNP (instead of InDel)
}
\description{
Determines whether the last vcf_parse/VCF_parse-call returned a SNP (instead of InDel)
}
\usage{
vcf_isSNP(vcff)
}
\arguments{
  \item{vcff}{VCF file handle}
}
\value{
None.
}
\details{
Use .Call("VCF_isSNP", vcff ) to eliminate the overhead of using the R wrapper function.
}
\author{
Ulrich Wittelsbuerger
}
\examples{
	vcffile <- vcf_open( system.file( "extdata", "ex.vcf.gz" , package="WhopGenome") )
	vcf_parseNextSNP( vcffile )
	vcf_getPos( vcffile )
	vcf_isSNP( vcffile )
}
\seealso{
vcf_isINDEL
}
