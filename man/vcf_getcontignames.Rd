\name{vcf_getcontignames}
\alias{vcf_getcontignames}
\title{
Return the contig/chromosome identifiers used in the VCF file
}
\description{
Return the contig/chromosome identifiers used in the VCF file
}
\usage{
vcf_getcontignames(vcff)
}
\arguments{
  \item{vcff}{VCF file handle}
}
\details{
vcf_setregion for example requires one of these identifiers to be able to successfully select a region for extraction.
Use .Call("VCF_getContigNames", vcff ) to eliminate the overhead of using the R wrapper function.
}
\value{
Vector with contig and/or chromosome identifiers.
}
\author{
Ulrich Wittelsbuerger
}
\examples{
	vcffile <- vcf_open( system.file( "extdata", "ex.vcf.gz" , package="WhopGenome") )
	vcf_getcontignames( vcffile )
	# [1] "Y"
}
\seealso{
vcf_setregion
}
