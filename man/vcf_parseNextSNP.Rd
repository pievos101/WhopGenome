\name{vcf_parseNextSNP}
\alias{vcf_parseNextSNP}
\alias{vcf_parseNextLine}
\title{
Read until next SNP or next line and buffer it
}
\description{
Read until next SNP or next line and buffer it. Use the vcf_getXXX functions to access specific fields of the line
}
\usage{
vcf_parseNextSNP(vcffh)
vcf_parseNextLine(vcffh)
}
\arguments{
  \item{vcffh}{VCF file handle}
}
\details{
Use .Call("VCF_parseNextSNP", vcffh ) and .Call("VCF_parseNextLine", vcffh ) respectively, to eliminate the overhead of using the R wrapper function.
}
\value{
None.
}
\author{
Ulrich Wittelsbuerger
}
\examples{
	vcffile <- vcf_open( system.file( "extdata", "ex.vcf.gz" , package="WhopGenome") )
	vcf_parseNextSNP( vcffile )
	vcf_getPos( vcffile )
}
\seealso{
vcf_isSNP, vcf_open, vcf_getPos
}
