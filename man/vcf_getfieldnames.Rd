\name{vcf_getfieldnames}
\alias{vcf_getfieldnames}
\title{
Return a vector with the field names used in the VCF file.
}
\description{
Return a vector with the field names used in the VCF file.
}
\usage{
vcf_getfieldnames(vcff)
}
\arguments{
  \item{vcff}{VCF file handle}
}
\details{
Use .Call("VCF_getFieldNames", vcff ) to eliminate the overhead of using the R wrapper function.
}
\value{
A vector of strings representing the field names present in the VCF file.
}
\examples{
	vcffile <- vcf_open( system.file( "extdata", "ex.vcf.gz" , package="WhopGenome") )
	vcf_getfieldnames( vcffile )
}
\author{
Ulrich Wittelsbuerger
}
