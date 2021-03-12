\name{vcf_valid}
\alias{vcf_valid}
\title{
Returns whether a VCF file handle is valid and usable.
}
\description{
Returns whether a VCF file handle is valid and usable.
}
\usage{
vcf_valid(vcffh)
}
\arguments{
  \item{vcffh}{VCF handle}
}
\value{
TRUE or FALSE
}
\examples{
	vcffile <- vcf_open( system.file( "extdata", "ex.vcf.gz" , package="WhopGenome") )
	vcf_valid( vcffile )
}
\author{
Ulrich Wittelsbuerger
}
