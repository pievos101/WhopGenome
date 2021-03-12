\name{vcf_getnumcontigs}
\alias{vcf_getnumcontigs}
\title{
Get the number of different contigs/chromosomes stored in the file
}
\description{
Get the number of different contigs/chromosomes stored in the file
}
\usage{
vcf_getnumcontigs(vcff)
}
\details{
Use .Call("VCF_getNumContig", vcff ) to eliminate the overhead of using the R wrapper function.
}
\arguments{
  \item{vcff}{VCF file handle}
}
\value{
The number of different contigs/chromosomes stored in the file.
}
\examples{
	vcffile <- vcf_open( system.file( "extdata", "ex.vcf.gz" , package="WhopGenome") )
	vcf_getnumcontigs( vcffile )
	# [1] 1
}
\author{
Ulrich Wittelsbuerger
}
