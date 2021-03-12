\name{vcf_getheaderline}
\alias{vcf_getheaderline}
\title{
Return one of the header lines of the VCF file
}
\description{
Return one of the header lines of the VCF file
}
\usage{
vcf_getheaderline(vcff, whichnum)
}
\arguments{
  \item{vcff}{VCF file handle}
  \item{whichnum}{Number of header line to retrieve}
}
\details{
Use .Call("VCF_getHeaderLine", vcff, whichnum ) to eliminate the overhead of using the R wrapper function.
}
\value{
A string containing the full header line.
}
\examples{
	vcffile <- vcf_open( system.file( "extdata", "ex.vcf.gz" , package="WhopGenome") )
	vcf_getheaderline( vcffile , as.integer(0) )
	vcf_getheaderline( vcffile , as.integer(1) )
}
\author{
Ulrich Wittelsbuerger
}
