\name{vcf_readLineDF}
\alias{vcf_readLineDF}
\title{
Read a line of data from the given VCF file and return it as a data frame
}
\description{
Read a line of data from the given VCF file and return it as a data frame
}
\usage{
vcf_readLineDF(vcffh)
}
\arguments{
  \item{vcffh}{VCF file handle}
}
\details{
Reads a line of data from the given VCF file, splits it up into its components (fields) and
fills a data.frame with the contents of the fields and names the entries according to the
header line of the VCF (e.g. CHROM, POS, ID, REF, ALT, ... ).
}
\value{
A data frame
}
\examples{
	vcffile <- vcf_open( system.file( "extdata", "ex.vcf.gz" , package="WhopGenome") )
	d <- vcf_readLineDF( vcffile )
}
\author{
Ulrich Wittelsbuerger
}
