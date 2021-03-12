\name{vcf_eor}
\alias{vcf_eor}
\title{
Determine whether all lines in the selected region have been read.
}
\description{
When reading SNP info within a region defined by VCF_setRegion, this function
returns TRUE/FALSE to indicate whether or not all lines within that region have been read.
}
\usage{
vcf_eor( vcffh )
}
\arguments{
  \item{vcffh}{Handle of a VCF file opened by VCF_open}
}
\details{
Use .Call("VCF_eor", vcffh ) to eliminate the overhead of using the R wrapper function.
}
\value{
TRUE if all SNPs inside the previously defined region have been read.
}
\examples{
	vcffile <- vcf_open( system.file( "extdata", "ex.vcf.gz" , package="WhopGenome") )
	vcf_setregion(vcffile, "Y", 1, 100000 )
	while( !vcf_eor(vcffile) )
	{
		vcf_readLineVec( vcffile )
	}
}
\author{
Ulrich Wittelsbuerger
}
