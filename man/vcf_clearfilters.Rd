\name{vcf_clearfilters}
\alias{vcf_clearfilters}
\title{
Removes all filter steps.
}
\description{
Removes all active filters, no pre-filtering of returned lines will take place.
There is no function to undo this step.
}
\usage{
vcf_clearfilters(vcffh)
}
\arguments{
  \item{vcffh}{VCF file handle}
}
\details{
Use .Call("VCF_clearFilters", vcffh ) to eliminate the overhead of using the R wrapper function.
}
\value{
None.
}
\author{
Ulrich Wittelsbuerger
}
\examples{
	##
	##	Example:
	##
	vcffile <- vcf_open( system.file( "extdata" , "ex.vcf.gz" , package="WhopGenome" ) )
	vcf_setregion(vcffile, "Y", 1, 100000 )
	vcf_addfilter( vcffile, "POS", "", "INT_CMP_OO",
					as.integer(49005), as.integer(49007), "DROP" )
	vcf_describefilters( vcffile )
	####
	####
	vcf_readLineVecFiltered( vcffile )
	vcf_readLineVecFiltered( vcffile )
	vcf_readLineVecFiltered( vcffile )
	#######
	#######
	vcf_clearfilters( vcffile )
	vcf_describefilters( vcffile )
	vcf_restartregion( vcffile )
	####
	####
	vcf_readLineVecFiltered( vcffile )
	vcf_readLineVecFiltered( vcffile )
	vcf_readLineVecFiltered( vcffile )
	##
	##
	vcf_close( vcffile )
}
