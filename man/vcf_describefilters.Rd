\name{vcf_describefilters}
\alias{vcf_describefilters}
\title{
Prints description of current filter rules
}
\description{
Prints a better understandable description of the filter rules currently active for the given VCF file.
}
\usage{
vcf_describefilters(vcffh)
}
\arguments{
  \item{vcffh}{VCF file handle}
}
\details{
Use .Call("VCF_describeFilterConfig", filename ) to eliminate the overhead of using the R wrapper function.
Note the different naming of the library function!
}
\value{
None.
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
\author{
Ulrich Wittelsbuerger
}
