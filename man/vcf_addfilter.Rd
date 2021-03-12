\name{vcf_addfilter}
\alias{vcf_addfilter}
\title{
Add a condition for SNP filtering from VCF files.
}
\description{
Add a condition for filtering SNPs based on any column in a given VCF file.
}
\usage{
vcf_addfilter(vcf, columnnam, fieldnam, cmptype, cmpvalue1, cmpvalue2 = 0, action)
}
\arguments{
  \item{vcf}{VCF file handle}
  \item{columnnam}{name of column containing the to-be-checked values}
  \item{fieldnam}{name of the subfield or "" to check}
  \item{cmptype}{Type of comparison to perform. See Details}
  \item{cmpvalue1}{Comparison reference value 1 or lower bound}
  \item{cmpvalue2}{Comparison reference value 2 or upper bound}
  \item{action}{Action to take if comparison matches : NOP, SKIP, KEEP or fails: SKIP_NOT, KEEP_NOT}
}
\details{
Parameter 'columnnam':
	Name of a VCF column, in which the data of interest is stored
Parameter 'fieldnam':
	For the INFO and samples columns, the key under which the interesting data is stored.
	Example:
	vcf_addfilter( vcffile , "INFO", "H2", "DOES_EXIST", 0, 0, "DROP_NOT" )
	would cause any subsequent calls to read functions that perform filtering to drop lines that do not have the "H2" key in the INFO column,
	which indicates that the SNP is not marked as being registered in HapMap2. The parameters <ref1> and <ref2> are not used by the "DOES_EXIST"
	operation.

Comparison types:
	\itemize{
	\item DOES_EXIST		Rule matches, if in column named by <columnnam> is a key with the same name as in <fieldnam>
	}
	for integer values:
	\itemize{
		\item INT_CMP			is value = ref1 ?
		\item INT_CMP_OO		is value in open range (ref1, ref2) 
		\item INT_CMP_OC		is value in half-closed range (ref1, ref2]
		\item INT_CMP_CO		is value in half-closed range [ref1, ref2)
		\item INT_CMP_CC		is value in closed range [ref1, ref2]
	}
	for floating point values:
	\itemize{
		\item FLT_CMP			is value = ref1 ?
		\item FLT_CMP_OO		is value in open range (ref1, ref2) 
		\item FLT_CMP_OC		is value in half-closed range (ref1, ref2]
		\item FLT_CMP_CO		is value in half-closed range [ref1, ref2)
		\item FLT_CMP_CC		is value in closed range [ref1, ref2]
	}
}
\value{
Success status: TRUE on success, FALSE if the rule could not be added.
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
