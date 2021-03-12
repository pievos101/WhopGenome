\name{vcf_rule.setrefvalues}
\alias{vcf_rule.setrefvalues}
\title{
Set reference values for a filtering rule's comparison operation.
}
\description{
Set the reference values 1 and 2 for the comparison operation of rule <ruleidx>.
Soem comparison operations need only the first <ref1> reference value and ignore <ref2>.
}
\usage{
vcf_rule.setrefvalues( vcffh, ruleidx, ref1, ref2 )
}
\arguments{
  \item{vcffh}{VCF file handle}
  \item{ruleidx}{name of column containing the to-be-checked values}
  \item{ref1}{name of the subfield or "" to check}
  \item{ref2}{Type of comparison to perform. See Details}
}
\details{
Certain VCF read functions support the fast pre-filtering mechanism of WhopGenome.
The pre-filtering mechanism is a list of rules that describe how and what to check in
SNP descriptions and is executed very quickly without using any R code.
Every rule specifies the column of data (e.g. INFO, POS, FILTER), the key in the column (e.g. AF in the INFO column),
the type of comparison , reference values to compare against and whether to keep or drop the line if the rule matches.
}
\value{
TRUE on success, FALSE if it failed.
}
\examples{
	##
	##	Example:
	##
	vcffile <- vcf_open( system.file( "extdata" , "ex.vcf.gz" , package="WhopGenome" ) )

	#
	#
	vcf_setregion(vcffile, "Y", 50000, 51000 )
	
	#
	# USELESS filter : # filter out SNPs with rule  "DROP if (0.0 < INFO:AF < 0.2)"
	#	pos 50001 has AF=0.285 , for which (0 < 0.285 < 0.2) is true 
	#
	vcf_addfilter( vcffile, "INFO", "AF", "FLT_CMP_OO", 0, 0.2, "DROP" )
	vcf_describefilters( vcffile )
	
	vcf_readLineVecFiltered( vcffile )	# pos 50001
	vcf_readLineVecFiltered( vcffile )	# pos 50002
	
	#
	#
	vcf_setregion(vcffile, "Y", 50000, 51000 )
	
	#CORRECT rule:
	# filter out SNP at pos 50001 with INFO:AF=0.285 with rule  "DROP if (0.2 < INFO:AF < 0.3)"
	#
	vcf_rule.setrefvalues( vcffile , 0 , 0.2, 0.3 )
	vcf_describefilters( vcffile )
	
	vcf_readLineVecFiltered( vcffile )	# pos 50002
	vcf_readLineVecFiltered( vcffile )	# pos 50003
}
\author{
Ulrich Wittelsbuerger
}
