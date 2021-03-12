\name{vcf_rule.setfield}
\alias{vcf_rule.setfield}
\title{
Set field or key of filtering rule.
}
\description{
Filtering rule number <ruleidx> should inspect the value stored under the key <field>.
This key is stored in the column defined for this rule (e.g. an INFO-column AF=0.34;RD=231;GQ=130 has keys AF,RD and GQ).
}
\usage{
vcf_rule.setfield( vcffh, ruleidx, field )
}
\arguments{
  \item{vcffh}{VCF file handle}
  \item{ruleidx}{number of rule in list}
  \item{field}{XXXX}
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
	# USELESS filter : # filter out SNPs with rule  "DROP if (0.0 < INFO:AA < 0.5)"
	#	AA= ancestral allele, is a floating point number!
	vcf_addfilter( vcffile, "INFO", "AA", "FLT_CMP_OO", 0, 0.5, "DROP" )
	vcf_describefilters( vcffile )
	
	vcf_readLineVecFiltered( vcffile )	# pos 50001
	vcf_readLineVecFiltered( vcffile )	# pos 50002
	
	#
	#
	vcf_setregion(vcffile, "Y", 50000, 51000 )
	
	#CORRECT rule:
	# filter out SNP at pos 50001 with INFO:AF=0.285 with rule  "DROP if (0.0 < INFO:AF < 0.5)"
	#
	vcf_rule.setfield( vcffile , 0 , "AF" )
	vcf_describefilters( vcffile )
	
	vcf_readLineVecFiltered( vcffile )	# pos 50002
	vcf_readLineVecFiltered( vcffile )	# pos 50003
	
}
\author{
Ulrich Wittelsbuerger
}
