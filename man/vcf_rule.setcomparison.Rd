\name{vcf_rule.setcomparison}
\alias{vcf_rule.setcomparison}
\title{
Set comparison operation for filtering rule.
}
\description{
For filtering rule <ruleidx> the comparison operation is set to <cmpop>,
which is one of the following strings:\cr
	string:		alternative:	meaning:\cr
	-------		------------	--------
	"HASKEY"	"DOES_EXIST"	key (specified as field) is present in column\cr
- integer comparisons:\cr
	"INT="		"INT_CMP"		ref1 = value\cr
	"INT()"		"INT_CMP_OO"	ref1  < value  < ref2\cr
	"INT(]"		"INT_CMP_OC"	ref1  < value <= ref2\cr
	"INT[)"		"INT_CMP_CO"	ref1 <= value  < ref2\cr
	"INT[]"		"INT_CMP_CC"	ref1 <= value <= ref2\cr
	
- floating point (real numbers):\cr
	"FLT=="		"FLT_CMP"		ref1 = value			FLT_CMP\cr
	"FLT()"		"FLT_CMP_OO"	ref1  < value  < ref2\cr
	"FLT(]"		"FLT_CMP_OC"	ref1  < value <= ref2\cr
	"FLT[)"		"FLT_CMP_CO"	ref1 <= value  < ref2\cr
	"FLT[]"		"FLT_CMP_CC"	ref1 <= value <= ref2\cr
}
\usage{
vcf_rule.setcomparison( vcffh, ruleidx, cmpop )
}
\arguments{
  \item{vcffh}{VCF file handle}
  \item{ruleidx}{number of rule in list}
  \item{cmpop}{One of the above strings, naming the comparison operation to perform}
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
	vcf_addfilter( vcffile, "POS", "", "INT_CMP_OO",
					as.integer(49005), as.integer(49007), "DROP" )
	vcf_describefilters( vcffile )
	
	vcf_rule.setcomparison( vcffile , 0, "INT_CMP_CC" )
	vcf_describefilters( vcffile )
}
\author{
Ulrich Wittelsbuerger
}
