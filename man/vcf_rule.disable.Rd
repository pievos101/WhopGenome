\name{vcf_rule.disable}
\alias{vcf_rule.disable}
\alias{vcf_rule.enable}
\title{
Disable and enable processing of a rule
}
\description{
Filtering rules can be enabled and disabled.
Disabled rules are ignored by any filtering VCF read function.
}
\usage{
vcf_rule.disable( vcffh, ruleidx )
vcf_rule.enable( vcffh, ruleidx )
}
\arguments{
  \item{vcffh}{VCF file handle}
  \item{ruleidx}{number of rule}
}
\details{
Certain VCF read functions support the fast pre-filtering mechanism of WhopGenome.
The pre-filtering mechanism is a list of rules that describe how and what to check in
SNP descriptions and is executed very quickly without using any R code.
Every rule specifies the column of data (e.g. INFO, POS, FILTER), the key in the column (e.g. AF in the INFO column),
the type of comparison , reference values to compare against and whether to keep or drop the line if the rule matches.
}
\value{
TRUE if succeeded, FALSE if not
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
	vcf_readLineVecFiltered( vcffile )

	vcf_rule.disable( vcffile, 0 )
	vcf_describefilters( vcffile )
	vcf_restartregion( vcffile )
	vcf_readLineVecFiltered( vcffile )
}
\author{
Ulrich Wittelsbuerger
}
