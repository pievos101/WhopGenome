\name{vcf_rule.setcolumn}
\alias{vcf_rule.setcolumn}
\title{
Set column a rule should examine.
}
\description{
The value that rule number <ruleidx> should inspect is stored in the column named <column>, e.g. "INFO" or "POS".
}
\usage{
vcf_rule.setcolumn( vcffh, ruleidx, column )
}
\arguments{
  \item{vcffh}{VCF file handle}
  \item{ruleidx}{Filter rule to change}
  \item{column}{name of column containing the to-be-checked values}
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
	
	vcf_rule.setcolumn( vcffile , 0, "ID" )
	vcf_describefilters( vcffile )
}
\author{
Ulrich Wittelsbuerger
}
