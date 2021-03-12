\name{vcf_rule.setaction}
\alias{vcf_rule.setaction}
\title{
Sets the kind of action to take when a rule matches (or does not match).
}
\description{
The value that rule number <ruleidx> should inspect is stored in the column named <column>, e.g. "INFO" or "POS".
}
\usage{
vcf_rule.setaction( vcffh, ruleidx, action )
}
\arguments{
  \item{vcffh}{VCF file handle}
  \item{ruleidx}{Filter rule to change}
  \item{action}{name of an action, see below}
}
\details{
Recognised values for 'action':

	NOP				do nothing
	SKIP			drop line on match, read next line
		DROP		
	KEEP			keep line on match, do not test further
	SKIP_NOT		drop line if not matching
		DROP_NOT
		SKIP_IF_NOT
		DROP_IF_NOT
	KEEP_NOT		keep line if not matching rule, do not test further

	Each action has also a 'disabled' variant, causing it to be ignored.

	NOP_DISABLED
	SKIP_DISABLED
		DROP_DISABLED
	KEEP_DISABLED
	SKIP_NOT_DISABLED
		DROP_NOT_DISABLED
	KEEP_NOT_DISABLED

The _NOT / _IF_NOT variants effectively invert the comparison operation.
(A == B) becomes (A != B), (1 <= A <= 100) becomes (A < 1 OR > 100).

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
