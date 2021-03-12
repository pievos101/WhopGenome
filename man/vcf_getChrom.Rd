\name{vcf_getChrom}
\alias{vcf_getChrom}
\alias{vcf_getPos}
\alias{vcf_getID}
\alias{vcf_getRef}
\alias{vcf_getAlt}
\alias{vcf_getQual}
\alias{vcf_getFilter}
\alias{vcf_getInfo}
\alias{vcf_getInfoField}
\alias{vcf_getFormat}
\alias{vcf_getSample}
\title{
Return a specific piece of information from the last line processed with vcf_parseNextSNP or vcf_parsenextline.
}
\description{
Return a specific piece of information from the last line processed with vcf_parseNextSNP or vcf_parsenextline.
}
\usage{
	vcf_getChrom( vcffh )
	vcf_getPos( vcffh )
	vcf_getID( vcffh )
	vcf_getRef( vcffh )
	vcf_getAlt( vcffh )
	vcf_getQual( vcffh )
	vcf_getFilter( vcffh )
	vcf_getInfo( vcffh )
	vcf_getInfoField( vcffh, fieldnam )
	vcf_getFormat( vcffh )
	vcf_getSample( vcffh, stridx )
}
\arguments{
  \item{vcffh}{VCF file handle}
  \item{fieldnam}{Name of a key of the key-value-pairs stored in the INFO subfield}
  \item{stridx}{Name of a sample}
}
\details{
Use .Call("VCF_getChrom", filename ) to eliminate the overhead of using the R wrapper function.
Replace getChrom by getPos, getID, getRef, getAlt, getQual, getFilter, getInfo, getInfoField, getSample and
add the respective function arguments in the order given above to call the respective other function.
}
\value{
None if the call failed, otherwise the respective data from the last read line is extracted.
}
\author{
Ulrich Wittelsbuerger
}
\examples{
	##
	##	Example:
	##
	vcffile <- vcf_open( system.file( "extdata", "ex.vcf.gz" , package="WhopGenome") )
	vcf_parseNextSNP( vcffile )
	vcf_getChrom( vcffile )
	vcf_getPos( vcffile )
	vcf_getID( vcffile )
	vcf_getAlt( vcffile )
	vcf_getQual( vcffile )
	vcf_getFilter( vcffile )
	vcf_getInfoField( vcffile, "AA" )
}
\seealso{
\code{\link{vcf_isSNP}}
}

