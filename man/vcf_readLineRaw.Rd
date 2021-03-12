\name{vcf_readLineRaw}
\alias{vcf_readLineRaw}
\alias{vcf_readLineRawFiltered}
\title{
Read a line of data from the given VCF file and return it as a string without postprocessing.
}
\description{
Read a line of data from the given VCF file and return it as a string without postprocessing.
}
\usage{
vcf_readLineRaw( vcffh )
%vcf_readLineRaw( vcffh, res )
vcf_readLineRawFiltered( vcffh )
%vcf_readLineRawFiltered( vcffh, res )
}
\arguments{
  \item{vcffh}{VCF file handle}
%  \item{res}{DEPRECATED : string variable receiving the data read from the file}
}
\details{
%NOTE : The 2-argument version is DEPRECATED and will cause warning messages!

vcf_readLineRawFiltered applies the filtering rules (see vcf_describefilters) and does not return any lines
that do not pass the filter rules.

Use .Call("VCF_readLineRaw", vcffh ) and .Call("VCF_readLineRawFiltered", vcffh ) respectively,
to eliminate the overhead of using the R wrapper function.

}
\value{

For the 1-argument versions:
	A raw string representing a line of data from the file or FALSE if no more lines to read

%For the 2-argument versions:
%	TRUE/FALSE whether another line could be read. The line itself is stored in the variable given for <res> 
}
\examples{
	vcffile <- vcf_open( system.file( "extdata", "ex.vcf.gz" , package="WhopGenome") )
	d <- vcf_readLineRaw( vcffile )
}
\author{
Ulrich Wittelsbuerger
}
