\name{vcf_close}
\alias{vcf_close}
\title{
Close a VCF file previously opened with vcf_open.
}
\description{
Closes the VCF file described by the given handle and prevents subsequent use.
}
\usage{
vcf_close(vcf_filehandle)
}
\arguments{
  \item{vcf_filehandle}{A VCF filehandle returned by vcf_open}
}
\details{
Use .Call("VCF_close", vcf_filehandle ) to eliminate the overhead of using the R wrapper function.
}
\value{
None
}
\author{
Ulrich Wittelsbuerger
}
\examples{
	##
	##	Example:
	##
	vcffile <- system.file( "extdata" , "ex.vcf.gz" , package="WhopGenome" )
	vcffile
	vcffh <- vcf_open( vcffile )
	vcffh
	
	##
	##
	vcf_close( vcffh )
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
vcf_open
}
