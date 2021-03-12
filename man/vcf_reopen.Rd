\name{vcf_reopen}
\alias{vcf_reopen}
\title{
Reopen a closed or stale VCF file handle.
}
\description{
Allows re-opening a previously opened VCF file.
}
\usage{
vcf_reopen(vcffh)
}
\arguments{
  \item{vcffh}{VCF file handle as returned by vcf_open}
}
\details{
If a file handle was closed (vcf_close) or became stale (e.g. after an R crash), it can be reactivated with this function.
Use .Call("VCF_reopen", vcffh ) to eliminate the overhead of using the R wrapper function.
}
\value{
Returns the reopened file handle.
}
\author{
Ulrich Wittelsbuerger
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
vcf_open
}
\examples{
	vcffile <- vcf_open( system.file( "extdata", "ex.vcf.gz" , package="WhopGenome") )
	vcffile
	vcf_close( vcffile )
	vcffile
	vcf_reopen( vcffile )
	vcffile
}
