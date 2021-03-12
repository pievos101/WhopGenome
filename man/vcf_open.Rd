\name{vcf_open}
\alias{vcf_open}
\title{
Open the specified VCF file and return a filehandle for subsequent access.
}
\description{
Open the specified VCF file and return a filehandle for subsequent access.
}
\usage{
vcf_open(filename)
}
\arguments{
  \item{filename}{A filename of a tabix-indexed and gzip-compressed VCF file}
}
\details{
Use .Call("VCF_open", filename ) to eliminate the overhead of using the R wrapper function.
}
\value{
A VCF file handle, used in most VCF functions
}
\examples{
	vcffile <- vcf_open( system.file( "extdata", "ex.vcf.gz" , package="WhopGenome") )
}
\author{
Ulrich Wittelsbuerger
}
