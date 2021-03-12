\name{vcf_getregion}
\alias{vcf_getregion}
\title{
Get description of currently selected chromosomal region.
}
\description{
Returns a textual description like 'chr3:10913000-20240100' representing the genomic range which is currently set.
}
\usage{
vcf_getregion( vcffh )
}
\arguments{
\item{vcffh}{Handle to the VCF for which the currently active region should be retrieved}
}
\details{
Returns the string describing the region which is currently set for the given VCF file.
The string has the form "<chromosome id>:<startpos>-<endpos>", e.g. "1:120300-130500",
where "1" is the identifier of the chromosome or contig stored in the file, 120300 is the
leftmost position in the sequence for which we want to get variation data and 130500 is the
rightmost position.
Because usually there is no variation data for every position, there is no guarantee that the
first reported SNP will be at position 120300.
Initially, before a region has been set by the user, the returned string is ":0-0".
}
\value{
NULL if vcffh is not a valid VCF filehandle as returned by vcf_open.
Otherwise, a region string.
}
\examples{
	vcffile <- vcf_open( system.file( "extdata", "ex.vcf.gz" , package="WhopGenome") )
	vcf_getregion( vcffile )
}
\author{
Ulrich Wittelsbuerger
}
