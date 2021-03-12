\name{vcf_selectsamples}
\alias{vcf_selectsamples}
\alias{vcf_getselectedsamples}
\alias{vcf_getsamples}
\title{
Set or query the active sample selection for a given VCF file or get the entire list of individuals.
}
\description{
Set (vcf_selectsamples) or query (vcf_getselectedsamples) which individuals are included in the returned results,
or get a list of selectable individuals.
}
\usage{
vcf_selectsamples( vcffh, sampleslist )
vcf_getselectedsamples( vcffh )
vcf_getsamples( vcffh )
}
\arguments{
\item{vcffh}{VCFhandle type as returned by vcf_open}
\item{sampleslist}{A vector containing the identifiers of the individuals}
}
\details{
When reading variants from VCF files, it is possible to restrict the returned results to a certain subset of the available
individuals (samples), e.g. members of a population or people with a certain trait.
With vcf_selectsamples the currently selected subset of individuals can be set for a given VCF file.
vcf_getselectedsamples returns the list of currently selected individuals and vcf_getsamples returns a list of all available identifiers in the file.

As with most other VCF functions, it is possible to call directly into the library to avoid some overhead.
Use .Call("VCF_getSampleNames", vcffh ) , .Call("VCF_getSelectedSamples", vcffh ) or .Call("VCF_selectSamples", vcffh, sampleslist ), respectively.
Note the different names!
}
\value{
A vector of strings representing the sample names selected or present in the VCF file.
}
\author{
Ulrich Wittelsbuerger
}
\examples{
	##
	##	Example:
	##
	vcffile <- vcf_open( system.file( "extdata" , "ex.vcf.gz" , package="WhopGenome" ) )
	allsamplenames <- vcf_getsamples( vcffile )
	vcf_selectsamples( vcffile , allsamplenames )
}
\seealso{
vcf_open
}
