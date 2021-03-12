\name{vcf_countSNPs}
\alias{vcf_countSNPs}
\alias{vcf_countBiallelicSNPs}
\title{
Count how many entries in the selected region 
}
\description{
Reads all data in the currently selected region of the given VCF file and counts how many loci with SNPs or biallelics SNPs respectively, are encountered.
}
\usage{
vcf_countSNPs( vcffh )
vcf_countBiallelicSNPs( vcffh )
}
\arguments{
\item{vcffh}{Handle to a VCF file, as returned by vcf_open}
}
\details{
For certain cases, like pre-allocating variables, it can be useful to know how many SNPs are present in a certain region.
In order to reduce the effort of this task and its impact on runtime to a minimum, the functions vcf_countSNPs and vcf_countBiallelicSNPs
were implemented.
Take note that they do not automatically 'restart' from the beginning of the selected region but continue from the current position. Use
vcf_restartregion to make sure that all SNPs in the currently set region are counted.
}
\value{
An integer number is returned: the number of SNPs or biallelic SNPs.
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
	vcf_countSNPs( vcffh )
}
\seealso{
vcf_restartregion
}
