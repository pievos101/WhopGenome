\name{vcf_buildindex}
\alias{vcf_buildindex}
\title{
Build Tabix-index required for processing VCF files.
}
\description{
Builds a Tabix-index for a VCF file that is already sorted and compressed.
}
\usage{
vcf_buildindex( filename )
}
\arguments{
\item{filename}{Name of VCF file}
}
\details{
Given the name of a VCF file, builds a Tabix-index file (automatically named <filename>.tbi) in the directory where the given VCF file is located.
Prerequisite is that the VCF file be sorted by chromosome and position as well as bgzip-compressed. Such files carry the extension .vcf.gz. Information
on how to sort data in VCF files can be found at < http://vcftools.sourceforge.net/docs.html >. Using bgzf_compress, you can thereafter compress the
file.
}
\value{
Returns TRUE if the index could be created or FALSE if not.
}
\author{
Ulrich Wittelsbuerger
}
\examples{
	##
	##	Example:
	##
}
\seealso{
\code{\link{tabix_build}}
\code{\link{bgzf_compress}}
}
