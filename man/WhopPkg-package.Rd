\name{WhopGenome-package}
\alias{WhopGenome-package}
\alias{WhopGenome}
\docType{package}
\title{
High-speed, high-specialisation population-scale whole-genome variation and sequence data access
}
\description{
WhopGenome provides read access to Variant Call Format files with maximum speed by means of C functions with
many specialised output formats and a configurable filtering engine. Allows indexing of FASTA files and any file
format using tab-separated columns, such as GFF, VCF and METAL, in preparation to high-speed access. Can read specified subsections of indexed FASTA files very fast.
It also provides many easy-to-use methods to access the UCSC Genome Browser SQL servers, the AmiGO gene ontology databases, PLINK .PED files and Bioconductor's organism
annotation databases.
}
\details{
\tabular{ll}{
Package: \tab WhopGenome\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-01-24\cr
License: \tab GPL-2\cr
}
	- Open a VCF file with
	handle <- vcf_open("filename")
	- Set a region of interest (chromosome/contig ID,start position, end position) with
	vcf_setregion(handle,"X",200000, 300000 )
	- Select (in this case the first 10) samples of interest:
	vcf_selectsamples( handle, vcf_getSampleNames(handle)[1:10] )
	- Read from the file via
	resvec <- vcf_readLineVec(handle)
}
\author{
Ulrich Wittelsbuerger \email{ulrich.wittelsbuerger@uni-duesseldorf.de}
}
\references{
	The 1000 Genomes Project \url{http://1000genomes.org/}

	The 1000 Genomes Project Consortium (2010), A map of human genome variation from population-scale sequencing.  Nature *467*, 1061-1073.

	Heng Li (2011), Tabix: Fast retrieval of sequence features from generic TAB-delimited files, Bioinformatics, doi: 10.1093/bioinformatics/btq671

	The Variant Call Format \url{http://www.1000genomes.org/wiki/Analysis/Variant\%20Call\%20Format/vcf-variant-call-format-version-41}

}

\keyword{ package }
\examples{
#vcfh <- .Call("VCF_open","/data/vcf/1000g/ALL.Chromosome1.consensus.vcf.gz",PACKAGE="WhopGenome")
}
