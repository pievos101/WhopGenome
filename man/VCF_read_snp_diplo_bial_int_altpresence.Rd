\name{VCF_read_snp_diplo_bial_int_altpresence}
\alias{VCF_read_snp_diplo_bial_int_altpresence}
\alias{VCF_read_snp_diplo_bial_int_nuclcodes}
\alias{VCF_read_snp_diplo_bial_str_01}
\alias{VCF_read_snp_diplo_bial_str_allelechars}
\alias{VCF_read_snp_diplo_bial_str_nuclcodes}
\title{
(OBSOLETE) Read batch of biallelic SNP data into matrices
}
\description{
OBSOLETE : please refer to documentation for "VCF_snpmat_diplo_bial_geno_filtered" to
find out about their replacements.\cr


Reads biallelic SNP data in different representations into pre-allocated matrices.
}
\usage{
VCF_read_snp_diplo_bial_int_altpresence( vcffh , mat )
VCF_read_snp_diplo_bial_int_nuclcodes( vcffh , mat )
VCF_read_snp_diplo_bial_str_01( vcffh , mat )
VCF_read_snp_diplo_bial_str_allelechars( vcffh , mat )
VCF_read_snp_diplo_bial_str_nuclcodes( vcffh , mat )
}
\arguments{
\item{vcffh}{VCF file handle as returned by vcf_open}
\item{mat}{A matrix of either integer or string type, corresponding to _str_ or _int_ named methods}
}
\details{
OBSOLETE : please refer to documentation for "VCF_snpmat_diplo_bial_geno_filtered" to
find out about their replacements.\cr

Prerequisites are:
- a valid, open VCF file handle, passed as vcffh
- a valid sample selection (vcf_getsamples,vcf_getselectedsamples,vcf_selectsmaples)
- a properly set region (vcf_setregion)
- and a result matrix, mat.

The matrix will be filled with allele data in one of 4 encodings and needs to be of either integer or character data type, both depending on the called function (VCF_..._int_... or VCF_..._str_...) .
Each column corresponds to a SNP locus and each row to a sample.
The number of matrix columns determines the maximum number of SNP loci that are parsed from the VCF.
Column names are set to the position of the SNP, the row names are named after the samples they represent.
There must be at least as many rows as selected samples. Unused rows will be filled with default (N) data.
If there are not enough SNPs to fill all columns, the unused columns will be numbered with -1 and filled with N or -1.

VCF data is required to be diploid.

Representations:
\itemize{
	\item	int_altpresence		: 0 if genotype is REF/REF, 1 if not
	\item	int_nuclcodes		: integers, two-digit numbers: 11=TT, 12=TC, 13=TG, 14=TA, 15=TN, 21=CT, etc. (1=T, 2=C, 3=G, 4=A, 5=N)
	\item	str_01				: string, either 00, 01, 10 or 11 : 00=ref/ref, 11=alt/alt, 10=alt/ref, 01=ref/alt
	\item	str_allelechars		: string, nucleotides of both chromosomes (no indication of reference allele)
	\item	str_nuclcodes		: string, two-digit numbers: 11=TT, 12=TC, 13=TG, 14=TA, 15=TN, 21=CT, etc. (1=T, 2=C, 3=G, 4=A, 5=N)
}
}
\value{
TRUE or FALSE
}
\author{
Ulrich Wittelsbuerger
}
\examples{
	warning("These functions are obsolete! Consult VCF_snpmat_diplo_bial_geno_filtered etc.")
}
\seealso{
\code{\link{VCF_snpmat_diplo_bial_geno_filtered}}
}
