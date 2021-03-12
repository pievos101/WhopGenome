\name{VCF_snpmat_diplo_bial_geno_filtered}
\alias{VCF_snpmat_diplo_bial_geno_filtered}
\alias{VCF_snpmat_diplo_bial_geno_filtered}
\alias{VCF_snpmat_diplo_anyal_geno_filtered}
\alias{VCF_snpmat_diplo_bial_geno_unfiltered}
\alias{VCF_snpmat_diplo_anyal_geno_unfiltered}
\alias{VCF_snpmat_diplo_bial_ishet_filtered}
\alias{VCF_snpmat_diplo_anyal_ishet_filtered}
\alias{VCF_snpmat_diplo_bial_ishet_unfiltered}
\alias{VCF_snpmat_diplo_anyal_ishet_unfiltered}
\alias{VCF_snpmat_diplo_bial_hasalt_filtered}
\alias{VCF_snpmat_diplo_bial_hasalt_unfiltered}
\alias{VCF_snpmat_diplo_anyal_hasalt_filtered}
\alias{VCF_snpmat_diplo_anyal_hasalt_unfiltered}
\alias{VCF_snpmat_diplo_bial_nucodes_filtered}
\alias{VCF_snpmat_diplo_bial_nucodes_unfiltered}
\alias{VCF_snpmat_diplo_anyal_nucodes_filtered}
\alias{VCF_snpmat_diplo_anyal_nucodes_unfiltered}
\alias{VCF_snpmat_anyplo_bial_nucodes_filtered}
\alias{VCF_snpmat_anyplo_bial_nucodes_unfiltered}
\alias{VCF_snpmat_anyplo_anyal_nucodes_filtered}
\alias{VCF_snpmat_anyplo_anyal_nucodes_unfiltered}
\alias{VCF_readIntoCodeMatrix}
\alias{read_snp_diplo_bial_int_altpresence}
\alias{read_snp_diplo_bial_int_nuclcodes}
\alias{read_snp_diplo_bial_str_allelechars}
\alias{read_snp_diplo_bial_str_01}
\alias{read_snp_diplo_bial_str_nuclcodes}
\title{
Read SNP matrices in one of various representations.
}
\description{
These functions read SNPs into matrices in a number of variations.
All VCF_snpmat functions read data into the provided integer matrices, except
for the <geno> format, which expects character/string-type matrices.
The functions return TRUE if the call was successfull, FALSE otherwise.

Each row corresponds to a sample, so make sure that the matrix you pass for <mat> has at least as many rows as selected samples.

Each column corresponds to a SNP. You can directly influence how many SNPs are read in at most by adjusting the number of columns of the matrices you pass.
These functions try to read as many SNPs as possible from the currently active region and fill unused columns with the value -2.

If the given matrices have dimnames, the column names are set to the genomic position (from the VCF_column "POS") of the SNPs.

}
\usage{
VCF_snpmat_diplo_bial_geno_filtered( vcffh, mat )
VCF_snpmat_diplo_anyal_geno_filtered( vcffh, mat )
VCF_snpmat_diplo_bial_geno_unfiltered( vcffh, mat )
VCF_snpmat_diplo_anyal_geno_unfiltered( vcffh, mat )

VCF_snpmat_diplo_bial_ishet_filtered( vcffh, mat )
VCF_snpmat_diplo_anyal_ishet_filtered( vcffh, mat )
VCF_snpmat_diplo_bial_ishet_unfiltered( vcffh, mat )
VCF_snpmat_diplo_anyal_ishet_unfiltered( vcffh, mat )

VCF_snpmat_diplo_bial_hasalt_filtered( vcffh, mat )
VCF_snpmat_diplo_bial_hasalt_unfiltered( vcffh, mat )
VCF_snpmat_diplo_anyal_hasalt_filtered( vcffh, mat )
VCF_snpmat_diplo_anyal_hasalt_unfiltered( vcffh, mat )

VCF_snpmat_diplo_bial_nucodes_filtered( vcffh, mat )
VCF_snpmat_diplo_bial_nucodes_unfiltered( vcffh, mat )
VCF_snpmat_diplo_anyal_nucodes_filtered( vcffh, mat )
VCF_snpmat_diplo_anyal_nucodes_unfiltered( vcffh, mat )
VCF_snpmat_anyplo_bial_nucodes_filtered( vcffh, mat )
VCF_snpmat_anyplo_bial_nucodes_unfiltered( vcffh, mat )
VCF_snpmat_anyplo_anyal_nucodes_filtered( vcffh, mat )
VCF_snpmat_anyplo_anyal_nucodes_unfiltered( vcffh, mat )


VCF_readIntoCodeMatrix( vcffh, mat )
read_snp_diplo_bial_int_altpresence( vcffh, mat )
read_snp_diplo_bial_int_nuclcodes( vcffh, mat )
read_snp_diplo_bial_str_allelechars( vcffh, mat )
read_snp_diplo_bial_str_01( vcffh, mat )
read_snp_diplo_bial_str_nuclcodes( vcffh, mat )
}
\arguments{
  \item{vcffh}{VCF file handle as returned by VCF_open}
  \item{mat}{Matrix to load data into}
}
\details{
The function names indicate what kind of data is read, how it is represented
and whether filtering rules are applied.
The names are constructed as follows:
	VCF_snpmat_[diploidy]_[allelicity]_[format]_[filtering]

For [diploidy] insert either
	diplo - SNPs from diploid data
	anyplo - SNPs of arbitrary ploidy.

For [allelicity] insert either
	bial - biallelic SNPs
	anyal - SNPs with an arbitrary number of alleles.


For [format] insert
	geno - genotype string ( typeof(mat) should be "character" ! )
	ishet - 1 or 0 depending on whether the genotype is heterozygous or not
	hasalt - 1 or 0 depending on whether the genotype features the alternate allele (either homo- or heterozygous).
	nucodes - nucleotide code, where ACTGN- are represented by a number between 1 and 6.

For [filtering] insert
	filtered - drop lines not matching filtering rules
	unfiltered - do not drop any lines

Example: the function
 VCF_snpmat_diplo_bial_nucodes_filtered
would read biallelic SNPs from diploid species data, turn their genotypes into numeric nucleotide codes and store them
in an integer matrix. Only SNPs that passed the currently active filtering rules

For [format] geno, provide a matrix of type "character".
For all other [format]s, provide a matrix of integer (not double!) type (typeof(mat) = "integer").



The following functions have become OBSOLETE:

	VCF_readIntoCodeMatrix
		- use VCF_snpmat_diplo_bial_nucodes_filtered() instead.

	read_snp_diplo_bial_int_altpresence
		- use VCF_snpmat_diplo_bial_hasalt_filtered() instead.

	read_snp_diplo_bial_int_nuclcodes
		- use VCF_snpmat_diplo_bial_nucodes_filtered() instead.

	read_snp_diplo_bial_str_allelechars( vcffh, mat )
		- use VCF_snpmat_diplo_bial_geno_filtered() instead.

	read_snp_diplo_bial_str_01( vcffh, mat )
		- use VCF_snpmat_diplo_bial_hasalt_filtered() with integer matrix.

	read_snp_diplo_bial_str_nuclcodes( vcffh, mat )
		- use VCF_snpmat_diplo_bial_nucodes_filtered() with integer matrix.

}
\value{
TRUE on success, FALSE if it failed.
}
\examples{
	##
	##	Example :
	##
	vcffile <- vcf_open( system.file( "extdata" , "ex.vcf.gz" , package="WhopGenome" ) )
	vcf_setregion(vcffile, "Y", 1, 100000 )
	
	sn <- vcf_getsamples( vcffile )
	vcf_selectsamples( vcffile , sn )
	
	m <- matrix( data=as.integer(0) , nrow = length(sn) , ncol = 4 )
	
	VCF_read_snp_diplo_bial_int_nuclcodes( vcffile , m )
	m

}
\author{
Ulrich Wittelsbuerger
}
\seealso{vcf_addfilter}
