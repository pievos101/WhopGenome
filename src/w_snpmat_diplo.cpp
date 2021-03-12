/*
**
**		WhopGenome
**
**		WHOle genome population genetics
**
**
**		
**		Read genotypes into matrices
**
**			- SNPs only
**
**
**
**

Variations of SNP-matrix functions:
===================================
				ploidy		allelicity	repre-	
				of			of			senta-	
				data:		sample-set:	tion:	
	---------------------------------------------------------------------------------------
	VCF_snpmat_	diplo_		bial_		nucode	_filtered
				anyplo_		anyal_		ishet
										hasalt
										nucstr
										
[nucleotide codes][heterozygous sample][nucleotides as string][sample has alt allele]

**
*/

//*
//*			INCLUDES
//*

#include	"w_common.h"
#include	"functor_common.h"



	//

using namespace std;


//*
//*			DEFINES
//*


//*
//*			STRUCTS
//*


//*
//*			CLASSES
//*



/*!
**
**
**
*/
class snpmat_read_nuccodes : public snpmat_read_info_int
{
public:

	snpmat_read_nuccodes( const char * cf ) : snpmat_read_info_int(cf){}

	int		nuccode_temp_lut[128];	//fast LUT to get nucleotide-codes from VCF-sample genotype-indices (i.e. like 0/1 , 1|1, 0|0 ...)
	char	nucleotides[128];		//DEBUG : the actual nucleotides making up REF and ALT
	int		numalleles;				//number of alleles (REF+ALT) in this line


	//! Read a sample's GT genotype information and construct a nucleotide code from it
	virtual	int		process_sample( const char* gt )
	{
//		ONDBG Rprintf("processLine [%s]\n",gt);
		
		int allelecode=0,
			allelenum
				;
		while( ! is_FORMAT_subfield_delimiter( *gt ) )
		{
//			ONDBG Rprintf("gt=%d [%c]\n",*gt,*gt);
			if( *gt == DOT )
			{
//				ONDBG Rprintf("DOT\n");
				allelecode = allelecode * 10 + NUC_N;
				gt++;
			}
			else if( is_ascii_digit( *gt ) )		//turn number into integer
			{
				
				//	turn number into integer
				//
//				ONDBG Rprintf("Parsing [%s]\n",gt );
				allelenum=0;
				while( is_ascii_digit( *gt ) )
				{
					allelenum = allelenum * 10 + ( *gt - '0' );
					gt++;
				}

				//
//				ONDBG Rprintf("	allelenum = %d (max %d)\n	nuccode = %d\n",allelenum,numalleles,nuccode_temp_lut[ allelenum ]);
				if( allelenum > numalleles )
				{
					Rprintf("	**** allelenum too large!\n");
					ptr[ cur_row ] = -1;	//WRITE ERROR CODE
					return 0;
				}
				
				//	turn integer into nucleotide code
				//
				allelenum = nuccode_temp_lut[ allelenum ];
				
				//	put (nucleotide code of allele on one chromosome) as new digit at the end of the (nucleotide code over all chromosomes for this sample)
				//
				allelecode = allelecode * 10 + ( allelenum );
				
				//
			}
			else if( ! is_GT_allele_divider( *gt ) )
			{
				Rprintf("ERROR : %c != div!\n", *gt );
				ptr[ cur_row ] = -1;	//WRITE ERROR CODE
				return 0;
			}
			else
				gt++;

		}//...
		
		//
		ptr[ cur_row ] = allelecode;
		
		// find out whether this line has variation in the genotypes of the selected samples
		//
		if( false == keep_column && cur_row > 0 )
		{
			keep_column = ( ptr[cur_row-1] != ptr[cur_row ] );
		}

		return 1;
	}///...process_GT


	//! Construct array of single-nucleotide codes for every allele in REF and ALT (e.g. REF=T ALT=A,C --> Array[] = {1,4,2})
	virtual	int		begin_line( void )
	{
		
		//
		numalleles=0;
		keep_column=!drop_uninformative_columns;	//if we want to drop uninformative columns, set keep_column to FALSE to make it check
		
		//
		const char*ref = f->getFieldPtr( REF );
		const char*alt = f->getFieldPtr( ALT );
		
		// look up nucleotide code for REF allele
		//
		nucleotides[0]=*ref;
		nuccode_temp_lut[ numalleles++ ] = nucleotide_mapping[ (unsigned char)*ref ];
		
		// look up nucleotide codes for all ALT alleles
		//
		for( int i=0; alt[i] != TAB ; i ++ )
		{
			if( alt[i] == COMMA )
				numalleles++;
			else
			{
				nuccode_temp_lut[ numalleles ] = nucleotide_mapping[ (unsigned char)alt[i] ];
				nucleotides[ numalleles ] = alt[i];
			}
		}
		
		//
		nucleotides[ numalleles+1 ]=0;

		return 1;

	}///...begin_line


	//!
	virtual	int		end_line( void )
	{
		return (keep_column?1:0);
	}///...end_line

};///...class snpmat_read_nuccodes



/*!
**
**
**
*/
class snpmat_read_nuccodes_str : public snpmat_read_info_str
{
public:

	snpmat_read_nuccodes_str( const char * cf ) : snpmat_read_info_str(cf){}

	int		nucleotidelut[128];	//fast LUT to get nucleotide-codes from VCF-sample genotype-indices (i.e. like 0/1 , 1|1, 0|0 ...)
	char	allelebuffer[128];		//
	int		numalleles;				//number of alleles (REF+ALT) in this line


	//! Read a sample's GT genotype information and construct a nucleotide code from it
	virtual	int		process_sample( const char* gt )
	{
//		ONDBG Rprintf("processLine [%s]\n",gt);
		
		int allelecode=0,
			bufpos=0
				;
		
		//	read genotype info
		//
		while( ! is_FORMAT_subfield_delimiter( *gt ) )
		{
			if( *gt == DOT )
			{
				allelebuffer[bufpos++] = 'N';
				gt++;
			}
			else if( is_ascii_digit( *gt ) )		//turn number into integer
			{
				
				//	turn number into integer
				//
				allelecode=0;
				while( is_ascii_digit( *gt ) )
				{
					allelecode = allelecode * 10 + ( *gt - '0' );
					gt++;
				}

				//
//				ONDBG Rprintf("	allelenum = %d (max %d)\n	nuccode = %d\n",allelenum,numalleles,nuccode_temp_lut[ allelenum ]);
				if( allelecode > numalleles )
				{
					Rprintf("	**** allelenum too large!\n");
					sptr[ cur_row ] = ERROR_CELL_VALUE;	//WRITE ERROR CODE
					return 0;
				}
				
				//	turn integer into nucleotide code
				//
				allelebuffer[bufpos++] = nucleotidelut[ allelecode ];
				
				//
			}
			else if( ! is_GT_allele_divider( *gt ) )
			{
				Rprintf("ERROR : %c != div!\n", *gt );
				sptr[ cur_row ] = ERROR_CELL_VALUE;	//WRITE ERROR CODE
				return 0;
			}
			else
				gt++;

		}//...while reading genotype
		
		//
		allelebuffer[bufpos]=0;
		sptr[ cur_row ] = mkChar( &allelebuffer[0] );
		
		// find out whether this line has variation in the genotypes of the selected samples
		//
		if( false == keep_column && cur_row > 0 )
		{
			keep_column = ( sptr[cur_row-1] != sptr[cur_row ] );
		}

		return 1;
	}///...process_GT


	//! Construct array of single-nucleotide codes for every allele in REF and ALT (e.g. REF=T ALT=A,C --> Array[] = {1,4,2})
	virtual	int		begin_line( void )
	{
		
		//
		numalleles=0;
		keep_column=!drop_uninformative_columns;	//if we want to drop uninformative columns, set keep_column to FALSE to make it check
		
		//
		const char*ref = f->getFieldPtr( REF );
		const char*alt = f->getFieldPtr( ALT );
		
		// look up nucleotide code for REF allele
		//
		nucleotidelut[ numalleles++ ] = *ref;
		
		// look up nucleotide codes for all ALT alleles
		//
		for( int i=0; alt[i] != TAB ; i ++ )
		{
			if( alt[i] == COMMA )
				numalleles++;
			else
			{
				nucleotidelut[ numalleles ] = alt[i];
			}
		}
		
		//
		nucleotidelut[ numalleles+1 ]=0;
		
		return 1;

	}///...begin_line


	//!
	virtual	int		end_line( void )
	{
		return (keep_column?1:0);
	}///...end_line

};///...class snpmat_read_nuccodes






/*!
**
**
**
*/
class snpmat_read_hasalt : public snpmat_read_info_int
{
public:
	snpmat_read_hasalt( const char * cf ) : snpmat_read_info_int(cf){}


	//! Write NULL into the matrix, if a sample genotype is homozygous reference, otherwise write a 1
	virtual	int		process_sample( const char* gt )
	{
		
		bool	hasalt = 0;

		while( ! is_FORMAT_subfield_delimiter( *gt ) )
		{
				
			//
			if( *gt >= '1' && *gt <= '9' )
			{
				hasalt=true;
				break;
			}

			gt++;

		}//...
		
		//
		ptr[ cur_row ] = ( hasalt ? 1 : 0 );
		
		// find out whether this line has variation in the genotypes of the selected samples
		//
		if( false == keep_column && cur_row > 0 )
		{
			keep_column = ( ptr[cur_row-1] != ptr[cur_row ] );
		}

		return 1;

	}///...process_GT


	//! Construct array of single-nucleotide codes for every allele in REF and ALT (e.g. REF=T ALT=A,C --> Array[] = {1,4,2})
	virtual	int		begin_line()
	{
		keep_column=!drop_uninformative_columns;	//if we want to drop uninformative columns, set keep_column to FALSE to make it check
		return 1;
	}///...begin_line


	//!
	virtual	int		end_line( )
	{
		//TODO find out whether to keep this column
		return (keep_column?1:0);
	}///...end_line

};///...class snpmat_read_hasalt






/*!
**
**
**
*/
class snpmat_read_ishet : public snpmat_read_info_int
{
public:
	snpmat_read_ishet( const char * cf ) : snpmat_read_info_int(cf){}


	//! Write NULL into the matrix, if a sample genotype is homozygous reference, otherwise write a 1
	virtual	int		process_sample( const char* gt )
	{
	//	ONDBG Rprintf("Parsing POS=%d [%s]\n",snppos,gt );
		int allelecode=-1,
			allelenum,
			res=0
				;
		while( ! is_FORMAT_subfield_delimiter( *gt ) )
		{
			if( *gt == DOT )
			{
				gt++;
			}
			else if( is_ascii_digit( *gt ) )
			{
				
				//	turn number into integer
				//
				allelenum=0;
				while( is_ascii_digit( *gt ) )
				{
					allelenum = allelenum * 10 + ( *gt - '0' );
					gt++;
				}

#ifdef VCF_FORMAT_CHECKS
				if( allelenum > numalleles )
				{
					Rprintf("	**** allelenum too large! incomplete matrix returned!\n");
					RPRINT_VCFLINE;
					ptr[ cur_row ] = -1;	//WRITE ERROR CODE
					return 0;
				}
#endif

				//
				//
			//	Rprintf("ishet: allelecode=%d allelenum=%d\n",allelecode,allelenum);
				if( allelecode == -1 )
				{
					allelecode=allelenum;
				}
				else if( allelecode!=allelenum )
				{
					res = 1;
					break;
				}

				//
			}
			else if( ! is_GT_allele_divider( *gt ) )
			{
				Rprintf("ERROR : expected allele divider but got '%c' !\n", *gt );
				ptr[ cur_row ] = -1;	//WRITE ERROR CODE
				return 0;
			}
			else
				gt++;

		}//...while not end of SAMPLE-GT field
		
		//
		ptr[ cur_row ] = res;

		// find out whether this line has variation in the genotypes of the selected samples
		//
		if( false == keep_column && cur_row > 0 )
		{
			//Rprintf("cmp : %d <> %d [r=%d]\n",ptr[cur_row-1] , ptr[cur_row ],cur_row);
			keep_column = ( ptr[cur_row-1] != ptr[cur_row ] );
		}

		return 1;

	}///...process_GT


	//! Construct array of single-nucleotide codes for every allele in REF and ALT (e.g. REF=T ALT=A,C --> Array[] = {1,4,2})
	virtual	int		begin_line( void )
	{
		keep_column=!drop_uninformative_columns;	//if we want to drop uninformative columns, set keep_column to FALSE to make it check
		return 1;
	}///...begin_line


	//!
	virtual	int		end_line( void )
	{
		return (keep_column?1:0);
	}///...end_line

};///...class snpmat_read_ishet


//*
//*			DATA
//*


//*
//*			EXTERNS
//*

SEXP _internal_VCF_snpmat_int( snpmat_read_info_int& inf, SEXP vcfptr, SEXP mat );

//*
//*			CODE
//*


							//
							//
							//
							//
							//		Genotype string
							//
							//
							//
							//

SEXP _internal_VCF_snpmat_str( snpmat_read_info_str& inf, SEXP vcfptr, SEXP mat );


/*!	Read matrix as heterozygous-YES/NO == 1/0
**
**
*/
EXPORT SEXP VCF_snpmat_diplo_bial_geno_filtered( SEXP vcfptr, SEXP mat )				//		geno		filter		bial
{
	snpmat_read_nuccodes_str		inf( __FUNCTION_NAME__ );
	inf.doFiltering = true;
	inf.biallelic_locus_required = true;
	inf.diploid_samples_only = true;
	return _internal_VCF_snpmat_str( inf , vcfptr, mat );
}

EXPORT SEXP VCF_snpmat_diplo_anyal_geno_filtered( SEXP vcfptr, SEXP mat )				//		geno		filter		anyal
{
	snpmat_read_nuccodes_str		inf( __FUNCTION_NAME__ );
	inf.doFiltering = true;
	inf.biallelic_locus_required = false;
	inf.diploid_samples_only = true;
	return _internal_VCF_snpmat_str( inf , vcfptr, mat );
}


EXPORT SEXP VCF_snpmat_diplo_bial_geno_unfiltered( SEXP vcfptr, SEXP mat )				//		geno		NOFILT		bial
{
	snpmat_read_nuccodes_str		inf( __FUNCTION_NAME__ );
	inf.doFiltering = false;
	inf.biallelic_locus_required = true;
	inf.diploid_samples_only = true;
	return _internal_VCF_snpmat_str( inf , vcfptr, mat );
}


EXPORT SEXP VCF_snpmat_diplo_anyal_geno_unfiltered( SEXP vcfptr, SEXP mat )				//		geno		NOFILT		anyal
{
	snpmat_read_nuccodes_str		inf( __FUNCTION_NAME__ );
	inf.doFiltering = false;
	inf.biallelic_locus_required = false;
	inf.diploid_samples_only = true;
	return _internal_VCF_snpmat_str( inf , vcfptr, mat );
}

							//
							//
							//
							//
							//		is Heterozygous
							//
							//
							//
							//





/*!	Read matrix as heterozygous-YES/NO == 1/0
**
**
*/
EXPORT SEXP VCF_snpmat_diplo_bial_ishet_filtered( SEXP vcfptr, SEXP mat )				//		isHet		filter		bial
{
	snpmat_read_ishet		inf( __FUNCTION_NAME__ );
	inf.doFiltering = true;
	inf.biallelic_locus_required = true;
	inf.diploid_samples_only = true;
	return _internal_VCF_snpmat_int( inf , vcfptr, mat );
}

EXPORT SEXP VCF_snpmat_diplo_anyal_ishet_filtered( SEXP vcfptr, SEXP mat )				//		isHet		filter		anyal
{
	snpmat_read_ishet		inf( __FUNCTION_NAME__ );
	inf.doFiltering = true;
	inf.biallelic_locus_required = false;
	inf.diploid_samples_only = true;
	return _internal_VCF_snpmat_int( inf , vcfptr, mat );
}


EXPORT SEXP VCF_snpmat_diplo_bial_ishet_unfiltered( SEXP vcfptr, SEXP mat )				//		isHet		NOFILT		bial
{
	snpmat_read_ishet		inf( __FUNCTION_NAME__ );
	inf.doFiltering = false;
	inf.biallelic_locus_required = true;
	inf.diploid_samples_only = true;
	return _internal_VCF_snpmat_int( inf , vcfptr, mat );
}


EXPORT SEXP VCF_snpmat_diplo_anyal_ishet_unfiltered( SEXP vcfptr, SEXP mat )			//		isHet		NOFILT		anyal
{
	snpmat_read_ishet		inf( __FUNCTION_NAME__ );
	inf.doFiltering = false;
	inf.biallelic_locus_required = false;
	inf.diploid_samples_only = true;
	return _internal_VCF_snpmat_int( inf , vcfptr, mat );
}










							//
							//
							//
							//
							//		has Alt Allele
							//
							//
							//
							//



/*!	Read matrix into "Sample Genotype has one of the Alternate Alleles" YES/NO == 1/0
**
**
*/
EXPORT SEXP VCF_snpmat_diplo_bial_hasalt_filtered( SEXP vcfptr, SEXP mat )
{
	snpmat_read_hasalt		inf( __FUNCTION_NAME__ );
	inf.doFiltering = true;
	inf.biallelic_locus_required = true;
	inf.diploid_samples_only = true;
	return _internal_VCF_snpmat_int( inf , vcfptr, mat );
}//...

EXPORT SEXP VCF_snpmat_diplo_bial_hasalt_unfiltered( SEXP vcfptr, SEXP mat )
{
	snpmat_read_hasalt		inf( __FUNCTION_NAME__ );
	inf.doFiltering = false;
	inf.biallelic_locus_required = true;
	inf.diploid_samples_only = true;
	return _internal_VCF_snpmat_int( inf , vcfptr, mat );
}//...

	//
	//		any-allelic
	//

EXPORT SEXP VCF_snpmat_diplo_anyal_hasalt_filtered( SEXP vcfptr, SEXP mat )
{
	snpmat_read_hasalt		inf( __FUNCTION_NAME__ );
	inf.doFiltering = true;
	inf.biallelic_locus_required = false;
	inf.diploid_samples_only = true;
	return _internal_VCF_snpmat_int( inf , vcfptr, mat );
}//...

EXPORT SEXP VCF_snpmat_diplo_anyal_hasalt_unfiltered( SEXP vcfptr, SEXP mat )
{
	snpmat_read_hasalt		inf( __FUNCTION_NAME__ );
	inf.doFiltering = false;
	inf.biallelic_locus_required = false;
	inf.diploid_samples_only = true;
	return _internal_VCF_snpmat_int( inf , vcfptr, mat );
}//...








							//
							//
							//
							//
							//		Nucleotide Codes
							//
							//
							//
							//



/*!	Read matrix of nucleotide codes for bi-allelic SNPs from diploid species
**
**
*/
EXPORT SEXP VCF_snpmat_diplo_bial_nucodes_filtered( SEXP vcfptr, SEXP mat )
{
	snpmat_read_nuccodes		inf( __FUNCTION_NAME__ );
	inf.doFiltering = true;
	inf.biallelic_locus_required = true;
	inf.diploid_samples_only = true;
	return _internal_VCF_snpmat_int( inf , vcfptr, mat );

}//...
EXPORT SEXP VCF_snpmat_diplo_bial_nucodes_unfiltered( SEXP vcfptr, SEXP mat )
{
	snpmat_read_nuccodes		inf( __FUNCTION_NAME__ );
	inf.doFiltering = false;
	inf.biallelic_locus_required = true;
	inf.diploid_samples_only = true;
	return _internal_VCF_snpmat_int( inf , vcfptr, mat );

}//...


EXPORT SEXP VCF_snpmat_diplo_anyal_nucodes_filtered( SEXP vcfptr, SEXP mat )
{
	snpmat_read_nuccodes		inf( __FUNCTION_NAME__ );
	inf.doFiltering = true;
	inf.biallelic_locus_required = false;
	inf.diploid_samples_only = true;
	return _internal_VCF_snpmat_int( inf , vcfptr, mat );

}//...
EXPORT SEXP VCF_snpmat_diplo_anyal_nucodes_unfiltered( SEXP vcfptr, SEXP mat )
{
	snpmat_read_nuccodes		inf( __FUNCTION_NAME__ );
	inf.doFiltering = false;
	inf.biallelic_locus_required = false;
	inf.diploid_samples_only = true;
	return _internal_VCF_snpmat_int( inf , vcfptr, mat );

}//...



/*
*/


EXPORT SEXP VCF_snpmat_anyplo_bial_nucodes_filtered( SEXP vcfptr, SEXP mat )
{
	snpmat_read_nuccodes		inf( __FUNCTION_NAME__ );
	inf.doFiltering = true;
	inf.biallelic_locus_required = true;
	inf.diploid_samples_only = false;
	return _internal_VCF_snpmat_int( inf , vcfptr, mat );

}//...
EXPORT SEXP VCF_snpmat_anyplo_bial_nucodes_unfiltered( SEXP vcfptr, SEXP mat )
{
	snpmat_read_nuccodes		inf( __FUNCTION_NAME__ );
	inf.doFiltering = false;
	inf.biallelic_locus_required = true;
	inf.diploid_samples_only = false;
	return _internal_VCF_snpmat_int( inf , vcfptr, mat );

}//...


EXPORT SEXP VCF_snpmat_anyplo_anyal_nucodes_filtered( SEXP vcfptr, SEXP mat )
{
	snpmat_read_nuccodes		inf( __FUNCTION_NAME__ );
	inf.doFiltering = true;
	inf.biallelic_locus_required = false;
	inf.diploid_samples_only = false;
	return _internal_VCF_snpmat_int( inf , vcfptr, mat );

}//...
EXPORT SEXP VCF_snpmat_anyplo_anyal_nucodes_unfiltered( SEXP vcfptr, SEXP mat )
{
	snpmat_read_nuccodes		inf( __FUNCTION_NAME__ );
	inf.doFiltering = false;
	inf.biallelic_locus_required = false;
	inf.diploid_samples_only = false;
	return _internal_VCF_snpmat_int( inf , vcfptr, mat );

}//...

		///
		///
		///
		///		Backwards-compatibility : previous function names
		///
		///
		///
		///
		///

		//
		//		integer matrix variants
		//


/*!	Reads a coding matrix for all biallelic SNPs of the given samples
**
**	- vcff : a VCFhandle
**
**
*/
EXPORT	SEXP	VCF_readIntoCodeMatrix( SEXP vcfptr, SEXP mat )
{
	//NOTE: used in PopGenome
	return VCF_snpmat_diplo_bial_nucodes_filtered( vcfptr, mat );
}


//!
EXPORT SEXP read_snp_diplo_bial_int_altpresence( SEXP vcfptr, SEXP mat )
{
	return VCF_snpmat_diplo_bial_hasalt_filtered( vcfptr, mat );
}

//!
EXPORT SEXP read_snp_diplo_bial_int_nuclcodes( SEXP vcfptr, SEXP mat )
{
	//NOTE: used in PopGenome
	return VCF_snpmat_diplo_bial_nucodes_filtered( vcfptr, mat );
}

		//
		//		string matrix variants
		//

//!
EXPORT SEXP read_snp_diplo_bial_str_allelechars( SEXP vcfptr, SEXP mat )
{
	return VCF_snpmat_diplo_bial_geno_filtered( vcfptr, mat );
}

//!
EXPORT SEXP read_snp_diplo_bial_str_01( SEXP vcfptr, SEXP mat )
{
	Rprintf("(!!) read_snp_diplo_bial_str_nuclcodes : redundant and obsolete - use VCF_snpmat_diplo_bial_hasalt_filtered with integer matrix instead!\n");
	return RBool::False();
}

//!
EXPORT SEXP read_snp_diplo_bial_str_nuclcodes( SEXP vcfptr, SEXP mat )
{
	Rprintf("(!!) read_snp_diplo_bial_str_nuclcodes : redundant and obsolete - use VCF_snpmat_diplo_bial_nucodes_filtered with integer matrix instead!\n");
	return RBool::False();
}

		//
		//
/*
//!
EXPORT SEXP read_snp_anyploid_multiallelic_int_nuclcodes( SEXP vcfptr, SEXP mat )
{
	TODO :
	- move code here
	- un-functor the design
}
*/
