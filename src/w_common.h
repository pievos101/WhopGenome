/*
**
**		WhopGen
**
**		WHOle genome population genetics with popGEN
**
**
**
**
**
**

	TODO:






**
*/

#ifndef		W_COMMON_H_INCLUDED
#define		W_COMMON_H_INCLUDED

//*
//*			INCLUDES
//*


	//	standard includes
	//
#include	<stdlib.h>				//string comparisons,...      
#include	<sys/time.h>			//gettimeofday

#define		STRICT_R_HEADERS		//to solve problems on Windows with Realloc, ERROR definitions made by R headers

	//	R includes
	//
#include	<R.h>
#include	<Rinternals.h>
#include	<R_ext/Rdynload.h>

	//
#include	<string>
#include	<vector>
#include	<strings.h>

//*
//*			DEFINES
//*


	// support some variations in VCF file format, e.g. GT being not first field, and "\" as genotype separator char besides "/" and "|" etc.
//#define			RELAXED_STANDARD		1

#define			Rboolean_TRUE	((Rboolean)1)

		//	--	debug mode related defines
		//
		//
#define			DEBUG			0

#if DEBUG
#	define		ONDBG			if( true )
#else
#	define		ONDBG			if( false )
#endif

#define			DBG				if( debug_level ) Rprintf

	//
	//
#define		EXPORT				extern "C"

	//
	//
#define		guard				try{
	
#define		unguards			}catch(const char*err){Rprintf("EXCEPTION CAUGHT:'%s' in file %s line %d\n",err,__FILE__,__LINE__);}catch(...){Rprintf("UNKNOWN EXCEPTION CAUGHT in file %s line %d\n",__FILE__,__LINE__);}


	//
	//


	//Nucleotide codes
	//		TCGAN-
	//
#define		NUC_T				1	/*also U*/
#define		NUC_C				2
#define		NUC_G				3
#define		NUC_A				4
#define		NUC_N				5
#define		GAP					6

#define		TAB					'\t'
#define		COMMA				','
#define		DOT					'.'
#define		COLON				':'
#define		SEMICOLON			';'
#define		SLASH				'/'
#define		BACKSLASH			'\\'
#define		BAR					'|'

	//
	//
#define		CHKSTR( nam, len )				CHKTYPE( nam, String, len )
#define		CHKTYPE( nam, type, len )		( is##type(nam) && (length(nam) == 1) )


#ifndef __FUNCTION_NAME__
    #ifdef WIN32   //WINDOWS
        #define __FUNCTION_NAME__   __FUNCTION__  
    #else          //*NIX
        #define __FUNCTION_NAME__   __func__ 
    #endif
#endif


//*
//*			STRUCTS
//*


	//	default SEXPs used to clear matrix column names and rows
//extern	SEXP	minus1_char, minus2_char;
#define			UNUSED_COLUMN_VALUE		mkChar("-1")
#define			UNUSED_CELL_VALUE		mkChar("-2")
#define			ERROR_CELL_VALUE		mkChar("-9999")


//*
//*			CLASSES
//*

class	vcff;




//*
//*			DATA
//*


extern	char	nucleotide_mapping[];


//*
//*			CODE
//*


	//
	//	internal functions related to filtering
	//

/*! Checks wether REF and ALT define a SNP (i.e. 1 single-nucleotide REF allele, up to 3 alt )
**
*/
inline bool	_internal_isSNP( const char * refptr, const char* altptr )
{
	//
	if( refptr[1] != '\t' )	//make sure its not a deletion
		return false;
		
	//	a SNP can have more than one alternate allele, but no multi-nucleotide alleles!
	//		=> test for letter followed by comma followed by letter [...] followed by \t
	
	//---------------------
	//examples:
	//
	//REF=A
	//
	//	and 
	//
	//ALT=C_		or
	//ALT=C,T_		or
	//ALT=C,T,G_
	//---------------------
	// ? can '.' or N alt-alleles happen? (i.e. deletion)
	
	int i=0;
	for( ; altptr[i] != 0 && altptr[i] != '\t'; i+=2 )
	{
		if( false == ((altptr[i] >= 'A' && altptr[i] <='Z') || (altptr[i]>='a' && altptr[i]<='z')) )
			break;
		if( altptr[i+1] == '\t' )
			return true;
		if( altptr[i+1] != ',' )
			break;
	}
	return false;
}




		//
		//	printing error message
		//

#define		RPRINT_ERROR		Rprintf("(!!) Error : %s : ", inf.caller_name ); Rprintf

									
#define		RPRINT_WARN			Rprintf("(?!) Warning : %s : ", inf.caller_name ); Rprintf

									
#define		RPRINT_VCFLINE		Rprintf("Current VCF line:\n[%s]\n", inf.f->getFieldPtr( 0 ) );

#define		DBG_RPRINT_WARN		if( debug_level ) Rprintf("(?!) Warning : %s : ", inf.caller_name ); if( debug_level )  Rprintf
#define		DBG_RPRINT_VCFLINE	if( debug_level ) Rprintf("Current VCF line:\n[%s]\n", inf.f->getFieldPtr( 0 ) );


//*
//*			EXTERNS
//*

extern	int		debug_level;

	//
	//
bool	filterLine( vcff* f );


	//	- utility functions -
	//
EXPORT	SEXP	WhopDebugLevel( void );
EXPORT	SEXP	SetWhopDebugLevel( SEXP lev );

	//	- debug log functions
	//
		void	df0( const char *fmt, ...);
#if DEBUG
#	define		df1		_df1
#	define		df2		_df2
#else
#	define		df1		if( false ) _df1
#	define		df2		if( false ) _df2
#endif
		void	_df1( const char *fmt, ...);
		void	_df2( const char *fmt, ...);

		//
		void	df(int level, const char *fmt, ...);

		//
		int		strcmp_cis( const char *a, const char *b );

	//	Open/Close
	//
EXPORT	SEXP	VCF_open( SEXP filename );
EXPORT	SEXP	VCF_close( SEXP vcfptr );
EXPORT	SEXP	VCF_eor( SEXP vcfptr );

	//
	//
EXPORT  SEXP VCF_buildindex( SEXP filename_sexp );

	//	Region
	//
EXPORT	SEXP	VCF_setRegion( SEXP vcfptr, SEXP tid, SEXP beg, SEXP end );
EXPORT	SEXP	VCF_getRegion( SEXP vcfptr );
EXPORT	SEXP	VCF_getCurrentRegionTid( SEXP vcfptr );
EXPORT	SEXP	VCF_getCurrentRegionBegin( SEXP vcfptr );
EXPORT	SEXP	VCF_getCurrentRegionEnd( SEXP vcfptr );
EXPORT	SEXP	VCF_restartRegion( SEXP vcfptr );

	//	Samples
	//
EXPORT	SEXP	VCF_getSampleNames( SEXP vcfptr );
EXPORT  SEXP    VCF_getSampleByName( SEXP vcfptr, SEXP str_samplename );
EXPORT	SEXP	VCF_getSelectedSamples( SEXP vcfptr );
EXPORT	SEXP	VCF_selectSamples( SEXP vcfptr, SEXP samplesvec );
EXPORT	SEXP	VCF_selectSampleByName( SEXP vcfptr, SEXP samplename );


	//	Reading
	//
EXPORT	SEXP	VCF_readLineRaw( SEXP vcfptr , SEXP str );
EXPORT	SEXP	VCF_readLineRawFiltered( SEXP vcfptr , SEXP str );
EXPORT	SEXP	VCF_readLineTSV( SEXP vcfptr );
EXPORT	SEXP	VCF_readLineTSVFiltered( SEXP vcfptr );
#if 0
	EXPORT	SEXP	VCF_readLineDF( SEXP vcfptr );
	EXPORT	SEXP	VCF_readLineDFFiltered( SEXP vcfptr );
#endif
EXPORT	SEXP	VCF_getBial( SEXP vcfptr, SEXP mat );

EXPORT	SEXP	VCF_parseNextLine( SEXP vcfptr );
EXPORT	SEXP	VCF_parseNextSNP( SEXP vcfptr );

		//
EXPORT	SEXP	VCF_getChrom( SEXP vcfptr );
EXPORT	SEXP	VCF_getPos( SEXP vcfptr );
EXPORT	SEXP	VCF_getID( SEXP vcfptr );
EXPORT	SEXP	VCF_getRef( SEXP vcfptr );
EXPORT	SEXP	VCF_getAlt( SEXP vcfptr );
EXPORT	SEXP	VCF_getQual( SEXP vcfptr );
EXPORT	SEXP	VCF_getFilter( SEXP vcfptr );
EXPORT	SEXP	VCF_getInfo( SEXP vcfptr );
EXPORT	SEXP	VCF_getInfoField( SEXP vcfptr, SEXP fieldnam );
EXPORT	SEXP	VCF_getFormat( SEXP vcfptr );
EXPORT	SEXP	VCF_getSample( SEXP vcfptr, SEXP stridx );

	//	Matrix readers
	//
EXPORT	SEXP	VCF_readBialMultilineCodeMatrix( SEXP vcfptr, SEXP mat );

EXPORT	SEXP	VCF_snpmat_diplo_bial_geno_filtered( SEXP vcfptr, SEXP mat );                                //              geno            filter          bial
EXPORT	SEXP	VCF_snpmat_diplo_anyal_geno_filtered( SEXP vcfptr, SEXP mat );                               //              geno            filter          anyal
EXPORT	SEXP	VCF_snpmat_diplo_bial_geno_unfiltered( SEXP vcfptr, SEXP mat );                              //              geno            NOFILT          bial
EXPORT	SEXP	VCF_snpmat_diplo_anyal_geno_unfiltered( SEXP vcfptr, SEXP mat );                             //              geno            NOFILT          anyal
EXPORT	SEXP	VCF_snpmat_diplo_bial_ishet_filtered( SEXP vcfptr, SEXP mat );                               //              isHet           filter          bial
EXPORT	SEXP	VCF_snpmat_diplo_anyal_ishet_filtered( SEXP vcfptr, SEXP mat );                              //              isHet           filter          anyal
EXPORT	SEXP	VCF_snpmat_diplo_bial_ishet_unfiltered( SEXP vcfptr, SEXP mat );                             //              isHet           NOFILT          bial
EXPORT	SEXP	VCF_snpmat_diplo_anyal_ishet_unfiltered( SEXP vcfptr, SEXP mat );                    //              isHet           NOFILT          anyal
EXPORT	SEXP	VCF_snpmat_diplo_bial_hasalt_filtered( SEXP vcfptr, SEXP mat );
EXPORT	SEXP	VCF_snpmat_diplo_bial_hasalt_unfiltered( SEXP vcfptr, SEXP mat );
EXPORT	SEXP	VCF_snpmat_diplo_anyal_hasalt_filtered( SEXP vcfptr, SEXP mat );
EXPORT	SEXP	VCF_snpmat_diplo_anyal_hasalt_unfiltered( SEXP vcfptr, SEXP mat );
EXPORT	SEXP	VCF_snpmat_diplo_bial_nucodes_filtered( SEXP vcfptr, SEXP mat );
EXPORT	SEXP	VCF_snpmat_diplo_bial_nucodes_unfiltered( SEXP vcfptr, SEXP mat );
EXPORT	SEXP	VCF_snpmat_diplo_anyal_nucodes_filtered( SEXP vcfptr, SEXP mat );
EXPORT	SEXP	VCF_snpmat_diplo_anyal_nucodes_unfiltered( SEXP vcfptr, SEXP mat );
EXPORT	SEXP	VCF_snpmat_anyplo_bial_nucodes_filtered( SEXP vcfptr, SEXP mat );
EXPORT	SEXP	VCF_snpmat_anyplo_bial_nucodes_unfiltered( SEXP vcfptr, SEXP mat );
EXPORT	SEXP	VCF_snpmat_anyplo_anyal_nucodes_filtered( SEXP vcfptr, SEXP mat );
EXPORT	SEXP	VCF_snpmat_anyplo_anyal_nucodes_unfiltered( SEXP vcfptr, SEXP mat );
		// obsolete ones
EXPORT	SEXP	VCF_readIntoNucleotideMatrix( SEXP vcfptr, SEXP mat );	//?
EXPORT	SEXP	VCF_readIntoCodeMatrix( SEXP vcfptr, SEXP mat );
EXPORT	SEXP	read_snp_diplo_bial_int_altpresence( SEXP vcfptr, SEXP mat );
EXPORT	SEXP	read_snp_diplo_bial_int_nuclcodes( SEXP vcfptr, SEXP mat );
EXPORT	SEXP	read_snp_diplo_bial_str_allelechars( SEXP vcfptr, SEXP mat );
EXPORT	SEXP	read_snp_diplo_bial_str_01( SEXP vcfptr, SEXP mat );
EXPORT	SEXP	read_snp_diplo_bial_str_nuclcodes( SEXP vcfptr, SEXP mat );
EXPORT	SEXP	read_snp_anyploid_multiallelic_int_nuclcodes( SEXP vcfptr, SEXP mat );

	//	Rules
	//
EXPORT	SEXP	VCF_setRuleAction( SEXP vcfptr, SEXP ruleidx_sexp, SEXP action );
EXPORT	SEXP	VCF_setRuleEnabled( SEXP vcfptr, SEXP ruleidx_sexp );
EXPORT	SEXP	VCF_setRuleDisabled( SEXP vcfptr, SEXP ruleidx_sexp );
EXPORT	SEXP	VCF_setRuleColumn( SEXP vcfptr, SEXP ruleidx_sexp, SEXP columnidx_sexp );
EXPORT	SEXP	VCF_setRuleRefValues( SEXP vcfptr, SEXP ruleidx_sexp, SEXP ref1, SEXP ref2 );
EXPORT	SEXP	VCF_setRuleField( SEXP vcfptr, SEXP ruleidx_sexp, SEXP field );
EXPORT	SEXP	VCF_setRuleCmpOp( SEXP vcfptr, SEXP ruleidx_sexp, SEXP cmpop );


	//	Misc. utility
	//

EXPORT	SEXP	VCF_countSNPs( SEXP vcfptr );
EXPORT	SEXP	VCF_countBiallelicSNPs( SEXP vcfptr );
EXPORT	SEXP	VCF_isSNP( SEXP vcfptr );
EXPORT	SEXP	VCF_isInDel( SEXP vcfptr );

	//	Filters
	//
EXPORT	SEXP	VCF_describeFilterConfig( SEXP vcfptr );
EXPORT	SEXP	VCF_addFilter( SEXP vcfptr, SEXP columnidx, SEXP fieldnam, SEXP cmptype, SEXP action, SEXP arg1, SEXP arg2 );
EXPORT	SEXP	VCF_clearFilters( SEXP vcfptr );

	//		Tabix for R
	//
EXPORT	SEXP	tabix_open( SEXP filename );
EXPORT	SEXP	tabix_close( SEXP tabix );
EXPORT	SEXP	tabix_reopen( SEXP tabix );
EXPORT	SEXP	tabix_setRegion( SEXP tabix, SEXP tid, SEXP begin, SEXP endpos );
EXPORT	SEXP	tabix_restartRegion( SEXP tabix );
EXPORT	SEXP	tabix_getRegion( SEXP tabix );
EXPORT	SEXP	tabix_readLine( SEXP tabix );
EXPORT  SEXP	tabix_build( SEXP filename_sexp, SEXP sc_sexp, SEXP bc_sexp, SEXP ec_sexp, SEXP meta_sexp, SEXP lineskip_sexp );

	//		FaIdx
	//
EXPORT	SEXP	FAI_query4( SEXP faiptr, SEXP seqn, SEXP beg, SEXP end );
EXPORT	SEXP	FAI_query2( SEXP faiptr, SEXP regionstr );
EXPORT	SEXP	FAI_build( SEXP filename );
EXPORT	SEXP	FAI_open( SEXP filename );
EXPORT	SEXP	FAI_close( SEXP faiptr );
EXPORT	SEXP	FAI_reopen( SEXP faiptr );
//EXPORT	SEXP	faidx_test( SEXP fname, SEXP reg );

	//
	//
EXPORT  SEXP BGZF_compress( SEXP in_filename_sexp, SEXP out_filename_sexp );

	//		GTF <nonfunctional>
	//
EXPORT	SEXP	gtftest( SEXP vcfptr );
	
	//
	//
#include	"w_tools.h"

	//
	//
#include	"w_rsupport.h"

	//
	//
#include	"w_tabix.h"

	//
	//
#include	"w_vcf.h"


#endif //W_COMMON_H_INCLUDED
