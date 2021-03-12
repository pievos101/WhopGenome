/*
**
**		WhopGen
**
**		WHOle genome population genetics with popGEN
**
**
**		main module - dll exports
**
**
**

	TODO:


**
*/

//*
//*			INCLUDES
//*

#include	"w_common.h"

	//
#include	<stdarg.h>			//for df1(), df0() debugging functions
#include	<ctype.h>			// for stricmp/strcasecmp implementation

	//	for internal benchmarks
	//
#include	<time.h>


	//
#include	<R_ext/Visibility.h>
#include	<R_ext/Boolean.h>

//*
//*			DEFINES
//*



#define		CHKTYPE( nam, type, len )		( is##type(nam) && (length(nam) == 1) )
#define		CHKSTR( nam, len )				CHKTYPE( nam, String, len )
#define		CHKINT( nam, len )				CHKTYPE( nam, Integer, len )


//*
//*			STRUCTS
//*


//*
//*			CLASSES
//*


#ifdef __linux__

/*!	
**
**
*/
class perftime
{
public:

	//
	perftime()
	{
		curtime.tv_sec =
		curtime.tv_nsec = 0;
	}
	
	//
	perftime( const perftime& P )
	{
		curtime = P.curtime;
	}
	
	//
	perftime(timespec inT )
	{
		curtime = inT;
	}
	
	//
	~perftime(){}
	
	//
	static timespec		diff(perftime pstart, perftime pend)
	{
		timespec start = pstart.curtime;
		timespec end = pend.curtime;
		timespec temp;
		
		if ((end.tv_nsec-start.tv_nsec)<0)
		{
			temp.tv_sec = end.tv_sec-start.tv_sec-1;
			temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
		}
		else
		{
			temp.tv_sec = end.tv_sec-start.tv_sec;
			temp.tv_nsec = end.tv_nsec-start.tv_nsec;
		}
		return temp;
	}//...diff

	//
	void	take( void )
	{
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &curtime);
	}
	
	//
	void	print( void )
	{
		perftime::print( curtime );
	}
	
	//
	static void	print( timespec t )
	{
		float alpha_nsec=( t.tv_nsec / 999999999.f );		// how much of a second has passed in nanoseconds ( range [0,1.0) )
		Rprintf("Elapsed time [%lds %ldns]\n",t.tv_sec,t.tv_nsec);
		Rprintf("Elapsed seconds (%f)\n", t.tv_sec + alpha_nsec);
	}
	
	//	
//protected:
	timespec	curtime;

};
#else
//#	error	Sorry, the baseline_{read|parse}.cpp programs are intended for linux only at the moment!
	//TODO issue compiler warning
#endif

//*
//*			DATA
//*

int			debug_level	= 0;


SEXP		vcfhandle_attrname_filename = R_NilValue;

//*
//*			EXTERNS
//*



//*
//*			CODE
//*

//!
int		strcmp_cis( const char *a, const char *b )
{
	unsigned char	ca,
					cb;
	do {
		ca = *a++;
		cb = *b++;
		ca = tolower(ca);
		cb = tolower(cb);
	}
	while( (ca == cb) && (ca) );
	return (int)(ca-cb);
}

//!
SEXP _internal_VcfGetAttrFilename( void )
{
	if( vcfhandle_attrname_filename == R_NilValue )
	{
		vcfhandle_attrname_filename = install("VCF.filename");
	}
	return vcfhandle_attrname_filename;
}

	//
	//		debug output control
	//



/*!
*/
EXPORT	SEXP	WhopDebugLevel( void )
{
	Rprintf("WhopGen : Current debug level is %d\n",debug_level);
	return R_NilValue;
}



/*!
*/
EXPORT	SEXP	SetWhopDebugLevel( SEXP lev )
{
	if( false == RNumeric::isInt( lev ) )
		return RBool::False();//R_NilValue;

	debug_level = RNumeric::getInt( lev );
	
	return RBool::True();//R_NilValue;
}




/*!
*/
void	df(int level, const char *fmt, ...)
{
	if( debug_level < level )
		return;
	char dfbuf[1024];
	va_list ap;
	va_start(ap, fmt);
	vsnprintf( dfbuf, sizeof(dfbuf),fmt, ap);
	va_end(ap);
	Rprintf("%s",&dfbuf[0] );
}




/*!
*/
void	df0( const char *fmt, ...)
{
	char dfbuf[1024];
	va_list ap;
	va_start(ap, fmt);
	vsnprintf( dfbuf, sizeof(dfbuf), fmt, ap);
	va_end(ap);
	Rprintf("%s",&dfbuf[0] );
}


#if DEBUG

/*!
*/
void	_df1( const char *fmt, ...)
{
	if( debug_level < 1 )
		return;
	char dfbuf[1024];
	va_list ap;
	va_start(ap, fmt);
	vsnprintf( dfbuf, sizeof(dfbuf), fmt, ap);
	va_end(ap);
	Rprintf("%s",&dfbuf[0] );
}




/*!
*/
void	_df2( const char *fmt, ...)
{
	if( debug_level < 2 )
		return;
	char dfbuf[1024];
	va_list ap;
	va_start(ap, fmt);
	vsnprintf( dfbuf, sizeof(dfbuf), fmt, ap);
	va_end(ap);
	Rprintf("%s",&dfbuf[0] );
}

#endif


									


		//
		//	finalizers
		//

/*!
**
**
**
*/
static void vcff_finalize(SEXP extPtr)
{
	df1("VCFF FINALIZE!\n");
	vcff* f = (vcff*)R_GetExtPtr( extPtr, "VCFhandle" );
	if( f )
	{
		df1("vcff_finalize : Finalizing VCFhandle!\n");
		VCF_close( extPtr );
		df1("vcff_finalize : Successfully finalized VCFhandle!\n");
	}
	else
	{
		df1("vcff_finalize : Could not finalize potential VCFhandle!\n");
	}

	//
}



		//
		//
		//




/*!	Open a Tabix-indexed VCF file and return a VCFhandle, a whopgen-managed EXTPTR SEXP
**
**
**
*/
EXPORT  SEXP VCF_open( SEXP filename )
{
	//
	if(! CHKSTR(filename,1) )
	{
		df0("VCF_open : filename is not a single string!");
		return R_NilValue;
	}
	
	//-----
	
	//open vcf with filename
	//	if failed: return NULL
	//
	vcff * f = new vcff( CHAR(STRING_ELT(filename, 0)) );
	if( f == 0 )
	{
		df0("VCF_open: Could not open file '%s' as tabix-indexed!\n", CHAR(STRING_ELT(filename, 0)));
		return R_NilValue;
	}
	
	//
	if( false == f->isValid() )
	{
		delete f;
		f=0;
		df0("VCF_open : Could not open file '%s' as tabix-indexed!\n",CHAR(STRING_ELT(filename, 0)));
		return R_NilValue;
	}
	
	//
	//
	//
#if 0
	int numseqs = f->getNumSequenceNames();
	for( int i=0; i < numseqs; i++ )
	{
		//
		const char * thisseqnam = f->getSequenceName(i);
		Rprintf("#%d=%s\n",i,thisseqnam);
		
		unsigned int	minr=0,maxr=2*1024*1024*1024;
		
		//
		for( unsigned int low = minr; low < maxr; low += 100000 )
		{
			const char * s = f->readNextLine();
			if( s )
			{
				while( *s != '\t' && *s != 0 )
					s++;
				int pos = atoi(s);
				
			}
		}
		
		//
	}
#endif

	//
	df1("(VCF_open) opened file '%s' is a VCF!\n",CHAR(STRING_ELT(filename, 0)));
	
	//
	//
	SEXP res;
	PROTECT(
		res = R_MakeExternalPtr( f, install("VCFhandle"), R_NilValue)
	);
	if( res == R_NilValue )
	{
		df0("VCF_open : could not create external pointer SEXP!\n");
		return res;
	}
	
	//
	R_RegisterCFinalizerEx(res, vcff_finalize, Rboolean_TRUE);
	
	//
	setAttrib( res, _internal_VcfGetAttrFilename(), filename );
	
	UNPROTECT( 1 );
	return res;

}


/*!	Close a VCFhandle, an EXTPTR SEXP of a whopgen-managed Tabix-indexed VCF file 
**
**
*/
EXPORT SEXP	VCF_close( SEXP vcfptr )
{
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		df0("VCF_close : parameter is not a VCFhandle or nil!\n");
		return RBool::False();
	}
	
	//
	R_ClearExternalPtr(vcfptr);
	
	//
	delete f;
	
	//
	return RBool::True();
}





/*!	Reopens a VCf file if the VCFhandle got stale
**
**
**
**
**
*/
EXPORT SEXP	VCF_reopen( SEXP vcfptr )
{
	//
	if( false == RType::IsExtPtr( vcfptr ) || ( strcasecmp( RExtPtr::getTag(vcfptr), "VCFhandle" ) != 0 ) )
	{
		df1("VCF_reopen : parameter is not an externalptr VCFhandle!\n");
		return RBool::False();
	}
	
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( f != 0 )
	{
		return RBool::True();
	}
	
	//
	SEXP filename = getAttrib( vcfptr, _internal_VcfGetAttrFilename() );
	
	//
	//open vcf with filename
	//	if failed: return False
	//
	f = new vcff( CHAR(STRING_ELT(filename, 0)) );
	if( f == 0 )
	{
		df0("VCF_reopen : Could not open file '%s' as tabix-indexed!\n",CHAR(STRING_ELT(filename, 0)));
		return RBool::False();
	}
	
	//	if loading failed, return with error
	//
	if( false == f->isValid() )
	{
		delete f;
		f=0;
		df0("VCF_reopen : Could not open file '%s' as tabix-indexed!\n",CHAR(STRING_ELT(filename, 0)));
		return RBool::False();
	}

	R_SetExternalPtrAddr( vcfptr, f );
	
	//
	return RBool::True();
}







/*!	Returns the next header-line, which was preparsed and stored in
**		 the VCFhandle object.
**	Returns NilValue, if last header-line was returned but also resets
**		the counter to begin at the first line again with next invocation.
**
**
*/
EXPORT SEXP	VCF_getHeaderLine( SEXP vcfptr, SEXP num )
{

	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		df0("VCF_getHeaderLine : Parameter 1 is not a VCFhandle!\n");
		return R_NilValue;
	}
	
	//
	if( false == isInteger( num ) )
	{
		df0("VCF_getHeaderLine : parameter 2 needs to be an integer!\n");
		return R_NilValue;
	}

	//
	int idx = INTEGER(num)[0];
	const char * s = f->getHeaderLine( idx );
	if( s != 0 )
	{
		RString res(s);
		return res.get();
	}

	//
	df1("No header line #%d to get!\n", idx );
	return R_NilValue;
}


	//
	//		Contig Identifiers / Chromosome Names (e.g. "1", "22", "chr2" )
	//
	//


/*!
**
**
**
*/
EXPORT SEXP	VCF_getContigNames( SEXP vcfptr )
{
	//
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("VCF_getContigNames : argument 1 is not a VCF!\n");
		return R_NilValue;
	}
	
	//
	int numseqs = f->getNumSequenceNames();
	//Rprintf("%d seqs\n",numseqs);
	
	//
	RString res;
	res.alloc( numseqs );
	for( int i=0; i < numseqs; i++ )
	{
		//Rprintf("#%d=%s\n",i,f->getSequenceName(i));
		res.set( f->getSequenceName(i), i );
	}
	
	//
	return res.get();
	
	//
}





/*!
**
**
**
*/
EXPORT SEXP	VCF_getNumContigs( SEXP vcfptr )
{
	//
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("VCF_getNumContigs : argument is not a VCF!\n");
		return R_NilValue;
	}
	
	//
	//
	RInteger	res( f->getNumSequenceNames() );
	return res.get();
}





		//
		//			VCF fields
		//
		//


/*!
**
**
**
*/
EXPORT  SEXP VCF_getFieldNames( SEXP vcfptr )
{
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("VCF_getFieldNames : argument 1 is not a VCF!\n");
		return R_NilValue;
	}
	
	//
	//
	RString st;
	st.alloc( f->getNumFields() );
	
	//
	for( unsigned int k=0; k < f->getNumFields(); k++ )
	{
		st.set( f->getFieldName( k ), k );
		//Rprintf("Token = '%s'\n",f->getFieldName( k ) );
	}
	return st.get();
	
	//
}





		//
		//			benchmarking tools
		//
		//


/*!
**
**
**
*/
EXPORT  SEXP int_benchmark_PARSE( SEXP vcfptr )
{
	//
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("VCF_int_benchmark_PARSE : argument 1 is not a VCF!\n");
		return R_NilValue;
	}

	//
	//
	int			numlines=0;

#ifdef __linux__
	perftime	tbegin,
				tend;
	tbegin.take();
#endif

	//
	while( f->parseNextLine() )
	{
		numlines++;
	}
#ifdef __linux__
	tend.take();
	perftime d = perftime::diff( tbegin, tend );
	perftime::print( d.curtime );
	Rprintf(" nsec and alpha = %ld and %f\n", d.curtime.tv_nsec, float( d.curtime.tv_nsec / 999999999.f ) );
	Rprintf(" lines/s = %f\n", float(numlines) / float(d.curtime.tv_sec + float( d.curtime.tv_nsec / 999999999.f )) );
#endif
	
	//
	Rprintf("num lines=%d\n",numlines);
	
	//
	return RBool::True();
	
	//
}


/*!
**
**
**
*/
EXPORT  SEXP int_benchmark_READ( SEXP vcfptr )
{
	//
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("int_benchmark_READ : argument 1 is not a VCF!\n");
		return R_NilValue;
	}

	//
	//
	int numlines=0;
	
	//
	//
	int current_line_len;
	const char * current_line;
#ifdef __linux__
	perftime	tbegin,tend;
	tbegin.take();
#endif
	while( (current_line = ti_read( f->getTabix(), f->getIter(), &current_line_len)) != 0 )
	{
		numlines++;
	}
#ifdef __linux__
	tend.take();
	perftime d = perftime::diff( tbegin, tend );
	perftime::print( d.curtime );
	Rprintf(" nsec and alpha = %ld and %f\n", d.curtime.tv_nsec, float( d.curtime.tv_nsec / 999999999.f ) );
	Rprintf(" lines/s = %f\n", float(numlines) / float(d.curtime.tv_sec + float( d.curtime.tv_nsec / 999999999.f )) );
#endif
	
	//
	Rprintf("num lines=%d\n",numlines);
	
	//
	return RBool::True();
	
	//
}


/*
 *
 * 
 * 
 * 
 * 
 * 
 */
R_CallMethodDef callMethods[]  = {
	//	Region
	//
	{"VCF_setRegion" , (DL_FUNC) & VCF_setRegion , 4 },		//VCF_setRegion( SEXP vcfptr, SEXP tid, SEXP beg, SEXP end )
	{"VCF_getCurrentRegionTid" , (DL_FUNC) & VCF_getCurrentRegionTid , 1 },		//VCF_getCurrentRegionTid( SEXP vcfptr )
	{"VCF_getCurrentRegionBegin" , (DL_FUNC) & VCF_getCurrentRegionBegin , 1 },		//VCF_getCurrentRegionBegin( SEXP vcfptr )
	{"VCF_getCurrentRegionEnd" , (DL_FUNC) & VCF_getCurrentRegionEnd , 1 },		//VCF_getCurrentRegionEnd( SEXP vcfptr )
	{"VCF_getRegion" , (DL_FUNC) & VCF_getRegion , 1 },		//VCF_getRegion( SEXP vcfptr )
	{"VCF_restartRegion" , (DL_FUNC) & VCF_restartRegion , 1 },		//VCF_restartRegion( SEXP vcfptr )
	{"VCF_eor" , (DL_FUNC) & VCF_eor , 1 },		//VCF_eor( SEXP vcfptr )
	
	//	Debugging tools
	//
	{"WhopDebugLevel" , (DL_FUNC) & WhopDebugLevel , 0 },		//WhopDebugLevel( void )
	{"SetWhopDebugLevel" , (DL_FUNC) & SetWhopDebugLevel , 1 },		//SetWhopDebugLevel( SEXP lev )
	
	//	VCF opening/closing/reopening
	//
	{"VCF_open" , (DL_FUNC) & VCF_open , 1 },		//VCF_open( SEXP filename )
	{"VCF_close" , (DL_FUNC) & VCF_close , 1 },		//VCF_close( SEXP vcfptr )
	{"VCF_reopen" , (DL_FUNC) & VCF_reopen , 1 },		//VCF_reopen( SEXP vcfptr )
	
	//	VCF misc
	//
	{"VCF_countBiallelicSNPs" , (DL_FUNC) & VCF_countBiallelicSNPs , 1 },		//VCF_countBiallelicSNPs( SEXP vcfptr )
	{"VCF_countSNPs" , (DL_FUNC) & VCF_countSNPs , 1 },						//VCF_countSNPs( SEXP vcfptr )
	{"VCF_buildindex" , (DL_FUNC) & VCF_buildindex , 1 },						//VCF_buildindex( SEXP filename_sexp )

	//	VCF information
	//
	{"VCF_getHeaderLine" , (DL_FUNC) & VCF_getHeaderLine , 2 },		//VCF_getHeaderLine( SEXP vcfptr, SEXP num )
	{"VCF_getContigNames" , (DL_FUNC) & VCF_getContigNames , 1 },		//VCF_getContigNames( SEXP vcfptr )
	{"VCF_getNumContigs" , (DL_FUNC) & VCF_getNumContigs , 1 },		//VCF_getNumContigs( SEXP vcfptr )
	{"VCF_getFieldNames" , (DL_FUNC) & VCF_getFieldNames , 1 },		//VCF_getFieldNames( SEXP vcfptr )
	
	//	Benchmarking helpers
	//
	{"int_benchmark_PARSE" , (DL_FUNC) & int_benchmark_PARSE , 1 },		//int_benchmark_PARSE( SEXP vcfptr )
	{"int_benchmark_READ" , (DL_FUNC) & int_benchmark_READ , 1 },		//int_benchmark_READ( SEXP vcfptr )
	
	//	FaIdx
	//
	{"FAI_query4" , (DL_FUNC) & FAI_query4 , 4 },		//FAI_query5( SEXP faiptr, SEXP seqn, SEXP beg, SEXP end )
	{"FAI_query2" , (DL_FUNC) & FAI_query2 , 2 },		//FAI_query3( SEXP faiptr, SEXP regionstr )
	{"FAI_build" , (DL_FUNC) & FAI_build , 1 },		//FAI_build( SEXP filename )
	{"FAI_open" , (DL_FUNC) & FAI_open , 1 },		//FAI_open( SEXP filename )
	{"FAI_close" , (DL_FUNC) & FAI_close , 1 },		//FAI_close( SEXP faiptr )
	{"FAI_reopen" , (DL_FUNC) & FAI_reopen , 1 },		//FAI_reopen( SEXP faiptr )
//	{"faidx_test" , (DL_FUNC) & faidx_test , 2 },		//faidx_test( SEXP fname, SEXP reg )
	
	//	Samples
	//
	{"VCF_selectSamples" , (DL_FUNC) & VCF_selectSamples , 2 },		//VCF_selectSamples( SEXP vcfptr, SEXP samplesvec )
	{"VCF_selectSampleByName" , (DL_FUNC) & VCF_selectSampleByName , 2 },		//VCF_selectSampleByName( SEXP vcfptr, SEXP samplename )
	{"VCF_getSelectedSamples" , (DL_FUNC) & VCF_getSelectedSamples , 1 },		//VCF_getSelectedSamples( SEXP vcfptr )
	{"VCF_getSampleNames" , (DL_FUNC) & VCF_getSampleNames , 1 },		//VCF_getSampleNames( SEXP vcfptr )
	
	//	Tests
	//
	{"VCF_isSNP" , (DL_FUNC) & VCF_isSNP , 1 },		//VCF_isSNP( SEXP vcfptr )
	{"VCF_isInDel" , (DL_FUNC) & VCF_isInDel , 1 },		//VCF_isInDel( SEXP vcfptr )
	
	//	Rules
	//
	{"VCF_clearFilters" , (DL_FUNC) & VCF_clearFilters , 1 },		//VCF_clearFilters( SEXP vcfptr )
	{"VCF_addFilter" , (DL_FUNC) & VCF_addFilter , 7 },		//VCF_addFilter( SEXP vcfptr, SEXP columnidx_sexp, SEXP fieldnam, SEXP cmptype, SEXP action, SEXP arg1, SEXP arg2 )
	{"VCF_describeFilterConfig" , (DL_FUNC) & VCF_describeFilterConfig , 1 },		//VCF_describeFilterConfig( SEXP vcfptr )
	{"VCF_setRuleAction" , (DL_FUNC) & VCF_setRuleAction , 3 },		//VCF_setRuleAction( SEXP vcfptr, SEXP ruleidx_sexp, SEXP action )
	{"VCF_setRuleEnabled" , (DL_FUNC) & VCF_setRuleEnabled , 2 },		//VCF_setRuleEnabled( SEXP vcfptr, SEXP ruleidx_sexp )
	{"VCF_setRuleDisabled" , (DL_FUNC) & VCF_setRuleDisabled , 2 },		//VCF_setRuleDisabled( SEXP vcfptr, SEXP ruleidx_sexp )
	{"VCF_setRuleColumn" , (DL_FUNC) & VCF_setRuleColumn , 3 },		//VCF_setRuleColumn( SEXP vcfptr, SEXP ruleidx_sexp, SEXP columnidx_sexp )
	{"VCF_setRuleRefValues" , (DL_FUNC) & VCF_setRuleRefValues , 4 },		//VCF_setRuleRefValues( SEXP vcfptr, SEXP ruleidx_sexp, SEXP ref1, SEXP ref2 )
	{"VCF_setRuleField" , (DL_FUNC) & VCF_setRuleField , 3 },		//VCF_setRuleField( SEXP vcfptr, SEXP ruleidx_sexp, SEXP field )
	{"VCF_setRuleCmpOp" , (DL_FUNC) & VCF_setRuleCmpOp , 3 },		//VCF_setRuleCmpOp( SEXP vcfptr, SEXP ruleidx_sexp, SEXP cmpop )
	
	//	Tabix generic
	//
	{"tabix_open" , (DL_FUNC) & tabix_open , 1 },		//tabix_open( SEXP filename )
	{"tabix_close" , (DL_FUNC) & tabix_close , 1 },		//tabix_close( SEXP tabix )
	{"tabix_reopen" , (DL_FUNC) & tabix_reopen , 1 },		//tabix_reopen( SEXP tabix )
	{"tabix_setRegion" , (DL_FUNC) & tabix_setRegion , 4 },		//tabix_setRegion( SEXP tabix, SEXP tid, SEXP begin, SEXP endpos )
	{"tabix_getRegion" , (DL_FUNC) & tabix_getRegion , 1 },		//tabix_getRegion( SEXP tabix )
	{"tabix_restartRegion" , (DL_FUNC) & tabix_restartRegion , 1 },		//tabix_restartRegion( SEXP tabix )
	{"tabix_readLine" , (DL_FUNC) & tabix_readLine , 1 },		//tabix_readLine( SEXP tabix )
	{"tabix_build" , (DL_FUNC) & tabix_build , 6 },		//tabix_build( SEXP filename_sexp, SEXP sc_sexp, SEXP bc_sexp, SEXP ec_sexp, SEXP meta_sexp, SEXP lineskip_sexp )
	
	//	Matrix readers
	//
	{"VCF_snpmat_diplo_bial_geno_filtered" , (DL_FUNC) & VCF_snpmat_diplo_bial_geno_filtered , 2 },				//VCF_snpmat_diplo_bial_geno_filtered( SEXP vcfptr, SEXP mat )                                //              geno            filter          bial
	{"VCF_snpmat_diplo_anyal_geno_filtered" , (DL_FUNC) & VCF_snpmat_diplo_anyal_geno_filtered , 2 },				//VCF_snpmat_diplo_anyal_geno_filtered( SEXP vcfptr, SEXP mat )                               //              geno            filter          anyal
	{"VCF_snpmat_diplo_bial_geno_unfiltered" , (DL_FUNC) & VCF_snpmat_diplo_bial_geno_unfiltered , 2 },			//VCF_snpmat_diplo_bial_geno_unfiltered( SEXP vcfptr, SEXP mat )                              //              geno            NOFILT          bial
	{"VCF_snpmat_diplo_anyal_geno_unfiltered" , (DL_FUNC) & VCF_snpmat_diplo_anyal_geno_unfiltered , 2 },			//VCF_snpmat_diplo_anyal_geno_unfiltered( SEXP vcfptr, SEXP mat )                             //              geno            NOFILT          anyal
	{"VCF_snpmat_diplo_bial_ishet_filtered" , (DL_FUNC) & VCF_snpmat_diplo_bial_ishet_filtered , 2 },				//VCF_snpmat_diplo_bial_ishet_filtered( SEXP vcfptr, SEXP mat )                               //              isHet           filter          bial
	{"VCF_snpmat_diplo_anyal_ishet_filtered" , (DL_FUNC) & VCF_snpmat_diplo_anyal_ishet_filtered , 2 },			//VCF_snpmat_diplo_anyal_ishet_filtered( SEXP vcfptr, SEXP mat )                              //              isHet           filter          anyal
	{"VCF_snpmat_diplo_bial_ishet_unfiltered" , (DL_FUNC) & VCF_snpmat_diplo_bial_ishet_unfiltered , 2 },			//VCF_snpmat_diplo_bial_ishet_unfiltered( SEXP vcfptr, SEXP mat )                             //              isHet           NOFILT          bial
	{"VCF_snpmat_diplo_anyal_ishet_unfiltered" , (DL_FUNC) & VCF_snpmat_diplo_anyal_ishet_unfiltered , 2 },		//VCF_snpmat_diplo_anyal_ishet_unfiltered( SEXP vcfptr, SEXP mat )                    //              isHet           NOFILT          anyal
	{"VCF_snpmat_diplo_bial_hasalt_filtered" , (DL_FUNC) & VCF_snpmat_diplo_bial_hasalt_filtered , 2 },			//VCF_snpmat_diplo_bial_hasalt_filtered( SEXP vcfptr, SEXP mat )
	{"VCF_snpmat_diplo_bial_hasalt_unfiltered" , (DL_FUNC) & VCF_snpmat_diplo_bial_hasalt_unfiltered , 2 },		//VCF_snpmat_diplo_bial_hasalt_unfiltered( SEXP vcfptr, SEXP mat )
	{"VCF_snpmat_diplo_anyal_hasalt_filtered" , (DL_FUNC) & VCF_snpmat_diplo_anyal_hasalt_filtered , 2 },			//VCF_snpmat_diplo_anyal_hasalt_filtered( SEXP vcfptr, SEXP mat )
	{"VCF_snpmat_diplo_anyal_hasalt_unfiltered" , (DL_FUNC) & VCF_snpmat_diplo_anyal_hasalt_unfiltered , 2 },		//VCF_snpmat_diplo_anyal_hasalt_unfiltered( SEXP vcfptr, SEXP mat )
	{"VCF_snpmat_diplo_bial_nucodes_filtered" , (DL_FUNC) & VCF_snpmat_diplo_bial_nucodes_filtered , 2 },			//VCF_snpmat_diplo_bial_nucodes_filtered( SEXP vcfptr, SEXP mat )
	{"VCF_snpmat_diplo_bial_nucodes_unfiltered" , (DL_FUNC) & VCF_snpmat_diplo_bial_nucodes_unfiltered , 2 },		//VCF_snpmat_diplo_bial_nucodes_unfiltered( SEXP vcfptr, SEXP mat )
	{"VCF_snpmat_diplo_anyal_nucodes_filtered" , (DL_FUNC) & VCF_snpmat_diplo_anyal_nucodes_filtered , 2 },		//VCF_snpmat_diplo_anyal_nucodes_filtered( SEXP vcfptr, SEXP mat )
	{"VCF_snpmat_diplo_anyal_nucodes_unfiltered" , (DL_FUNC) & VCF_snpmat_diplo_anyal_nucodes_unfiltered , 2 },	//VCF_snpmat_diplo_anyal_nucodes_unfiltered( SEXP vcfptr, SEXP mat )
	{"VCF_snpmat_anyplo_bial_nucodes_filtered" , (DL_FUNC) & VCF_snpmat_anyplo_bial_nucodes_filtered , 2 },		//VCF_snpmat_anyplo_bial_nucodes_filtered( SEXP vcfptr, SEXP mat )
	{"VCF_snpmat_anyplo_bial_nucodes_unfiltered" , (DL_FUNC) & VCF_snpmat_anyplo_bial_nucodes_unfiltered , 2 },	//VCF_snpmat_anyplo_bial_nucodes_unfiltered( SEXP vcfptr, SEXP mat )
	{"VCF_snpmat_anyplo_anyal_nucodes_filtered" , (DL_FUNC) & VCF_snpmat_anyplo_anyal_nucodes_filtered , 2 },		//VCF_snpmat_anyplo_anyal_nucodes_filtered( SEXP vcfptr, SEXP mat )
	{"VCF_snpmat_anyplo_anyal_nucodes_unfiltered" , (DL_FUNC) & VCF_snpmat_anyplo_anyal_nucodes_unfiltered , 2 },	//VCF_snpmat_anyplo_anyal_nucodes_unfiltered( SEXP vcfptr, SEXP mat )
	
		// obsolete
	{"VCF_readIntoCodeMatrix" , (DL_FUNC) & VCF_readIntoCodeMatrix , 2 },									//VCF_readIntoCodeMatrix( SEXP vcfptr, SEXP mat )
	{"read_snp_diplo_bial_int_altpresence" , (DL_FUNC) & read_snp_diplo_bial_int_altpresence , 2 },		//read_snp_diplo_bial_int_altpresence( SEXP vcfptr, SEXP mat )
	{"read_snp_diplo_bial_int_nuclcodes" , (DL_FUNC) & read_snp_diplo_bial_int_nuclcodes , 2 },			//read_snp_diplo_bial_int_nuclcodes( SEXP vcfptr, SEXP mat )
	{"read_snp_diplo_bial_str_allelechars" , (DL_FUNC) & read_snp_diplo_bial_str_allelechars , 2 },		//read_snp_diplo_bial_str_allelechars( SEXP vcfptr, SEXP mat )
	{"read_snp_diplo_bial_str_01" , (DL_FUNC) & read_snp_diplo_bial_str_01 , 2 },						//read_snp_diplo_bial_str_01( SEXP vcfptr, SEXP mat )
	{"read_snp_diplo_bial_str_nuclcodes" , (DL_FUNC) & read_snp_diplo_bial_str_nuclcodes , 2 },			//read_snp_diplo_bial_str_nuclcodes( SEXP vcfptr, SEXP mat )
	
		// ?
	{"read_snp_anyploid_multiallelic_int_nuclcodes" , (DL_FUNC) & read_snp_anyploid_multiallelic_int_nuclcodes , 2 },		//read_snp_anyploid_multiallelic_int_nuclcodes( SEXP vcfptr, SEXP mat )
	//{"read_snp_anyploid_multiallelic_int_nuclcodes" , (DL_FUNC) & read_snp_anyploid_multiallelic_int_nuclcodes , 2 },		//read_snp_anyploid_multiallelic_int_nuclcodes( SEXP vcfptr, SEXP mat )
	{"VCF_readBialMultilineCodeMatrix" , (DL_FUNC) & VCF_readBialMultilineCodeMatrix , 2 },		//VCF_readBialMultilineCodeMatrix( SEXP vcfptr, SEXP mat )
	
//	{"helper_read_intmatrix_functored" , (DL_FUNC) & helper_read_intmatrix_functored , XXX },		//helper_read_intmatrix_functored( SEXP vcfptr, MatrixLoaderBaseClass &fn, SEXP mat )
	{"VCF_readIntoNucleotideMatrix" , (DL_FUNC) & VCF_readIntoNucleotideMatrix , 2 },		//VCF_readIntoNucleotideMatrix( SEXP vcfptr, SEXP mat )
	{"VCF_getBial" , (DL_FUNC) & VCF_getBial , 2 },		//VCF_getBial( SEXP vcfptr, SEXP mat )
	
	//
	//
	{"BGZF_compress" , (DL_FUNC) & BGZF_compress , 2 },						//BGZF_compress( SEXP in_filename_sexp, SEXP out_filename_sexp )
	
	//	Single element read functions
	//
	{"VCF_parseNextSNP" , (DL_FUNC) & VCF_parseNextSNP , 1 },			//VCF_parseNextSNP( SEXP vcfptr )
	{"VCF_parseNextLine" , (DL_FUNC) & VCF_parseNextLine , 1 },		//VCF_parseNextLine( SEXP vcfptr )
	{"VCF_getChrom" , (DL_FUNC) & VCF_getChrom , 1 },					//VCF_getChrom( SEXP vcfptr )
	{"VCF_getPos" , (DL_FUNC) & VCF_getPos , 1 },						//VCF_getPos( SEXP vcfptr )
	{"VCF_getID" , (DL_FUNC) & VCF_getID , 1 },						//VCF_getID( SEXP vcfptr )
	{"VCF_getRef" , (DL_FUNC) & VCF_getRef , 1 },						//VCF_getRef( SEXP vcfptr )
	{"VCF_getAlt" , (DL_FUNC) & VCF_getAlt , 1 },						//VCF_getAlt( SEXP vcfptr )
	{"VCF_getQual" , (DL_FUNC) & VCF_getQual , 1 },					//VCF_getQual( SEXP vcfptr )
	{"VCF_getFilter" , (DL_FUNC) & VCF_getFilter , 1 },				//VCF_getFilter( SEXP vcfptr )
	{"VCF_getInfo" , (DL_FUNC) & VCF_getInfo , 1 },					//VCF_getInfo( SEXP vcfptr )
	{"VCF_getInfoField" , (DL_FUNC) & VCF_getInfoField , 2 },			//VCF_getInfoField( SEXP vcfptr, SEXP fieldnam )
	{"VCF_getFormat" , (DL_FUNC) & VCF_getFormat , 1 },				//VCF_getFormat( SEXP vcfptr )
	{"VCF_getSampleByName" , (DL_FUNC) & VCF_getSampleByName , 2 },	//VCF_getSampleByName( SEXP vcfptr, SEXP str_samplename )
	{"VCF_getSample" , (DL_FUNC) & VCF_getSample , 2 },				//VCF_getSample( SEXP vcfptr, SEXP stridx )
	
	//	Line readers
	//
	{"VCF_readLineRaw" , (DL_FUNC) & VCF_readLineRaw , 2 },					//VCF_readLineRaw( SEXP vcfptr , SEXP str=NA )
	{"VCF_readLineRawFiltered" , (DL_FUNC) & VCF_readLineRawFiltered , 2 },	//VCF_readLineRawFiltered( SEXP vcfptr, SEXP str=NA )
	{"VCF_readLineTSV" , (DL_FUNC) & VCF_readLineTSV , 1 },					//VCF_readLineTSV( SEXP vcfptr )
	{"VCF_readLineTSVFiltered" , (DL_FUNC) & VCF_readLineTSVFiltered , 1 },	//VCF_readLineTSVFiltered( SEXP vcfptr )
#if 0
	{"VCF_readLineDF" , (DL_FUNC) & VCF_readLineDF , 1 },						//VCF_readLineDF( SEXP vcfptr )
	{"VCF_readLineDFFiltered" , (DL_FUNC) & VCF_readLineDFFiltered , 1 },		//VCF_readLineDFFiltered( SEXP vcfptr )
#endif

	//--
	{NULL, NULL, 0}
};




EXPORT void /*attribute_visible*/ R_init_WhopGenome(DllInfo *info)
{
	/* Register routines, allocate resources. */
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, (Rboolean)FALSE);
}

EXPORT void /*attribute_visible*/ R_unload_WhopGenome(DllInfo *info)
{
	/* Release resources. */
}



