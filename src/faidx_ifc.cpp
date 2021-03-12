/*
**
**		WhopGen
**
**		WHOle genome population genetics with PopGENome
**
**
**		Loading FaIdx indexed fasta files
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
#include	"tabix/faidx.h"




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


class faifile {
public:
	
	//-
	//-		CTORs + DTOR
	//-

	faifile() : faidx(0), num_sequences(0) {
		ONDBG Rprintf("(!!) fai CTOR\n");
	}
	faifile(const char * filename ){
		ONDBG Rprintf("(!!) fai CTOR '%s'\n",filename);
		open(filename);
	}
	~faifile(){
		ONDBG Rprintf("(!!) fai DTOR\n");
		close();
	}
	
	//-
	//-		METHODS
	//-
	
	//

	//!
	bool			open( const char * filename, bool autocreate_index=false )
	{
		ONDBG Rprintf("(!!) opening fai '%s'\n",filename);
		faidx = fai_load(filename);
		if( faidx == 0 )
		{
			ONDBG Rprintf("(!!) faiload fail first-chance\n");
			if(
				false == autocreate_index
				|| (false ==  create(filename) )
				|| (0 == (faidx=fai_load(filename)))
			  ) 
			{
				ONDBG Rprintf("(!!) faiload fail final\n");
				return false;//true;
			}
		}
		
		ONDBG Rprintf("(!!) faiload success\n");
		num_sequences = faidx_fetch_nseq(faidx);
		ONDBG Rprintf("(!!) faiload %d seqs\n",num_sequences);
		return true;
	}
	
	//!
	bool			create( const char * filename )
	{
		ONDBG Rprintf("(!!) building fasta-index for '%s'\n",filename);
		int r = fai_build(filename);
		if( r < 0 )
		{
			ONDBG Rprintf("(!!) failed to build fasta index for (%s)!\n",filename);
			return false;
		}
		return true;
	}
	
	//!
	void			close( void )
	{
		ONDBG Rprintf("(!!) closing fasta-index\n");
		if( faidx )
		{
			fai_destroy( faidx );
			faidx=0;
			num_sequences=0;
		}
	}
	
	//!
	bool			isValid( void )
	{
		return( (0!=faidx) && (0<num_sequences) );
	}

	//
	//
	faidx_t			*faidx;
	unsigned int	num_sequences;

private:

	//
};


//*
//*			DATA
//*


SEXP		faihandle_attrname_filename = R_NilValue;


//*
//*			EXTERNS
//*

EXPORT SEXP	FAI_close( SEXP faiptr );

//*
//*			CODE
//*


//!
SEXP _internal_FaiGetAttrFilename( void )
{
	if( faihandle_attrname_filename == R_NilValue )
	{
		faihandle_attrname_filename = install("FAI.filename");
	}
	return faihandle_attrname_filename;
}



/*!
**
**
**
*/
static void fai_finalize(SEXP extPtr)
{
	df1("FAIFILE FINALIZE!\n");
	faifile* f = (faifile*)R_GetExtPtr( extPtr, "FAIhandle" );
	if( f )
	{
		df1("fai_finalize : Finalizing FAIhandle!\n");
		FAI_close( extPtr );
		df1("fai_finalize : Successfully finalized FAIhandle!\n");
	}
	else
	{
		df1("fai_finalize : Could not finalize potential FAIhandle!\n");
	}

	//
}





/*!	
**
**
*/
EXPORT SEXP	FAI_query4( SEXP faiptr, SEXP seqn, SEXP beg, SEXP end )
{
	//
	faifile * f = (faifile*)R_GetExtPtr( faiptr , "FAIhandle" );
	if( 0 == f )
	{
		df0("FAI_query4 : parameter 1 is not a FAIhandle or nil!\n");
		return RBool::False();
	}

	//
	if( false == RString::isStr( seqn ) )
	{
		Rprintf("FAI_query4 : argument 2, 'seqname', is not a string!\n");
		return RBool::False();
	}

	//
	//
	const char* seqnam = RString::get(seqn);
	int frompos = RNumeric::getInt(beg);
	int topos = RNumeric::getInt(end);
	if( frompos <= 0 || topos <= 0 )
	{
		Rprintf("FAI_query : unexpected values for parameters 3, start(%d), and 4, end(%d)\n",frompos,topos);
		return RBool::False();
	}
	
	//
	int seqlen = 0;
	char * fastaseq = faidx_fetch_seq(f->faidx, (char*)seqnam, frompos, topos, &seqlen);
	if( 0 == fastaseq )
	{
		//ONDBG Rprintf("FAI_query : no result\n");
		//SET_STRING_ELT( resstr, 0, mkChar("") );
		return RBool::False();
	}
	
	//
	RString res(fastaseq);
	free( fastaseq );

	//
	return res.get();//RBool::True();
}





/*!	
**
**
*/
EXPORT SEXP	FAI_query2( SEXP faiptr, SEXP regionstr )
{
	//
	faifile * f = (faifile*)R_GetExtPtr( faiptr , "FAIhandle" );
	if( 0 == f )
	{
		df0("FAI_query2 : parameter 1 is not a FAIhandle or nil!\n");
		return RBool::False();
	}

	//
	if( false == RString::isStr( regionstr ) )
	{
		Rprintf("FAI_query2 : argument 2, 'regionstr', is not a string!\n");
		return RBool::False();
	}

	//
	//
	const char* regstr = RString::get(regionstr);
	
	//
	int seqlen = 0;
	const char * fastaseq = fai_fetch(f->faidx, (char*)regstr, &seqlen);
	if( 0 == fastaseq )
	{
		//ONDBG Rprintf("FAI_query : no result\n");
		//SET_STRING_ELT( resstr, 0, mkChar("") );
		return RBool::False();
	}
	
	//
	
	//
	RString res(fastaseq);
	free( (void*)fastaseq );

	//
	return res.get();
	
	//SET_STRING_ELT( resstr, 0, mkChar(fastaseq) );
	//return RBool::True();
}





/*!	
**
**
*/
EXPORT SEXP	FAI_build( SEXP filename )
{

	//
	if( false == RString::isStr( filename ) )
	{
		Rprintf("FAI_build : argument 1, 'filename', is not a string!\n");
		return RBool::False();
	}

	//
	//
	const char* filename_cstr = RString::get(filename);
	if( 0 == filename_cstr )
	{
		Rprintf("FAI_query : argument 5, 'resultstring', is not a string!\n");
		return RBool::False();
	}
	
	//
	int r = fai_build(filename_cstr);
	if( r < 0 )
	{
		Rprintf("(!!) failed to build fasta index for (%s)!\n",filename_cstr);
		return RBool::False();
	}

	//
	return RBool::True();
}




/*!	Open a Faidx-indexed fasta file and return a FAIhandle, a whopgen-managed EXTPTR SEXP
**
**
*/
EXPORT  SEXP FAI_open( SEXP filename )
{
	//
	if(! CHKSTR(filename,1) )
	{
		df0("FAI_open : filename is not a single string!");
		return R_NilValue;
	}
	
	//-----
	
	//open FAI with filename
	//	if failed: return NULL
	//
	faifile * f = new faifile( CHAR(STRING_ELT(filename, 0)) );
	if( f == 0 )
	{
		df0("FAI_open: Could not open indexed fasta file '%s'!\n", CHAR(STRING_ELT(filename, 0)));
		return R_NilValue;
	}
	
	//
	if( false == f->isValid() )
	{
		delete f;
		f=0;
		df0("FAI_open : Could not open file '%s' as tabix-indexed!\n",CHAR(STRING_ELT(filename, 0)));
		return R_NilValue;
	}

	//
	df1("(FAI_open) opened file '%s' is a FAI!\n",CHAR(STRING_ELT(filename, 0)));
	
	//
	//
	SEXP res;
	PROTECT(
		res = R_MakeExternalPtr( f, install("FAIhandle"), R_NilValue)
	);
	if( res == R_NilValue )
	{
		df0("FAI_open : could not create external pointer SEXP!\n");
		return res;
	}
	
	//
	R_RegisterCFinalizerEx(res, fai_finalize, Rboolean_TRUE);
	
	//
	setAttrib( res, _internal_FaiGetAttrFilename(), filename );
	
	UNPROTECT( 1 );
	return res;

}



/*!	Close a FAIhandle, an EXTPTR SEXP of a whopgen-managed faidx-indexed FASTA file 
**
**
*/
EXPORT SEXP	FAI_close( SEXP faiptr )
{
	//
	faifile * f = (faifile*)R_GetExtPtr( faiptr , "FAIhandle" );
	if( 0 == f )
	{
		df0("FAI_close : parameter is not a FAIhandle or nil!\n");
		return RBool::False();
	}
	
	//
	R_ClearExternalPtr(faiptr);
	
	//
	delete f;
	
	//
	return RBool::True();
}





/*!	Reopens a FAI file if the FAIhandle got stale
**
**
*/
EXPORT SEXP	FAI_reopen( SEXP faiptr )
{
	//
	if( false == RType::IsExtPtr( faiptr ) || ( strcasecmp( RExtPtr::getTag(faiptr), "FAIhandle" ) != 0 ) )
	{
		df1("FAI_reopen : parameter is not an externalptr FAIhandle!\n");
		return RBool::False();
	}
	
	//
	faifile * f = (faifile*)R_GetExtPtr( faiptr , "FAIhandle" );
	if( f != 0 )
	{
		return RBool::True();
	}
	
	//
	SEXP filename = getAttrib( faiptr, _internal_FaiGetAttrFilename() );
	
	//
	//open fai with filename
	//	if failed: return False
	//
	f = new faifile( CHAR(STRING_ELT(filename, 0)) );
	if( f == 0 )
	{
		df0("FAI_reopen : Could not open file '%s' as faidx-indexed!\n",CHAR(STRING_ELT(filename, 0)));
		return RBool::False();
	}
	
	//	if loading failed, return with error
	//
	if( false == f->isValid() )
	{
		delete f;
		f=0;
		df0("FAI_reopen : Could not open file '%s' as faidx-indexed!\n",CHAR(STRING_ELT(filename, 0)));
		return RBool::False();
	}

	R_SetExternalPtrAddr( faiptr, f );
	
	//
	return RBool::True();
}


//---------------------------------------------------------			DEPRECATED:


/*!	
**
**
*/
EXPORT SEXP	FAI_query5( SEXP faiptr, SEXP seqn, SEXP beg, SEXP end, SEXP resstr )
{
	//
	faifile * f = (faifile*)R_GetExtPtr( faiptr , "FAIhandle" );
	if( 0 == f )
	{
		df0("FAI_query5 : parameter 1 is not a FAIhandle or nil!\n");
		return RBool::False();
	}

	//
	if( false == RString::isStr( seqn ) )
	{
		Rprintf("FAI_query5 : argument 2, 'seqname', is not a string!\n");
		return RBool::False();
	}

	//
	if( false == RString::isStr( resstr ) )
	{
		Rprintf("FAI_query5 : argument 5, 'resstr', is not a string!\n");
		return RBool::False();
	}

	//
	//
	const char* seqnam = RString::get(seqn);
	int frompos = RNumeric::getInt(beg);
	int topos = RNumeric::getInt(end);
	if( frompos <= 0 || topos <= 0 )
	{
		Rprintf("FAI_query : unexpected values for parameters 3, start(%d), and 4, end(%d)\n",frompos,topos);
		return RBool::False();
	}
	
	//
	int seqlen = 0;
	char * fastaseq = faidx_fetch_seq(f->faidx, (char*)seqnam, frompos, topos, &seqlen);
	if( 0 == fastaseq )
	{
		//ONDBG Rprintf("FAI_query : no result\n");
		SET_STRING_ELT( resstr, 0, mkChar("") );
		return RBool::False();
	}
	
	//
	SET_STRING_ELT( resstr, 0, mkChar(fastaseq) );
	free( fastaseq );

	//
	return RBool::True();
}





/*!	
**
**
*/
EXPORT SEXP	FAI_query3( SEXP faiptr, SEXP regionstr, SEXP resstr )
{
	//
	faifile * f = (faifile*)R_GetExtPtr( faiptr , "FAIhandle" );
	if( 0 == f )
	{
		df0("FAI_query3 : parameter 1 is not a FAIhandle or nil!\n");
		return RBool::False();
	}

	//
	if( false == RString::isStr( regionstr ) )
	{
		Rprintf("FAI_query3 : argument 2, 'regionstr', is not a string!\n");
		return RBool::False();
	}

	//
	if( false == RString::isStr( resstr ) )
	{
		Rprintf("FAI_query3 : argument 3, 'resstr', is not a string!\n");
		return RBool::False();
	}

	//
	//
	const char* regstr = RString::get(regionstr);
	
	//
	int seqlen = 0;
	const char * fastaseq = fai_fetch(f->faidx, (char*)regstr, &seqlen);
	if( 0 == fastaseq )
	{
		//ONDBG Rprintf("FAI_query : no result\n");
		SET_STRING_ELT( resstr, 0, mkChar("") );
		return RBool::False();
	}
	
	//
	SET_STRING_ELT( resstr, 0, mkChar(fastaseq) );
	free( (void*)fastaseq );	
	return RBool::True();
}
