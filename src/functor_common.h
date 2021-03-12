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
*/

//*
//*			INCLUDES
//*

#include	"w_common.h"

//*
//*			DEFINES
//*



//*
//*			STRUCTS
//*



//*
//*			CLASSES
//*


/*! Abstract base class of all functors used to load matrix representations from VCF files
*/
class MatrixLoaderBaseClass
{
public:
	MatrixLoaderBaseClass() : vcf(0) {
		EmptyMatrixColumnNameSEXP = mkChar("unused");
		EmptyMatrixCellSEXP = R_NilValue;
		EmptyMatrixCellInteger = -1;
	}

	//-
	//-		DATA
	//-

	//
	char		errormessage[256];
	
	//
	char		*refptr;
	char		*altptr;
	
	//
	bool		bHadRef,			//whether processSampleGT_ produced a meaningful column of differing alleles
				bHadAlt				//		(dependent on sample selection, this may not always be the case)
					;

	//
	SEXP		EmptyMatrixColumnNameSEXP;
	SEXP		EmptyMatrixCellSEXP;	//what to fill in for unused elements in the result matrix
	int			EmptyMatrixCellInteger;
	
	
	vcff		*vcf;					// set by caller : vcff class instance we're working with
	
	
	//-
	//-		VIRTUAL METHODS
	//-
	
	virtual		bool		findNextLine( vcff * f, char** refptr, char** altptr )=0;
	
	virtual		bool		processSampleGTi( const char * samplePtr, int* resptr )=0;	
	virtual		bool		processSampleGTs( const char * samplePtr, SEXP* resptr )=0;
	
	//-
	//-		FIXED METHODS
	//-


	//-
	//-		default impl to find a biallelic SNP line
	//-
	//-
	bool		findNextLine_BiallelicSNP( vcff * f, char** refptr, char**altptr )
	{
	
		bool	bLineParsed=false;
		
		//
		//	read lines until a SNP line is found (or the chromosomal region's end is reached)
		//
		//while( false == (bLineParsed = f->parseNextLine()) )//
		bLineParsed = f->parseNextLine();
		while( bLineParsed )
		{
			//
			//
			/*
				- make sure its a bi-allelic SNP line
				* 
				* REFerence is a single nucleotide in case of a SNP
				* ALTernative is also a single nucleotide in case of a SNP
				*
			*/
			*refptr = (char*)f->getFieldPtr( REF );		//REFerence allele
			*altptr = (char*)f->getFieldPtr( ALT );		//ALTernative allele
			
			if( *refptr == 0 )		return false;
			
			//second char must be a \tab to ensure there is only a single REFerence nucleotide
			//
			if( *refptr && (*refptr)[1] == '\t' )
			{
				//second char must be a \tab to ensure there is only a single ALTernative nucleotide
				//
				if( *altptr && (*altptr)[1] == '\t' )
				{
					
					//
					//	make sure the line contains a GT field ( as first subfield in the SAMPLE-columns as per VCF specifications )
					//
					char * formatptr = (char*)f->getFieldPtr( FORMAT );
					if( formatptr[0] == 'G' && formatptr[1] == 'T' && ( formatptr[2] == ':' || formatptr[2] == '\t' ) )
					{
						break;
					}
					
					//this line's FORMAT field does either
					//			a) NOT specify a GT field at all [ => allowed according to VCF specifications ]
					//	or		b) NOT specify it as first subfield in the SAMPLE-columns [ => VCF specifications violation ]
					//
					// ... so continue looking for valid lines
						
				}
				//else
				//{
				//	//ALTernative allele is something else than a single nucleotide
				//	//df1("Not a biallelic SNP! (REF=%9s)\n",refptr);
				//}
			}
			//else
			//{
			//	//REFerence is not something else than a single nucleotide
			//	//df1("Not a SNP! (REF=%2s)\n",refptr);
			//}
			
			
			//	try next line
			//
			bLineParsed = f->parseNextLine();

		}//while( could read another line from VCF )
		
		return bLineParsed;

	}//...findNextLine_BiallelicSNP()
	
	
	//-
	//-		default impl to find SNPs without regard for the number of alleles
	//-			(actual POLY-morphism instead of dimorphism)
	//-
	bool		findNextLine_AnySNP( vcff * f, char** refptr, char**altptr )
	{
	
		bool	bLineParsed=false;
		
		//
		//	read lines until a SNP line is found (or the chromosomal region's end is reached)
		//
		while( (bLineParsed = f->parseNextLine()) == true )
		{
			
			// make sure we have a single-nucleotide reference allele
			//
			*refptr = (char*)f->getFieldPtr( REF );		//REFerence allele
			//Rprintf("fnany: ref(%s)\n",*refptr );
			if( *refptr==0 || (*refptr)[1] != '\t' )
			{
				continue;
			}
			
			// make sure all alternative alleles are single nucleotides
			//
			*altptr = (char*)f->getFieldPtr( ALT );		//ALTernative allele
			//Rprintf("fnany: alt(%s)\n",*altptr );
			if( *altptr == 0 )
			{
				continue;
			}
			
			int i=1;
			//Rprintf("fnany1: (%c)\n",(*altptr)[i] );
			while( (*altptr)[i] == ',' )
			{
				i+=2;
				//Rprintf("fnanyX: (%c)\n",(*altptr)[i] );
			}
						
			if( (*altptr)[i] != '\t' )
			{
				//Rprintf("fnanyN: (%c) -- CONT\n",(*altptr)[i] );
				continue;
			}
			
			
			//Rprintf("fnanyE: (%c) OK\n",(*altptr)[i] );
			
			break;

		}//while( could read another line from VCF )
		
		
		//Rprintf("fnany: (%s)\n",bLineParsed?"OK":"NO" );
		return bLineParsed;

	}//...findNextLine_AnySNP()
	
};//...class MatrixLoaderBaseClass



//*
//*			DATA
//*


//*
//*			CODE
//*

//*
//*			EXTERNS
//*
