/*
**
**		WhopGen
**
**		WHOle genome population genetics with PopGENome
**
**
**		Whopgen - Reading data into matrices
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



//*
//*			DATA
//*




//*
//*			EXTERNS
//*

extern	SEXP	nucleotide_A;
extern	SEXP	nucleotide_C;
extern	SEXP	nucleotide_G;
extern	SEXP	nucleotide_T;
extern	SEXP	nucleotide_N;

extern	char	nucleotide_mapping[];


//*
//*			CODE
//*



/*! Reads SNP lines from the VCF and processes them using the given functor
**
**	- 	Writes data directly into the given R-matrix in <mat>
**
**
*/
EXPORT SEXP helper_read_intmatrix_functored( SEXP vcfptr, MatrixLoaderBaseClass &fn, SEXP mat )
{

	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("helper_read_intmatrix_functored :: Parameter not a VCFhandle EXTPTR!\n");
		return RBool::False();
	}
	fn.vcf = f;
	
	df2("A\n");
	
	//	any samples selected ?
	//
	if( f->num_wanted_samples < 1 )
	{
		Rprintf("helper_read_intmatrix_functored :: No samples selected!\n");
		return RBool::False();
	}
	
	df2("B\n");

	//	without a FORMAT field, it is not clear how to discover the genotype subfield per sample!
	//
	unsigned int	samplefieldindex = f->getFirstSampleFieldIndex();
	if( samplefieldindex <= FORMAT )
	{
		Rprintf("helper_read_intmatrix_functored :: VCF does not have a FORMAT field!\n");
		return RBool::False();
	}
	
	df2("C\n");
	
	//	is mat a matrix?
	//
	RMatrix m(mat);
	if( false == m.isValid() )
	{
		Rprintf("helper_read_intmatrix_functored :: Second parameter is not an integer matrix!\n");
		return RBool::False();
	}
	
	df2("D\n");

	//	biallelic matrices are integer-typed
	//		save 4 bytes per entry over standard 'double' datatype
	//
	if( m.getType() != INTSXP )
	{
		df0("helper_read_intmatrix_functored :: Parameter not a integer matrix!\n");
		return RBool::False();
	}
	
	df2("E\n");
	
	
	//	enough rows to hold all samples?
	//
	unsigned int	nrow = m.numRows();
	if( f->num_wanted_samples > (unsigned)nrow )
	{
		Rprintf("helper_read_intmatrix_functored :: %d samples selected but matrix offers only rows for %d samples!\n",f->num_wanted_samples,nrow);
		return  RBool::False();
	}
	
	//
	//
	//
	int				nonbialcols=0,
					bialcols=0
						;
	unsigned int	ncol = m.numCols();
	int*		 	ptr = m.getIntPtr();
	
	//
	//
	if( ptr == 0 )
	{
		Rprintf("helper_read_intmatrix_functored :: ERROR : Could not get access to the matrix in form of an int*!\n");
		return  RBool::False();
	}
	
	//
	char			textbuffer[256];		//for sprintf()-ing SNP positions for matrix-column names
	char			*fieldptr=0;
	unsigned int	per_column = 0;		//vars here to find out what to clear
	unsigned int	per_row = 0;		//	when too little data exists
	unsigned int	column_stepsize = nrow;
	int				snppos=-1;
	
	SEXP colnamvec = m.getColNames();
	if( R_NilValue == colnamvec )
		Rprintf("helper_read_intmatrix_functored :: WARNING : matrix has no column names vector! not setting SNP positions in matrix!\n");

	//
	//
	for( ; per_column < ncol ; per_column ++ )
	{
		
		//-
		//
		//	Get a valid biallelic SNP line from the file
		//
		//-

		//char	*refptr=0;
		//char	*altptr=0;
		fn.refptr=0;
		fn.altptr=0;
		if( false == fn.findNextLine( f, &fn.refptr, &fn.altptr ) )
		{
			//df1("No more lines!\n");
			//break;
			return RBool::False();
		}

		//--------------
		//
		//	POSTCONDITION : got a line from VCF and it is a valid SNP
		//
		//--------------
		
		//
		fieldptr = (char*)f->getFieldPtr( POS );
		if( 0 == fieldptr )
		{
			df1("POS field not found!\n");
			//break;
			return RBool::False();
		}
		
		//
		snppos = atoi( fieldptr );
		if( 0 >= snppos )
		{
			df1("POS field specifies 0 or negative value!\n");
			//break;
			return RBool::False();
		}

		//--------------
		//
		//	POSTCONDITIONS :
		//		- got a line from VCF and it is a valid SNP
		//		- successfully identified the GT subfield
		//
		//--------------
		
		
		//	check whether filtering needs to be done
		//
		if( /*bFilterActivated && */(false == filterLine( f )) )
		{
			df1("helper_read_intmatrix_functored :: Line (pos %d)has been filtered away\n",snppos);
			per_column--;
			continue;
		}

		//-
		//
		//	for each selected individual, get the SNP genotype information
		//		and store a 0 or 1 in the biallelic matrix
		//
		//-
		fn.bHadRef = false;
		fn.bHadAlt = false;
		
		//
		for( per_row = 0; per_row < f->num_wanted_samples ; per_row ++ )
		{
			//	- get field of sample
			//
			fieldptr = (char*)f->getFieldPtr( f->wanted_samples[per_row] );
			if( fieldptr == 0 )
			{
				Rprintf("helper_read_intmatrix_functored :: ERROR when trying to get sample %d (matrix row %d) in file!\n",f->wanted_samples[per_row],per_row);
				Rprintf("	per_row =%d\nwanted_sample[per_row]=%d\n",per_row, f->wanted_samples[per_row] );
				Rprintf("	baseindex=%d, field = %d\n",samplefieldindex, (samplefieldindex + f->wanted_samples[per_row]) );
				Rprintf("	numparsedfields=%d\n",f->numParsedFields());
				return RBool::False();
			}
			
			/*
			**	- parse the GT subfield ( regexp: [0-9]+[\/\|][0-9]+ )
			*/
			if( false == fn.processSampleGTi( fieldptr, &ptr[per_row] ) )
			{
				Rprintf("(!!) Error during sample processing (pos=%d,matrix col=%d,row=%d): '%s'\n",snppos,per_column,per_row, &fn.errormessage[0] );
				//break;
				return RBool::False();
			}

			//
		}//...for( all rows == all samples )
		
		
		//----------------------------------------
		//
		//	clear any unused matrix rows of this column, if it had both alleles
		//
		if( true==fn.bHadAlt && true==fn.bHadRef )
		{
			
			//	fill unused matrix elements in current column
			//
			for( ; per_row < nrow ; per_row ++ )
			{
				ptr[per_row] = -2;
			}
			
			//
			ptr += column_stepsize;
			
			bialcols++;
			
			//	Set Column Names = SNP Positions
			//
			if( R_NilValue != colnamvec )
			{
				snprintf(textbuffer,sizeof(textbuffer)-2,"%d",snppos);
				SET_STRING_ELT( colnamvec, per_column, mkChar(textbuffer) );
			}
		
			//FIXME any rows beyond the ones we need : zero them out here
			//	or as a separate step after this double for-loop
		}
		//
		//	if the last VCF-line's individuals did not provide both alleles
		//		reuse current column and try the next line
		//
		else
		{

			if( R_NilValue != colnamvec )
			{
				SET_STRING_ELT( colnamvec, per_column, fn.EmptyMatrixColumnNameSEXP );
			}
			
			per_column--;	//re-use the matrix column for the next SNP
			nonbialcols++;	//count how many non-biallelic columns we found
		}
		//
		//
		//
		//----------------------------------------
		
		//
	}//...for( each column in the matrix )

	return ( (0<per_column)?RBool::True():RBool::False() );
}




			//
			//
			//
			//
			//				N u m e r i c   M a t r i x   R e p r e s e n t a t i o n s
			//
			//
			//
			//






/*!	Write into a SEXP matrix biallelic SNPS as multi-digit nucleotide codes (two-digit for diploid)
**		- see mapping from ASCII to nucleotide codes in file w_codemat.cpp (ca. line 71)
**
**
*/
EXPORT SEXP read_snp_anyploid_multiallelic_int_nuclcodes( SEXP vcfptr, SEXP mat )
{
	
	class MatrixLoader_AnyploidMultiallelicSNP_Int_NucleotideCodes : public MatrixLoaderBaseClass
	{
	public:
	
		MatrixLoader_AnyploidMultiallelicSNP_Int_NucleotideCodes()
		{
			EmptyMatrixCellSEXP = R_NilValue;
			EmptyMatrixCellInteger = 66;
		}

		//-
		//-		VIRTUAL METHODS
		//-
		
		//!
		virtual		bool	findNextLine( vcff * f, char** refptr, char** altptr )
		{
			bool r = findNextLine_AnySNP( f, refptr, altptr );
			return r;
		}
		
		
		//!
		virtual		bool	processSampleGTs( const char * samplePtr, SEXP*ptr )
		{
			return false;
		}
		
		//!
		virtual		bool	processSampleGTi( const char * samplePtr, int*ptr )
		{
//			Rprintf("ENTER (%s)\n",samplePtr);
			//	read all genotypes into a multi-digit code
			//
			int 	gtcode = 0,						// stores numeric code to write into the matrix for this sample 
					allelenum = 0,					// temporary variable for turning the ASCII allele number to int
					ismultiallelicposition = 0		// to find out whether this is an informative position (i.e. more than one allele present in at least one of the selected samples)
						;
			
			char	allelenuc = '?',
					*tmp
						;
			
			//
			const char* fcopy = samplePtr;		// string of the form "0/1/1/0/2" or "1|0" or "0|1/1"
			
			//
			//	process entire GT field
			//
			while( fcopy[0] != 0 && fcopy[0] != ':' && fcopy[0] != '\t' )
			{
//				Rprintf("(%c)",fcopy[0]);
				allelenum=0;
				
				//
				//	read allele number (supporting numbers > 9)
				//
				while( fcopy[0] >= '0' && fcopy[0] <= '9' )
				{
					
					// simple ASCII-decimal-number-string to C++ int conversion
					allelenum = (allelenum * 10) + (fcopy[0] - '0');
					
					//
					ismultiallelicposition |= (1 << allelenum);
					
					fcopy++;

				}//...while character is not a in-GT-field chromosome-divider ( | phased  / unphased or \ obsolete )
				

				//
				//	end of allele number is marked by one of \t | \ / or :
				//
				if( fcopy[0] == '|' || fcopy[0] == '/' || fcopy[0] == '\\'  )
				{
					fcopy++;
				}
				else if( ! is_FORMAT_subfield_delimiter( fcopy[0] ) )	//fcopy[0] != '\t' && fcopy[0] != '\n' && fcopy[0] != ':' && fcopy[0] != 0 )
				{
					//Rprintf("ERROR : unexpected character '%c' in Genotype field at position %d\n",fcopy[1],snppos);
					snprintf( errormessage, sizeof(errormessage)/sizeof(errormessage[0]), "ERROR : unexpected character '%c' in Genotype (GT) field\n", fcopy[1] );
					df0("	=> Syntax error in GT field (%s)!\n",samplePtr);
					return false;
				}//...if character is not a decimal digit
				
//				Rprintf("(AN=%d)",allelenum);
				

				//	turn indices into the REF/ALT fields into the actual nucleotide-codes
				//
				if( allelenum == 0 )
					allelenuc = *refptr;
				else
				{
					tmp=altptr;
					while( --allelenum > 0 )
					{
						tmp+=2;	//skip nucleotide code and separator (,)
					}
					allelenuc = *tmp;
				}
//				Rprintf("(AC=%c)",allelenuc);
				gtcode = ( gtcode * 10 ) + nucleotide_mapping[ allelenuc&0xFF ];

			}//...while character is not end of GT-field

			// write numerical genotype code
			//
			*ptr = gtcode;

			//
			//
			if( ismultiallelicposition&-2 )
			{
				bHadAlt = true;	//memorise that we got samples with alternate alleles
			}
			
			if( ismultiallelicposition&1 )
			{
				bHadRef = true;	//memorise that we got samples with reference alleles
			}
			
//			Rprintf("(MA=%02x) [ha=%c,hr=%c][gtcode=%d]\n",ismultiallelicposition,bHadAlt?'Y':'n',bHadRef?'Y':'n', gtcode);
			
			return true;

		}//...processSampleGT


	} f ;
	//...class MatrixLoader_AnyploidMultiallelicSNP_Int_NucleotideCodes
	
	return helper_read_intmatrix_functored( vcfptr, f, mat );

}








