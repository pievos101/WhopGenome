/*
**
**		WhopGenome
**
**		WHOle genome population genetics
**
**
**
**
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



//*
//*			CODE
//*




//!
inline	bool	ALT_is_snp( const char * ptr )
{
	ptr++;
	while( *ptr != '\t' )
	{
		if( *ptr != ',' )
			return false;
		ptr+=2;
	}
	return true;
}




			//
			//
			//
			//
			//
			//
			//
			//
			//
			//
			//




/*
**
**
*/
bool			snpmat_init_validate_int( snpmat_read_info_int &inf, SEXP vcfptr, SEXP mat )
{
//	DBG( "[%s]\n" , __FUNCTION_NAME__ );

	inf.f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == inf.f )
	{
		RPRINT_ERROR( "Parameter not a VCFhandle EXTPTR!" );
		return false;
	}
	
	//	any samples selected ?
	//
	if( inf.f->num_wanted_samples < 1 )
	{
		RPRINT_ERROR( "No samples selected!" );
		return false;
	}

	//	without a FORMAT field, it is not clear how to discover the genotype subfield per sample!
	//
	unsigned int	samplefieldindex = inf.f->getFirstSampleFieldIndex();
	if( samplefieldindex <= FORMAT )
	{		
		RPRINT_ERROR( "VCF does not have a FORMAT field!" );
		return false;
	}
	
	//	is mat a matrix?
	//
	RMatrix m(mat);
	if( false == m.isValid() )
	{
		RPRINT_ERROR( "Second parameter is not an integer matrix!" );
		return false;
	}

	//	biallelic matrices are integer-typed
	//		save 4 bytes per entry over standard 'double' datatype
	//
	if( m.getType() != INTSXP )
	{
		RPRINT_ERROR( "Parameter not a integer matrix!" );
		return false;
	}
	
	
	//	enough rows to hold all samples?
	//
	inf.nrow = m.numRows();
	if( inf.f->num_wanted_samples > (unsigned)inf.nrow )
	{
		RPRINT_ERROR( "%d samples selected but matrix offers only rows for %d samples!\n",inf.f->num_wanted_samples,inf.nrow );
		return false;
	}
	
	//
	//
	inf.ncol = m.numCols();
	inf.ptr = m.getIntPtr();					//VARIABLE PART
	
	//
	//
	//
	if( 0 == inf.ptr )
	{
		RPRINT_ERROR( "Could not get access to the matrix in form of an int*!" );					//VARIABLE PART in MESSAGE
		return false;
	}
	
	//
	//
	//
	inf.colnamvec = m.getColNames();
	if( R_NilValue == inf.colnamvec )
	{
		//Rprintf("helper_read_intmatrix_functored :: WARNING : matrix has no column names vector! not setting SNP positions in matrix!\n");
		RPRINT_WARN( "matrix has no column names vector! not setting SNP positions in matrix!\n" );
	}

	return true;

}///..snpmat_init_validate_int







			//
			//
			//
			//
			//			INT,			FILTERED
			//
			//
			//
			//
			//
			//






/*!	Run the main SNP matrix read loop
**	+ FILTERING
**	+ BIALLELIC SNPS ONLY
**	+ INTEGER MATRIX
**
**
**
*/
bool		snpmat_run_loop_int( snpmat_read_info_int& inf )
{
//	DBG( "[%s]\n" , __FUNCTION_NAME__ );

	//
	char			*fieldptr=0;	
	unsigned int	column_stepsize = inf.nrow;

	//
	//
	//
	for( ; inf.cur_col < inf.ncol ; inf.cur_col ++ )
	{
		//-
		//
		//	Get a valid biallelic SNP line from the file
		//
		//-
		const char* refptr;
		const char* altptr;
		bool	bLineParsed;
		while( true == (bLineParsed = inf.f->parseNextLine()) )
		{
			
			refptr = (char*)inf.f->getFieldPtr( REF );
			altptr = (char*)inf.f->getFieldPtr( ALT );

			//	- make sure its a bi-allelic SNP line
			//
			if( refptr && refptr[1] == '\t' )			//TODO use macro/inline
			{
				if( inf.biallelic_locus_required )
				{
					if( altptr && altptr[1] == '\t' )	//both REF and ALT have only single-char alleles?
					{
						break;	//..yes, process this line
					}
				}
				//
				//
				else
				{
					if( ALT_is_snp( altptr ) )		//ALT has any number of comma-separated single-char alleles ?
						break;	// ..yes, process this line
				}
			}///...if REF is single-char allele

		}///...while( could read another line from VCF )

		//
		//
		//	if there are no more lines in the VCF, break the loop early
		//
		//
		if( bLineParsed == false )
		{
			//df1("No more lines!\n");
			break;
		}
		
		//
		//
		//	filter line by rule-set
		//
		//
		if( (true == inf.doFiltering) && (false == filterLine( inf.f )) )
		{
			df1("get_bial :: Line (pos %d) has been filtered away\n",inf.snppos);
			inf.cur_col--;
			continue;
		}

		//
		//
		//	get POS as integer
		//
		//
		fieldptr = (char*)inf.f->getFieldPtr( POS );
		if( 0 == fieldptr )
		{
			RPRINT_ERROR( "Could not find POS column in line!" );
			RPRINT_VCFLINE;
			return false;
		}
		inf.snppos = atoi( fieldptr );

		//-
		//
		//	for each selected individual, get the SNP genotype information
		//		and store a 0 or 1 in the biallelic matrix
		//
		//-
		
		inf.keep_column = 0;		//undecided whether to keep this column or not

		//
		inf.begin_line();
		
		//
		//
		//	for each sample, read genotype data, convert it to integer data and write into matrix
		//
		//
		for( inf.cur_row = 0; inf.cur_row < inf.f->num_wanted_samples ; inf.cur_row ++ )
		{
			
			//	- get field of sample
			//
			fieldptr = (char*)inf.f->getFieldPtr( inf.f->wanted_samples[inf.cur_row] );			
			if( 0 == fieldptr )
			{
				RPRINT_ERROR( "Could not access sample column %d\n" , inf.f->wanted_samples[inf.cur_row] );
				RPRINT_VCFLINE;
				return false;
			}
			
			//
			//
			//
			int process_gt_result = inf.process_sample( fieldptr );
			if( 0 == process_gt_result )
			{
				RPRINT_ERROR( "Failed to process line!" );
				RPRINT_VCFLINE;
				return false;
			}
			

		}///...for all samples
	
		//
		//	clear any unused matrix rows
		//
		if( 1 == inf.end_line() )
		{
			
			//
			for( ; inf.cur_row < inf.nrow ; inf.cur_row ++ )
			{
				inf.ptr[inf.cur_row] = inf.empty_cell_value_int;
			}
			
			//
			inf.ptr += column_stepsize;
			
			
			//	Set Column Names = SNP Positions
			//
			if( R_NilValue != inf.colnamvec )
			{
				char posbuffer[256];
				snprintf(posbuffer,sizeof(posbuffer)-2,"%d",inf.snppos);
				SET_STRING_ELT( inf.colnamvec, inf.cur_col, mkChar(posbuffer) );
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
			//
			inf.cur_col--;
		}

	}///...for( each column in the matrix )
	
	//
	bool completely_filled_matrix = !(inf.cur_col < inf.ncol);
	
	//
	//	Reset/clear unused columns
	//
	//TODO in loop_str, use _str for cell values
	for( unsigned int rcol = inf.cur_col; rcol < inf.ncol; rcol ++ )
	{
		//
		for( unsigned int rrow=0; rrow < inf.nrow ; rrow ++ )
		{
			inf.ptr[rrow] = inf.empty_cell_value_int;
		}
		inf.ptr += column_stepsize;
			
		//	Clear column name
		//
		if( R_NilValue != inf.colnamvec )
		{
			SET_STRING_ELT( inf.colnamvec, rcol, inf.empty_col_name );
		}
	}

	//
	return completely_filled_matrix;


}///...snpmat_run_loop_int




/*!	Run the main SNP matrix read loop
**	+ FILTERING
**	+ BIALLELIC SNPS ONLY
**	+ INTEGER MATRIX
**	+ ONLY DIPLOID SAMPLE GENOTYPES
**
**
*/
bool		snpmat_run_loop_int_diploidonly( snpmat_read_info_int& inf )
{
//	DBG( "[%s]\n" , __FUNCTION_NAME__ );

	//
	char			*fieldptr=0;	
	unsigned int	column_stepsize = inf.nrow;

	//
	//
	//
	for( ; inf.cur_col < inf.ncol ; inf.cur_col ++ )
	{
		//-
		//
		//	Get a valid biallelic SNP line from the file
		//
		//-
		const char* refptr;
		const char* altptr;
		bool	bLineParsed;
		while( true == (bLineParsed = inf.f->parseNextLine()) )
		{
			
			refptr = (char*)inf.f->getFieldPtr( REF );
			altptr = (char*)inf.f->getFieldPtr( ALT );

			//	- make sure its a bi-allelic SNP line
			//
			if( refptr && refptr[1] == '\t' )			//TODO use macro/inline
			{
				if( inf.biallelic_locus_required )
				{
					if( altptr && altptr[1] == '\t' )	//both REF and ALT have only single-char alleles?
					{
						break;	//..yes, process this line
					}
				}
				//
				//
				else
				{
					if( ALT_is_snp( altptr ) )		//ALT has any number of comma-separated single-char alleles ?
						break;	// ..yes, process this line
				}
			}///...if REF is single-char allele

		}///...while( could read another line from VCF )

		//
		//
		//	if there are no more lines in the VCF, break the loop early
		//
		//
		if( bLineParsed == false )
		{
			//df1("No more lines!\n");
			break;
		}
		
		//
		//
		//	filter line by rule-set
		//
		//
		if( (true == inf.doFiltering) && (false == filterLine( inf.f )) )
		{
			df1("get_bial :: Line (pos %d) has been filtered away\n",inf.snppos);
			inf.cur_col--;
			continue;
		}

		//
		//
		//	get POS as integer
		//
		//
		fieldptr = (char*)inf.f->getFieldPtr( POS );
		if( 0 == fieldptr )
		{
			RPRINT_ERROR( "Could not find POS column in line!" );
			RPRINT_VCFLINE;
			return false;
		}
		inf.snppos = atoi( fieldptr );

		//-
		//
		//	for each selected individual, get the SNP genotype information
		//		and store a 0 or 1 in the biallelic matrix
		//
		//-
		
		inf.keep_column = 0;		//undecided whether to keep this column or not

		//
		inf.begin_line();
		
		//
		//
		//	for each sample, read genotype data, convert it to integer data and write into matrix
		//
		//
		for( inf.cur_row = 0; inf.cur_row < inf.f->num_wanted_samples ; inf.cur_row ++ )
		{

			//	- get field of sample
			//
			fieldptr = (char*)inf.f->getFieldPtr( inf.f->wanted_samples[inf.cur_row] );			
			if( 0 == fieldptr )
			{
				RPRINT_ERROR( "Could not access sample column %d\n" , inf.f->wanted_samples[inf.cur_row] );
				RPRINT_VCFLINE;
				return false;
			}

			//	check if sample genotype is diploid
			//
			if( !(is_GT_allele_divider(fieldptr[1])) || !(is_FORMAT_subfield_delimiter(fieldptr[3])) )
			{
				DBG_RPRINT_WARN( "Non-diploid data! (%s)\n" , fieldptr );
				DBG_RPRINT_VCFLINE;
				break;//continue;
			}

			//
			//
			int process_gt_result = inf.process_sample( fieldptr );
			if( 0 == process_gt_result )
			{
				RPRINT_ERROR( "Failed to process line!" );
				RPRINT_VCFLINE;
				return false;
			}

		}///...for all samples

		//	clear any unused matrix rows
		//
		if( 1 == inf.end_line() )
		{

			//
			for( ; inf.cur_row < inf.nrow ; inf.cur_row ++ )
			{
				inf.ptr[inf.cur_row] = inf.empty_cell_value_int;
			}

			//
			inf.ptr += column_stepsize;

			//	Set Column Names = SNP Positions
			//
			if( R_NilValue != inf.colnamvec )
			{
				char posbuffer[256];
				snprintf(posbuffer,sizeof(posbuffer)-2,"%d",inf.snppos);
				SET_STRING_ELT( inf.colnamvec, inf.cur_col, mkChar(posbuffer) );
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
			//
			inf.cur_col--;
		}

	}///...for( each column in the matrix )
	
	//
	//printf("	inf.cur_col %d < inf.ncol %d\n",inf.cur_col,inf.ncol);
	bool completely_filled_matrix = !(inf.cur_col < inf.ncol);
	//printf("	T/F ? (%s)\n", completely_filled_matrix?"T":"F" );

	//
	//	Reset/clear unused columns
	//
	//TODO in loop_str, use _str for cell values
	for( unsigned int rcol = inf.cur_col; rcol < inf.ncol; rcol ++ )
	{
		//
		for( unsigned int rrow=0; rrow < inf.nrow ; rrow ++ )
		{
			inf.ptr[rrow] = inf.empty_cell_value_int;
		}
		inf.ptr += column_stepsize;
			
		//	Clear column name
		//
		if( R_NilValue != inf.colnamvec )
		{
			SET_STRING_ELT( inf.colnamvec, rcol, inf.empty_col_name );
		}
	}

	//
	return completely_filled_matrix;

}///...snpmat_run_loop_int_diploidonly








			//
			//
			//				INT,		UNFILTERED
			//
			//
			//








/*!	Run the main SNP matrix read loop
**	+ NO FILTER!
**	+ BIALLELIC SNPS ONLY
**	+ INTEGER MATRIX
**
**
**
*/
bool		snpmat_run_loop_int_nofilter( snpmat_read_info_int& inf )
{
//	DBG( "[%s]\n" , __FUNCTION_NAME__ );

	//
	char			*fieldptr=0;	
	unsigned int	column_stepsize = inf.nrow;

	//
	//
	//
	for( ; inf.cur_col < inf.ncol ; inf.cur_col ++ )
	{


		//-
		//	Get a valid biallelic SNP line from the file
		//-
		const char* refptr;
		const char* altptr;
		bool	bLineParsed;
		while( true == (bLineParsed = inf.f->parseNextLine()) )
		{
		
//			DBG("%d:%d [%s]\n",inf.cur_col,inf.ncol,inf.f->getFieldPtr( CHROM ));
			
			refptr = (char*)inf.f->getFieldPtr( REF );
			altptr = (char*)inf.f->getFieldPtr( ALT );

//			DBG("refptr %p altptr %p\n",refptr,altptr);
//			DBG("ref[%s]\n",refptr);

			//	- make sure its a bi-allelic SNP line
			//
			if( refptr && refptr[1] == '\t' )			//TODO use macro/inline
			{
//				DBG("ref_SNP\n");
				if( inf.biallelic_locus_required )
				{
					if( altptr && altptr[1] == '\t' )	//both REF and ALT have only single-char alleles?
					{
						break;	//..yes, process this line
					}
				}
				//
				//
				else
				{
					if( ALT_is_snp( altptr ) )		//ALT has any number of comma-separated single-char alleles ?
						break;	// ..yes, process this line
				}
			}///...if REF is single-char allele
			else
			{
//				DBG("not SNP ref (%d=%c)\n",refptr?refptr[1]:-1,refptr?refptr[1]:'?');
			}

		}///...while( could read another line from VCF )

		//
		//
		//	if there are no more lines in the VCF, break the loop early
		//
		//
		if( bLineParsed == false )
		{
			//df1("No more lines!\n");
//			DBG("nomorelines\n");
			break;
		}

		//
		//
		//	get POS as integer
		//
		//
		fieldptr = (char*)inf.f->getFieldPtr( POS );
		if( 0 == fieldptr )
		{
			RPRINT_ERROR( "Could not find POS column in line!" );
			RPRINT_VCFLINE;
			return false;
		}
		inf.snppos = atoi( fieldptr );
//		DBG("snppos=%d\n",inf.snppos);

		//-
		//
		//	for each selected individual, get the SNP genotype information
		//		and store a 0 or 1 in the biallelic matrix
		//
		//-
		
		inf.keep_column = 0;		//undecided whether to keep this column or not

		//
		inf.begin_line();
		
		//
		//
		//	for each sample, read genotype data, convert it to integer data and write into matrix
		//
		//
		for( inf.cur_row = 0; inf.cur_row < inf.f->num_wanted_samples ; inf.cur_row ++ )
		{
			
			//	- get field of sample
			//
			fieldptr = (char*)inf.f->getFieldPtr( inf.f->wanted_samples[inf.cur_row] );			
			if( 0 == fieldptr )
			{
				RPRINT_ERROR( "Could not access sample column %d\n" , inf.f->wanted_samples[inf.cur_row] );
				RPRINT_VCFLINE;
				return false;
			}
			
			//
			//
			//
			int process_gt_result = inf.process_sample( fieldptr );
			if( 0 == process_gt_result )
			{
				RPRINT_ERROR( "Failed to process line!" );
				RPRINT_VCFLINE;
				return false;
			}
			

		}///...for all samples
	
		//
		//	clear any unused matrix rows
		//
		if( 1 == inf.end_line() )
		{
			
			//
			for( ; inf.cur_row < inf.nrow ; inf.cur_row ++ )
			{
				inf.ptr[inf.cur_row] = inf.empty_cell_value_int;
			}
			
			//
			inf.ptr += column_stepsize;
			
			
			//	Set Column Names = SNP Positions
			//
			if( R_NilValue != inf.colnamvec )
			{
				char posbuffer[256];
				snprintf(posbuffer,sizeof(posbuffer)-2,"%d",inf.snppos);
				SET_STRING_ELT( inf.colnamvec, inf.cur_col, mkChar(posbuffer) );
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
			//
			inf.cur_col--;
		}

	}///...for( each column in the matrix )
	
	//
	bool completely_filled_matrix = !(inf.cur_col < inf.ncol);
	
	//
	//	Reset/clear unused columns
	//
	//TODO in loop_str, use _str for cell values
	for( unsigned int rcol = inf.cur_col; rcol < inf.ncol; rcol ++ )
	{
		//
		for( unsigned int rrow=0; rrow < inf.nrow ; rrow ++ )
		{
			inf.ptr[rrow] = inf.empty_cell_value_int;
		}
		inf.ptr += column_stepsize;
			
		//	Clear column name
		//
		if( R_NilValue != inf.colnamvec )
		{
			SET_STRING_ELT( inf.colnamvec, rcol, inf.empty_col_name );
		}
	}

	//
	return completely_filled_matrix;

}///...snpmat_run_loop_int_nofilter








/*!	Run the main SNP matrix read loop
**	+ NO FILTER!
**	+ BIALLELIC SNPS ONLY
**	+ INTEGER MATRIX
**	+ ONLY DIPLOID SAMPLE GENOTYPES
**
**
*/
bool		snpmat_run_loop_int_diploidonly_nofilter( snpmat_read_info_int& inf )
{
//	DBG( "[%s]\n" , __FUNCTION_NAME__ );

	//
	char			*fieldptr=0;	
	unsigned int	column_stepsize = inf.nrow;

	//
	//
	//
	for( ; inf.cur_col < inf.ncol ; inf.cur_col ++ )
	{
		//-
		//	Get a valid biallelic SNP line from the file
		//-
		const char* refptr;
		const char* altptr;
		bool	bLineParsed;
		while( true == (bLineParsed = inf.f->parseNextLine()) )
		{
			
			refptr = (char*)inf.f->getFieldPtr( REF );
			altptr = (char*)inf.f->getFieldPtr( ALT );

			//	- make sure its a bi-allelic SNP line
			//
			if( refptr && refptr[1] == '\t' )			//TODO use macro/inline
			{
				if( inf.biallelic_locus_required )
				{
					if( altptr && altptr[1] == '\t' )	//both REF and ALT have only single-char alleles?
					{
						break;	//..yes, process this line
					}
				}
				//
				//
				else
				{
					if( ALT_is_snp( altptr ) )		//ALT has any number of comma-separated single-char alleles ?
						break;	// ..yes, process this line
				}
			}///...if REF is single-char allele

		}///...while( could read another line from VCF )

		//
		//
		//	if there are no more lines in the VCF, break the loop early
		//
		//
		if( bLineParsed == false )
		{
			//df1("No more lines!\n");
			break;
		}

		//
		//
		//	get POS as integer
		//
		//
		fieldptr = (char*)inf.f->getFieldPtr( POS );
		if( 0 == fieldptr )
		{
			RPRINT_ERROR( "Could not find POS column in line!" );
			RPRINT_VCFLINE;
			return false;
		}
		inf.snppos = atoi( fieldptr );

		//-
		//
		//	for each selected individual, get the SNP genotype information
		//		and store a 0 or 1 in the biallelic matrix
		//
		//-
		
		inf.keep_column = 0;		//undecided whether to keep this column or not

		//
		inf.begin_line();
		
		//
		//
		//	for each sample, read genotype data, convert it to integer data and write into matrix
		//
		//
		for( inf.cur_row = 0; inf.cur_row < inf.f->num_wanted_samples ; inf.cur_row ++ )
		{

			//	- get field of sample
			//
			fieldptr = (char*)inf.f->getFieldPtr( inf.f->wanted_samples[inf.cur_row] );			
			if( 0 == fieldptr )
			{
				RPRINT_ERROR( "Could not access sample column %d\n" , inf.f->wanted_samples[inf.cur_row] );
				RPRINT_VCFLINE;
				return false;
			}

			//	check if sample genotype is diploid
			//
			if( !(is_GT_allele_divider(fieldptr[1])) || !(is_FORMAT_subfield_delimiter(fieldptr[3])) )
			{
				DBG_RPRINT_WARN( "Non-diploid data! (%s)\n" , fieldptr );
				DBG_RPRINT_VCFLINE;
				continue;
			}

			//
			//
			int process_gt_result = inf.process_sample( fieldptr );
			if( 0 == process_gt_result )
			{
				RPRINT_ERROR( "Failed to process line!" );
				RPRINT_VCFLINE;
				return false;
			}

		}///...for all samples

		//	clear any unused matrix rows
		//
		if( 1 == inf.end_line() )
		{

			//
			for( ; inf.cur_row < inf.nrow ; inf.cur_row ++ )
			{
				inf.ptr[inf.cur_row] = inf.empty_cell_value_int;
			}

			//
			inf.ptr += column_stepsize;

			//	Set Column Names = SNP Positions
			//
			if( R_NilValue != inf.colnamvec )
			{
				char posbuffer[256];
				snprintf(posbuffer,sizeof(posbuffer)-2,"%d",inf.snppos);
				SET_STRING_ELT( inf.colnamvec, inf.cur_col, mkChar(posbuffer) );
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
			//
			inf.cur_col--;
		}

	}///...for( each column in the matrix )
	
	//
	bool completely_filled_matrix = !(inf.cur_col < inf.ncol);

	//
	//	Reset/clear unused columns
	//
	//TODO in loop_str, use _str for cell values
	for( unsigned int rcol = inf.cur_col; rcol < inf.ncol; rcol ++ )
	{
		//
		for( unsigned int rrow=0; rrow < inf.nrow ; rrow ++ )
		{
			inf.ptr[rrow] = inf.empty_cell_value_int;
		}
		inf.ptr += column_stepsize;
			
		//	Clear column name
		//
		if( R_NilValue != inf.colnamvec )
		{
			SET_STRING_ELT( inf.colnamvec, rcol, inf.empty_col_name );
		}
	}

	//
	return completely_filled_matrix;

}///...snpmat_run_loop_int_diploidonly_nofilter










/*!
**
*/
SEXP _internal_VCF_snpmat_int( snpmat_read_info_int& inf, SEXP vcfptr, SEXP mat )
{
//	DBG( "[%s]\n" , __FUNCTION_NAME__ );
	//	initialise matrix read
	//
	if( false == snpmat_init_validate_int( inf , vcfptr, mat ) )
	{
		Rprintf(" init snpmat-loop failed!\n");
		return RBool::False();
	}

	//	run the main loop
	//
	bool runresult;
	if( inf.doFiltering )
	{
		if( inf.diploid_samples_only )
			runresult = snpmat_run_loop_int_diploidonly( inf );
		else
			runresult = snpmat_run_loop_int( inf );
	}
	else
	{
		if( inf.diploid_samples_only )
			runresult = snpmat_run_loop_int_diploidonly_nofilter( inf );
		else
			runresult = snpmat_run_loop_int_nofilter( inf );
	}

	//
	return (runresult ? RBool::True() : RBool::False() );

}//...



