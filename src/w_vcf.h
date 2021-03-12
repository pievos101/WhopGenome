/*
**
**		WhopGen
**
**		WHOle genome population genetics with popGEN
**
**
**

	TODO:

**
*/

//*
//*			INCLUDES
//*


	//

//*
//*			DEFINES
//*

//#define	RELAXED_STANDARD

	// any character that divides the chromosomes in a GT subfield ( a slash / for unphased, a bar | for phased, the obsolete backslash \ sign )
	//
//#ifdef	RELAXED_STANDARD
//#	define		is_GT_allele_divider( chr )		( ((chr) == '/') || ((chr) == '|') || ((chr) == '\\') )
//#else
//#	define		is_GT_allele_divider( chr )		( ((chr) == '/') || ((chr) == '|') )
#	define		is_GT_allele_divider( chr )		( ((chr) == '/') || ((chr) == '|') || ((chr) == '\\') )
//#endif

	//
#define		is_column_divider( ch )			( (ch) == '\t' )
#define		is_column_end( ch )				( (ch) == '\t' || (ch) == '\n' || (ch) == 0 )

	// any character that terminates a FORMAT/SAMPLE-subfield ( a colon : dividing subfields, a tab dividing columns, a newline or 0-byte ending lines )
#define		is_FORMAT_subfield_delimiter( ch )		 ( ((ch) == '\t') || ((ch) == '\n') || ((ch) == ':') || ((ch) == 0) )

	// any character that terminates a INFO-subfield ( a semicolon ; dividing subfields, a tab dividing columns, a newline or 0-byte ending lines )
#define		is_INFO_subfield_delimiter( ch )		 ( ((ch) == '\t') || ((ch) == '\n') || ((ch) == ';') || ((ch) == 0) )

	//
#define		is_alphanumeric( ch )			( ( ((ch) >= '0') && ((ch) <= '9') ) || (((ch) >= 'a') && ((ch) <= 'z')) || (((ch) >= 'A') && ((ch) <= 'Z')) )
#define		is_ascii_letter( ch )			( (((ch) >= 'a') && ((ch) <= 'z')) || (((ch) >= 'A') && ((ch) <= 'Z')) )

	// 01,2,3,4,5,6,7,8,9
#define		is_ascii_digit( ch )			( ((ch) >= '0') && ((ch) <= '9') )

	// 
#define		is_ascii_hexdigit( ch )			( ( ((ch) >= '0') && ((ch) <= '9') ) || (((ch) >= 'a') && ((ch) <= 'f')) || (((ch) >= 'A') && ((ch) <= 'F')) )
#define		ascii_digit_to_int( ch )		( (ch) - '0' )

	//	parameters for the filtering system : VCF column to inspect
	//
enum filterable_columns {
	CHROM		=	0,
	POS			=	1,
	ID			=	2,
	REF			=	3,
	ALT			=	4,
	QUAL		=	5,
	FILTER		=	6,
	INFO		=	7,
	FORMAT		=	8,
	NUM_FILTERABLE_COLUMNS
};

	//	parameters for the filtering system : how to compare data in VCF column with reference values of filter rule
	//
enum compare_how_enum {
	DOES_EXIST	=	0,	//subfield exists?
	INT_CMP,			//read int and compare
	INT_CMP_OO,
	INT_CMP_OC,
	INT_CMP_CO,
	INT_CMP_CC,
	FLT_CMP,			//
	FLT_CMP_OO,
	FLT_CMP_OC,
	FLT_CMP_CO,
	FLT_CMP_CC,
	NUM_CMP_OPS
};

	//	parameters for the filtering system : how to treat a rule match in terms of processing the line further
	//
enum do_what_enum {
	DW_NOP		=	0x00,
	
	DW_SKIP		=	0x01,			//skip if false, dont test further
	DW_KEEP		=	0x02,			//keep if true, dont test further
	DW_SKIP_NOT	=	0x03,			//skip if true, dont test further
	DW_KEEP_NOT	=	0x04,			//keep if false, dont test further
	
	DW_NUM_NONDISABLED_ACTIONS,
	
		//	rules can be selectively disabled for testing purposes .. numeric codes have bit 7 set (OR with 0x80)
		//
	DW_NOP_DISABLED		=	0x80,	//actually a non-sense meaning
	DW_SKIP_DISABLED	=	0x81,	//DISABLED : skip if false, dont test further
	DW_KEEP_DISABLED	=	0x82,	//DISABLED : keep if true, dont test further
	DW_SKIP_NOT_DISABLED=	0x83,	//DISABLED : skip if true, dont test further
	DW_KEEP_NOT_DISABLED=	0x84,	//DISABLED : keep if false, dont test further

	NUM_DO_WHATS
};

#define		MAX_NUM_RULES			20
#define		MAX_NUM_SUBFIELDS		3

//*
//*			STRUCTS
//*


//	filter rules
//
struct filter_entry_t
{
	int		filtered_column;		// 0-8 for CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
	char	filtered_subfield[32];	// name of subfield in column (only valid for INFO or FORMAT)
	int		compare_how;			// one of some constants describing comparison operation
	int		do_what;				// how to consider the line for further processing (skip, keep, ...)
		// reference values for comparison
	union {
		struct {
			int		i1,i2;
		} integers;
		struct {
			float	f1,f2;
		} floats;
	} ref;
};

//*
//*			CLASSES
//*



/*!	Class wrapping working with VCF files
**
**
**
*/
class vcff : public whop_tabix
{
public:
	
		//-
		//-		CTORs + DTOR
		//-

	vcff();
	vcff(const char * filename );
	~vcff();
	
		//-
		//-		METHODS
		//-
	
	//	file
	//
	bool			open( const char * filename );


	//	fields
	//
	unsigned int	getNumFields( void );
	const char*	getFieldName( unsigned int idx );
	
	//	samples
	//
	unsigned int	getNumSamples( void );
	unsigned int	getFirstSampleFieldIndex( void );
	unsigned int	getSampleIndexByName( const char * samplename );
	
	//!
	void			resetSampleSelection( void )
	{
		num_wanted_samples=0;
	}

	//!
	bool			selectSample( unsigned int samplefieldidx )
	{
		
		//	sanity check : sample 
		if( samplefieldidx >= num_fields ) {
			df0("(ERROR) selectSample : sample field index %d > total number of fields %d\n",samplefieldidx,num_fields);
			return false;
		}
		
		//	sanity check : number of wanted samples should not exceed number of available samples
		if( num_wanted_samples >= num_samples ) {
			df0("(ERROR) selectSample : number of wanted samples %d >= total number of available sampes %d\n",num_wanted_samples,num_samples);
			return false;
		}

		//	disallow redundant selection of a sample
		for( unsigned int i=0; i < num_wanted_samples; i ++ )
		{
			if( wanted_samples[i] == samplefieldidx ) {
				df1("(WARNING) selectSample : sample field %d already listed as wanted!\n");
				return false;
			}
		}
				
		//
		wanted_samples[ num_wanted_samples ] = samplefieldidx;
		wanted_samples[ num_wanted_samples+1 ] = -1;
		num_wanted_samples++;
		return true;
	}
	
	//	reading
	//
	unsigned int	lastReadPos( void ){
		return _last_readline_pos;
	}

	//	debugging
	//
	bool			testfunc( void );

	
		//-
		//-		DATA
		//-
	
	//	Filter rules
	//
	int							num_rules_used;
	int							num_fieldnames_used;
	struct filter_entry_t		ruleset[ MAX_NUM_RULES ];
	
	//	Parse-related (transient) data
	//

		//	samples to extract
	
	unsigned int				*wanted_samples;	//
	unsigned int				num_wanted_samples;


	//-
	//-		DATA
	//-

protected:
	
	//	VCF specific data
	//
	unsigned long				num_fields;			//number of fields in VCF/TSV file
	unsigned long				num_samples;		//number of samples in VCF
	std::vector<std::string>	field_names;		//CHROM,POS,ID,REF,ALT,... from the "#CHROM	POS	ID	REF	ALT..." line
	unsigned long				sample_begin_index;
	unsigned long				_last_readline_pos;

private:

	//
};




/*!
**
*/
class snpmat_read_info_int {
public:
	snpmat_read_info_int( const char * callfuncname ) : 
		caller_name(callfuncname), errormessage(0),
		nrow(0), ncol(0), ptr(0), colnamvec(R_NilValue),
		empty_col_name(UNUSED_COLUMN_VALUE), empty_cell_value_int(-2), error_cell_value_int(-1),
		f(0),
		drop_uninformative_columns(true),keep_column(false),
		cur_col(0), cur_row(0),
		biallelic_locus_required(0), diploid_samples_only(0), doFiltering(0)
	{}
	
	virtual ~snpmat_read_info_int(){}
	
	//
	//
	//
	virtual		int		begin_line( void ) = 0;
	virtual		int		end_line( void ) = 0;
	virtual		int		process_sample( const char* gt ) = 0;
	
	//	Additional information for the user
	//
	const char			*caller_name;
	char				*errormessage;
	
	//	Matrix
	//
	unsigned int		nrow;
	unsigned int		ncol;
	
	int					*ptr;
	
	SEXP				colnamvec;
	
	SEXP				empty_col_name;			//name to mark empty columns with
	int					empty_cell_value_int;	//value to mark empty cells of INTEGER matrices
	int					error_cell_value_int;	//value to mark empty cells of INTEGER matrices
	
	//	VCF
	//
	class vcff			*f;

	//
	//
	bool				drop_uninformative_columns;
	bool				keep_column;					// 0 = dont know/column just started , 1 = yes definitely, -1 = definitely NOT
	
	//	current loop variables
	//
	unsigned int		snppos;
	unsigned int		cur_col;
	unsigned int		cur_row;
	
	//	Loop conditional processing
	//
	bool				biallelic_locus_required;		//	only consider biallelic SNP loci (1 REF 1 ALT) instead of (1 REF x ALT)
	bool				diploid_samples_only;			//	only keep columns in the matrix showing actual variation (may be sample selection dependent!)

	bool				doFiltering;					//	run line by the filter-rules

};





/*!
**
*/
class snpmat_read_info_str {
public:
	snpmat_read_info_str( const char * callfuncname ) : 
		caller_name(callfuncname), errormessage(0),
		nrow(0), ncol(0), sptr(0), colnamvec(R_NilValue),
		empty_col_name(UNUSED_COLUMN_VALUE), empty_cell_value_str(UNUSED_CELL_VALUE), error_cell_value_str(ERROR_CELL_VALUE),
		f(0),
		drop_uninformative_columns(true),keep_column(false),
		cur_col(0), cur_row(0),
		biallelic_locus_required(0), diploid_samples_only(0), doFiltering(0)
	{}
	
	virtual ~snpmat_read_info_str(){}
	
	//
	//
	//
	virtual		int		begin_line( void ) = 0;
	virtual		int		end_line( void ) = 0;
	virtual		int		process_sample( const char* gt ) = 0;
	
	//	Additional information for the user
	//
	const char			*caller_name;
	char				*errormessage;
	
	//	Matrix
	//
	unsigned int		nrow;
	unsigned int		ncol;
	
	SEXP				*sptr;
	
	SEXP				colnamvec;
	
	SEXP				empty_col_name;			//name to mark empty columns with
	SEXP				empty_cell_value_str;	//value to mark empty cells of STRSXP matrices
	SEXP				error_cell_value_str;
	
	//	VCF
	//
	class vcff			*f;

	//
	//
	bool				drop_uninformative_columns;
	bool				keep_column;					// 0 = dont know/column just started , 1 = yes definitely, -1 = definitely NOT
	
	//	current loop variables
	//
	unsigned int		snppos;
	unsigned int		cur_col;
	unsigned int		cur_row;
	
	//	Loop conditional processing
	//
	bool				biallelic_locus_required;		//	only consider biallelic SNP loci (1 REF 1 ALT) instead of (1 REF x ALT)
	bool				diploid_samples_only;			//	only keep columns in the matrix showing actual variation (may be sample selection dependent!)

	bool				doFiltering;					//	run line by the filter-rules

};


//*
//*			PROTOTYPES
//*

	//
	//

EXPORT SEXP VCF_snpmat_diplo_bial_hasalt_filtered( SEXP vcfptr, SEXP mat );
EXPORT SEXP VCF_snpmat_diplo_bial_hasalt_unfiltered( SEXP vcfptr, SEXP mat );
EXPORT SEXP VCF_snpmat_diplo_anyal_hasalt_filtered( SEXP vcfptr, SEXP mat );
EXPORT SEXP VCF_snpmat_diplo_anyal_hasalt_unfiltered( SEXP vcfptr, SEXP mat );

EXPORT SEXP VCF_snpmat_diplo_bial_nucodes_filtered( SEXP vcfptr, SEXP mat );
EXPORT SEXP VCF_snpmat_diplo_bial_nucodes_unfiltered( SEXP vcfptr, SEXP mat );
EXPORT SEXP VCF_snpmat_diplo_anyal_nucodes_filtered( SEXP vcfptr, SEXP mat );
EXPORT SEXP VCF_snpmat_diplo_anyal_nucodes_unfiltered( SEXP vcfptr, SEXP mat );

EXPORT SEXP VCF_snpmat_diplo_bial_ishet_filtered( SEXP vcfptr, SEXP mat );
EXPORT SEXP VCF_snpmat_diplo_anyal_ishet_filtered( SEXP vcfptr, SEXP mat );
EXPORT SEXP VCF_snpmat_diplo_bial_ishet_unfiltered( SEXP vcfptr, SEXP mat );
EXPORT SEXP VCF_snpmat_diplo_anyal_ishet_unfiltered( SEXP vcfptr, SEXP mat );


