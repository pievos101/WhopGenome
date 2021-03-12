/*
**
**		WhopGen
**
**		WHOle genome population genetics with popGEN
**
**
**		VCF Filtering functions
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

	//	R includes
	//
#include	<R.h>
#include	<Rinternals.h>
#include	<R_ext/Rdynload.h>



	//

using namespace std;


//*
//*			DEFINES
//*

#define	filterLineDF		if( false == true ) Rprintf

//*
//*			STRUCTS
//*


//*
//*			CLASSES
//*

//*
//*			DATA
//*




//-----

/*

//! Short descriptions of the comparison operators, matching the ENUMs
static const char * cmptype_to_string[] = {
	"EXISTS",
	"INT==",
	"INT(,)",
	"INT(,]",
	"INT[,)",
	"INT[,]",
	"FLT==",
	"FLT()",
	"FLT(]",
	"FLT[)",
	"FLT[]",
	0
};
*/

//-----
//! Names of VCF column
static const char * column_names[] = {
	"CHROM",
	"POS",
	"ID",
	"REF",
	"ALT",
	"QUAL",
	"FILTER",
	"INFO",
	"FORMAT",
	0
};

//*
//*			EXTERNS
//*



//*
//*			CODE
//*



/*!
**
*/
bool	get_subfield( vcff* f, unsigned int rulenum, char**ptr )
{
	//TODO validate f, rulenum, ptr
	if( 0 == f || MAX_NUM_RULES < rulenum || 0 == ptr )
	{
		Rprintf("(!!) get_subfield : error : invalid arguments (f=%x, rulenum=%d/%x, ptr=%x)\n",f,rulenum,rulenum,ptr);
		return false;
	}
	
	//
	*ptr=0;

	// get pointer to start of string for column number <filtered_column> in the table
	//
	const char * columnptr = f->getFieldPtr( f->ruleset[ rulenum ].filtered_column );
	if( 0 == columnptr )
	{
		Rprintf("(!!) get_subfield : error : stringptr to field %d is at %x\n", f->ruleset[ rulenum ].filtered_column ,columnptr);
		return false;
	}

	//	get pointer to sub-field inside column to look for (INFO: AP,RD,SQ etc.)
	//		can be empty if the entire column's text is to be considered
	//
	const char *fnam = f->ruleset[ rulenum ].filtered_subfield;
	
	//
	//
	if( *fnam == 0 )
	{
		*ptr = (char*)columnptr;
		return true;
	}

	//
	//Rprintf("fieldnam='%s'\n",fnam);

	//
	//
	int cnt=500;
	while( cnt-- )
	{
		int i=0;
		//FIXME : R-code vcf_addfilter() uses toupper() on submitted fieldnam, this should probably better be handled in the C function
		for( ; fnam[i] == toupper(*columnptr) ; i++ )
		{
			columnptr++;
		}
//		Rprintf("last cmp was : %c <-> %c\n",fnam[i],*info);

		//
		//
		if( fnam[i] == 0 )
		{
			//<fnam[i] == 0 : has fieldname ended too ?>
			if( *columnptr == '=' )
			{
				//
				*ptr = (char*)columnptr+1;

				//success, get ptr behind =
				//
			//	Rprintf("success:[[");
			//	for( int i=0; i < 40 && a.data.fldptr[i] != '\t'; i++ )
			//		Rprintf("%c",a.data.fldptr[i] );
			//	Rprintf("]]\n");

				//
				return true;
			}
			else
			if( *columnptr == ';' || *columnptr == '\t' || *columnptr == 0 )	//fieldname is last entry and has no params ?( e.g. AC=200;FLAG_XYZ\t )
			{
				//success, set ptr to 0
			//	Rprintf("success: no data associated[[");
			//	for( int i=0; i < 40 && info[i] != '\t'; i++ )
			//		Rprintf("%c",columnptr[i] );
			//	Rprintf("]]\n");
				*ptr = 0;
				return true;
			}
			//else
			//{
			//	// fnam at end, but fieldname not yet
			//}
		}

		//fail, find next ';'
		//
//		Rprintf("	subfield didnt match, find next:(");
		for( ; *columnptr != ';' ; columnptr++ )
		{
//			Rprintf("%c",*columnptr);
			if( *columnptr == '\t' || *columnptr == 0 )
				return false;
		}
		columnptr++;
//		Rprintf("\n");

		//
	}

	//
	return false;
}

//---------------------------------------------

/*!	Applies the ruleset stored in the given VCFhandle and returns TRUE if this line passes
**
**
**
*/
bool	filterLine( vcff* f )
{
	int		numfes = f->num_rules_used;
	char	*subfieldptr;
	int		iv;
	float	fv;
//	filterLineDF("[%d RULES]\n",numfes);


	//
	//
	for( int i=0; i < numfes; i++ )			//for every configured rule...
	{
		//
		//
		const filter_entry_t & fe = f->ruleset[i];
//		filterLineDF("Rule[%d]\n",i);
		
		//	if action = DISABLED kind, just continue
		//
		if( fe.do_what & 0x80 )
		{
			filterLineDF("	disabled action[%02x]\n",fe.do_what);
			continue;
		}
		
		if( fe.do_what >= DW_NUM_NONDISABLED_ACTIONS )
		{
			Rprintf("(!!) Error : invalid action %d > maximum %d!\n",fe.do_what,DW_NUM_NONDISABLED_ACTIONS);
			continue;
		}

		// get ptr to data of field
		//
		subfieldptr=0;
		bool cmpresult = get_subfield( f, i , &subfieldptr );
		if( (f->ruleset[i].compare_how != DOES_EXIST) && (cmpresult==false) )
		{
			df1("	could not get field!\n");
			//TODO : what to do if a subfield of INFO could not be found ?
			//	depends on f->fieldset : test for MUSTEXIST bit
			//
			continue;
		}
		filterLineDF("	field=[%s]\n",subfieldptr);

		//
		//	got ptr to subfield and do the comparison
		//
		switch( f->ruleset[i].compare_how )
		{
			//	directly related to result of get_subfield() above
			//
//			case DOES_EXIST:
//				cmpresult = true;
//				break;
			//	integer comparisons
			//
			case INT_CMP:
				iv = atoi( subfieldptr );
				cmpresult = (f->ruleset[i].ref.integers.i1 == iv );
				break;
			//
			//
			case INT_CMP_OO:
				iv = atoi( subfieldptr );
				cmpresult = ((f->ruleset[i].ref.integers.i1 < iv ) && (iv < f->ruleset[i].ref.integers.i2));
				break;
			//
			//
			case INT_CMP_OC:
				iv = atoi( subfieldptr );
				cmpresult = ((f->ruleset[i].ref.integers.i1 < iv ) && (iv <= f->ruleset[i].ref.integers.i2));
				break;
			//
			//
			case INT_CMP_CO:
				iv = atoi( subfieldptr );
				cmpresult = ((f->ruleset[i].ref.integers.i1 <= iv ) && (iv < f->ruleset[i].ref.integers.i2));
				break;
			//
			//
			case INT_CMP_CC:
				iv = atoi( subfieldptr );
				cmpresult = ((f->ruleset[i].ref.integers.i1 <= iv ) && (iv <= f->ruleset[i].ref.integers.i2));
				break;
			
			
			//
			//
			//	float comparisons
			//
			//


			//
			//
			case FLT_CMP:
				fv = atof( subfieldptr );
				cmpresult = (f->ruleset[i].ref.floats.f1 == fv );
				break;
			//
			//
			case FLT_CMP_OO:
				fv = atof( subfieldptr );
				cmpresult = ((f->ruleset[i].ref.floats.f1 < fv ) && (fv < f->ruleset[i].ref.floats.f2));
				break;
			//
			//
			case FLT_CMP_OC:
				fv = atof( subfieldptr );
				cmpresult = ((f->ruleset[i].ref.floats.f1 < fv ) && (fv <= f->ruleset[i].ref.floats.f2));
				break;
			//
			//
			case FLT_CMP_CO:
				fv = atof( subfieldptr );
				cmpresult = ((f->ruleset[i].ref.floats.f1 <= fv ) && (fv < f->ruleset[i].ref.floats.f2));
				break;
			//
			//
			case FLT_CMP_CC:
				fv = atof( subfieldptr );
				cmpresult = ((f->ruleset[i].ref.floats.f1 <= fv ) && (fv <= f->ruleset[i].ref.floats.f2));
				break;
				
			
			//
			//
			//
			//
			default:
				df1("	compare-how : unknown op %d!",f->ruleset[i].compare_how);
				break;
		}//switch( compare_how )

		filterLineDF("	cmpresult=%s\n",cmpresult?"TRUE":"FALSE");

		//
		//	react to the result of the comparison
		//
		switch( f->ruleset[i].do_what )
		{

			//	skip line if this filter matches
			//
			case DW_SKIP_NOT://	=	0x81,			//skip if true, dont test further
				filterLineDF("	action=skipnot..");
				cmpresult=!cmpresult;
			case DW_SKIP://		=	0x01,			//skip if false, dont test further
				filterLineDF("	action=skip\n");
				if( cmpresult )
				{
					filterLineDF("	DROPPING line!\n");
					return false;
				}
				break;

			//	keep line if this filter matches
			//
			case DW_KEEP_NOT://	=	0x82,			//keep if false, dont test further
				filterLineDF("	action=keepnot..");
				cmpresult=!cmpresult;
			case DW_KEEP://		=	0x02,			//keep if true, dont test further
				filterLineDF("	action=keep\n");
				if( cmpresult )
				{
					filterLineDF("	LINE ACCEPTED\n");
					return true;
				}
				break;
			//
			case DW_NOP:
				filterLineDF("	NOP\n");
				break;
/*
			case DW_NOP_DISABLED:
			case DW_SKIP_NOT_DISABLED:
			case DW_SKIP_DISABLED:
			case DW_KEEP_NOT_DISABLED:
			case DW_KEEP_DISABLED:
				//rule is disabled
				break;
*/
			default:
				filterLineDF("UNKNOWN ACTION %d!\n", f->ruleset[i].do_what );
				break;
		}

		//
	}//...for all rules

//	filterLineDF("	LINE ACCEPTED\n");
	return true;
}








/*!	Remove all filters
**
**
**
*/
EXPORT SEXP VCF_clearFilters( SEXP vcfptr )
{

	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("Parameter not a VCFhandle EXTPTR!\n");
		return R_NilValue;
	}

	//
	f->num_rules_used = 0;
	f->num_fieldnames_used = 0;

	//
	for( int i=0; i < MAX_NUM_RULES; i++ )
	{
		f->ruleset[i].filtered_column = -1;
		memset( & f->ruleset[i].filtered_subfield[0] , 0, sizeof(f->ruleset[i].filtered_subfield) );
		f->ruleset[i].compare_how = 0;
		f->ruleset[i].do_what = 0;
		f->ruleset[i].ref.integers.i1 =
		f->ruleset[i].ref.integers.i2 = 0;
	}

	//
	return R_NilValue;
}


/*!	Add a filter
**
**
*/
EXPORT SEXP VCF_addFilter( SEXP vcfptr, SEXP columnidx_sexp, SEXP fieldnam, SEXP cmptype, SEXP action, SEXP arg1, SEXP arg2 )
{

	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("(!!) Error : Parameter 1 not a VCFhandle EXTPTR!\n");
		return RBool::False();
	}

	//
	//
	if( f->num_rules_used >= MAX_NUM_RULES )
	{
		Rprintf("(!!) Error : already reached maximum number of %d rules!\n",MAX_NUM_RULES);
		return RBool::False();
	}
	
	//
	//
	int columnidx = INTEGER( columnidx_sexp )[0];

	//
	//
	const char * fieldname = RString::get( fieldnam , 0 );
	if( fieldname == 0 && ( (columnidx==7) || (columnidx==8) ) )
	{
		Rprintf("(!!) Error : Parameter 2, fieldname, is empty but required for columns INFO (7) and FORMAT (8)!\n");
		return RBool::False();
	}

	//	check comparison type
	//
	int cmp = INTEGER( cmptype )[0];
	if( cmp < DOES_EXIST || cmp > NUM_CMP_OPS )
	{
		Rprintf("(!!) Error : Parameter 3, cmptype, value %d not within range [0,%d]!\n",cmp,NUM_CMP_OPS);
		return RBool::False();
	}

	//	check action type
	//
	unsigned int act = INTEGER( action )[0];
	if( (act&0x80) >= DW_NUM_NONDISABLED_ACTIONS  )
	{
		Rprintf("(!!) Error : Parameter 4, action (value=%d) is not valid!\n",act);
		return RBool::False();
	}

	//
	//
	filter_entry_t * fe = &f->ruleset[f->num_rules_used];
	
	//
	fe->filtered_column = columnidx;

	//
	fe->compare_how = cmp;			//how to compare

	//
	//
	if( cmp >= INT_CMP && cmp <= INT_CMP_CC )	//1...5
	{
		fe->ref.integers.i1 = INTEGER( arg1 )[0];
		fe->ref.integers.i2 = INTEGER( arg2 )[0];
	}
	else if( cmp >= FLT_CMP && cmp <= FLT_CMP_CC )	//6...10
	{
		fe->ref.floats.f1 = REAL( arg1 )[0];
		fe->ref.floats.f2 = REAL( arg2 )[0];
	}

	strncpy( &fe->filtered_subfield[0] , fieldname , sizeof(fe->filtered_subfield) );

	fe->do_what = act;

	f->num_rules_used++;

	//
	return RBool::True();
}




/*!	Print a human-readable description of the current filter-rule set for the given VCF file
**
*/
EXPORT SEXP VCF_describeFilterConfig( SEXP vcfptr )
{

	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("Parameter not a VCFhandle EXTPTR!\n");
		return RBool::False();
	}

	//
	//
	Rprintf(	"Filtering rules:\n"
				"Rule#	If in column	in subfield	check that	and if so,\n"
				"---------\n"
				);
	
	//
	//
	for( int i=0; i < f->num_rules_used; i ++ )
	{
		
		const filter_entry_t & fe = f->ruleset[i];
		
		//
		Rprintf("#%d\t",i);
		
		//
		if( fe.filtered_column >= CHROM && fe.filtered_column <= FORMAT )
		{
			Rprintf("c=%2d n=(%-10s)\t", fe.filtered_column, column_names[ fe.filtered_column ] );
			Rprintf("sf=%-11s\t", &fe.filtered_subfield[0] );
		}
		else
		{
			Rprintf("c=%-10d\t", fe.filtered_column );
			Rprintf("sf=%-11s\t", &fe.filtered_subfield[0] );
		}
		
		//	describe the type of comparison to do and against which reference values
		//
		switch( fe.compare_how )
		{
			case DOES_EXIST:
				Rprintf("exists\t");
				break;
			case INT_CMP:
				Rprintf("== integer %d\t",fe.ref.integers.i1);
				break;
			case INT_CMP_OO:
				Rprintf("in range (%d,%d)\t",fe.ref.integers.i1,fe.ref.integers.i2);
				break;
			case INT_CMP_OC:
				Rprintf("in range (%d,%d]\t",fe.ref.integers.i1,fe.ref.integers.i2);
				break;
			case INT_CMP_CO:
				Rprintf("in range [%d,%d)\t",fe.ref.integers.i1,fe.ref.integers.i2);
				break;
			case INT_CMP_CC:
				Rprintf("in range [%d,%d]\t",fe.ref.integers.i1,fe.ref.integers.i2);
				break;
			//
			//
			case FLT_CMP:
				Rprintf("== float %f\t",fe.ref.floats.f1);
				break;
			case FLT_CMP_OO:
				Rprintf("in range (%f,%f)\t",fe.ref.floats.f1,fe.ref.floats.f2);
				break;
			case FLT_CMP_OC:
				Rprintf("in range (%f,%f]\t",fe.ref.floats.f1,fe.ref.floats.f2);
				break;
			case FLT_CMP_CO:
				Rprintf("in range [%f,%f)\t",fe.ref.floats.f1,fe.ref.floats.f2);
				break;
			case FLT_CMP_CC:
				Rprintf("in range [%f,%f]\t",fe.ref.floats.f1,fe.ref.floats.f2);
				break;
			default:
				Rprintf("UNKNOWN COMPARISON MODE %d!\t",fe.compare_how);
		}//..switch ( compare_how )


		//	describe what to do, when the comparison is true
		//
		switch( fe.do_what )
		{
			case DW_SKIP_NOT:
				Rprintf("DO NOT SKIP (active)");
				break;
			case DW_SKIP:
				Rprintf("SKIP (active)");
				break;
			case DW_KEEP_NOT:
				Rprintf("DO NOT KEEP (active)");
				break;
			case DW_KEEP:
				Rprintf("KEEP (active)");
				break;
			//
			case DW_SKIP_NOT_DISABLED:
				Rprintf("DO NOT SKIP (disabled)");
				break;
			case DW_SKIP_DISABLED:
				Rprintf("SKIP (disabled)");
				break;
			case DW_KEEP_NOT_DISABLED:
				Rprintf("DO NOT KEEP (disabled)");
				break;
			case DW_KEEP_DISABLED:
				Rprintf("then KEEP (disabled)");
				break;
			//
			case DW_NOP:
				Rprintf("NOP");
				break;
			default:
				Rprintf("UNKNOWN:%d!",fe.do_what);
				break;
		}

		//
		Rprintf("\n");
	}
	Rprintf("---------\n");

	//
	return R_NilValue;
}





/*!
**
**
*/
EXPORT	SEXP	VCF_isSNP( SEXP vcfptr )
{
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f ){	Rprintf("Parameter not a VCFhandle EXTPTR!\n");	return R_NilValue;	}
	if( 0 == f->getFieldPtr( REF ) ){	Rprintf("No previously parsed line available!\n");	return R_NilValue;	}
	if( 0 == f->getFieldPtr( ALT ) ){	Rprintf("No previously parsed line available!\n");	return R_NilValue;	}
	if( _internal_isSNP( f->getFieldPtr( REF ), f->getFieldPtr( ALT ) ))
		return RBool::True();
	return RBool::False();
}




/*!
**
**
*/
EXPORT	SEXP	VCF_isInDel( SEXP vcfptr )
{
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f ){	Rprintf("Parameter not a VCFhandle EXTPTR!\n");	return R_NilValue;	}
	if( 0 == f->getFieldPtr( REF ) ){	Rprintf("No previously parsed line available!\n");	return R_NilValue;	}
	if( 0 == f->getFieldPtr( ALT ) ){	Rprintf("No previously parsed line available!\n");	return R_NilValue;	}
	if( ! _internal_isSNP( f->getFieldPtr( REF ), f->getFieldPtr( ALT ) ))
		return RBool::True();
	return RBool::False();
}





		//*
		//*			Filtering rule modifications
		//*
		//*




/*!	Set a filtering rule's action
**
**
*/
EXPORT SEXP VCF_setRuleAction( SEXP vcfptr, SEXP ruleidx_sexp, SEXP action )
{

	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("(!!) Error : Parameter 1 not a VCFhandle EXTPTR!\n");
		return RBool::False();
	}
	
	//	if rule index is not an integer, return false or handle it somehow
	//
	if( RNumeric::isInt( ruleidx_sexp ) == false )
	{
		Rprintf("Expecting an integer as parameter 2 <ruleidx> in VCF_setRuleAction!\n");
		return RBool::False();
	}

	//	make sure the rule# is valid ( >= 0 and <= currently set number of rules)
	//
	int ruleidx = INTEGER( ruleidx_sexp )[0];
	if( f->num_rules_used < ruleidx || 0 > ruleidx )
	{
		Rprintf("(!!) Error : VCF_setRuleAction : invalid rule-number to change\n");
		return RBool::False();
	}
	
	//	get pointer to wanted rule-definition
	//
	filter_entry_t * fe = &f->ruleset[ruleidx];
	
	//	convert and validate parameter to modify
	//
	int newaction = INTEGER(action)[0];
	switch( newaction )
	{
		case DW_SKIP_NOT:
		case DW_SKIP:
		case DW_KEEP_NOT:
		case DW_KEEP:
		case DW_SKIP_NOT_DISABLED:
		case DW_SKIP_DISABLED:
		case DW_KEEP_NOT_DISABLED:
		case DW_KEEP_DISABLED:
		case DW_NOP:
			fe->do_what = newaction;
			break;
		default:
			Rprintf("(!!) Error : unknown action %d!",newaction);
			return RBool::False();
	}
	
	return( RBool::True() );
}





/*!	Enable a filtering rule
**
**
*/
EXPORT SEXP VCF_setRuleEnabled( SEXP vcfptr, SEXP ruleidx_sexp )
{

	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("(!!) Error : Parameter 1 not a VCFhandle EXTPTR!\n");
		return RBool::False();
	}
	
	//	if rule index is not an integer, return false or handle it somehow
	//
	if( RNumeric::isInt( ruleidx_sexp ) == false )
	{
		Rprintf("Expecting an integer as parameter 2 <ruleidx> in VCF_setRuleEnabled!\n");
		return RBool::False();
	}

	//	make sure the rule# is valid ( >= 0 and <= currently set number of rules)
	//
	int ruleidx = INTEGER( ruleidx_sexp )[0];
	if( f->num_rules_used < ruleidx || 0 > ruleidx )
	{
		Rprintf("(!!) Error : VCF_setRuleAction : invalid rule-number to change\n");
		return RBool::False();
	}
	
	//	get pointer to wanted rule-definition
	//
	filter_entry_t * fe = &f->ruleset[ruleidx];
	
	//	convert and validate parameter to modify
	//
	switch( fe->do_what )
	{
		case DW_SKIP_NOT_DISABLED:
			fe->do_what = DW_SKIP_NOT;
			break;			
		case DW_SKIP_DISABLED:
			fe->do_what = DW_SKIP;
			break;
		case DW_KEEP_NOT_DISABLED:
			fe->do_what = DW_KEEP_NOT;
			break;
		case DW_KEEP_DISABLED:
			fe->do_what = DW_KEEP;
			break;
	}
	
	return( RBool::True() );
}






/*!	Disable a filtering rule
**
**
*/
EXPORT SEXP VCF_setRuleDisabled( SEXP vcfptr, SEXP ruleidx_sexp )
{

	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("(!!) Error : Parameter 1 not a VCFhandle EXTPTR!\n");
		return RBool::False();
	}
	
	//	if rule index is not an integer, return false or handle it somehow
	//
	if( RNumeric::isInt( ruleidx_sexp ) == false )
	{
		Rprintf("Expecting an integer as parameter 2 <ruleidx> in VCF_setRuleDisabled!\n");
		return RBool::False();
	}

	//	make sure the rule# is valid ( >= 0 and <= currently set number of rules)
	//
	int ruleidx = INTEGER( ruleidx_sexp )[0];
	if( f->num_rules_used < ruleidx || 0 > ruleidx )
	{
		Rprintf("(!!) Error : VCF_setRuleAction : invalid rule-number to change\n");
		return RBool::False();
	}
	
	//	get pointer to wanted rule-definition
	//
	filter_entry_t * fe = &f->ruleset[ruleidx];
	
	//	change action to the disabled state
	//
	switch( fe->do_what )
	{
		case DW_SKIP_NOT:
			fe->do_what = DW_SKIP_NOT_DISABLED;
			break;			
		case DW_SKIP:
			fe->do_what = DW_SKIP_DISABLED;
			break;
		case DW_KEEP_NOT:
			fe->do_what = DW_KEEP_NOT_DISABLED;
			break;
		case DW_KEEP:
			fe->do_what = DW_KEEP_DISABLED;
			break;
	}
	
	return( RBool::True() );
}






/*!	For a given filtering rule, set the VCF column it should check
**
**
*/
EXPORT SEXP VCF_setRuleColumn( SEXP vcfptr, SEXP ruleidx_sexp, SEXP columnidx_sexp )
{

	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("(!!) Error : Parameter 1 not a VCFhandle EXTPTR!\n");
		return RBool::False();
	}
	
	//	if rule index is not an integer, return false or handle it somehow
	//
	if( RNumeric::isInt( ruleidx_sexp ) == false )
	{
		Rprintf("Expecting an integer as parameter 2 <ruleidx> in VCF_setRuleColumn!\n");
		return RBool::False();
	}
	
	//	if column index is not an integer, return false or handle it somehow
	//
	if( RNumeric::isInt( columnidx_sexp ) == false )
	{
		Rprintf("Expecting an integer as parameter 3 <column> in VCF_setRuleColumn!\n");
		return RBool::False();
	}

	//	make sure the rule# is valid ( >= 0 and <= currently set number of rules)
	//
	int ruleidx = INTEGER( ruleidx_sexp )[0];
	if( f->num_rules_used < ruleidx || 0 > ruleidx )
	{
		Rprintf("(!!) Error : VCF_setRuleColumn : invalid rule-number to change\n");
		return RBool::False();
	}
	
	//	get pointer to wanted rule-definition
	//
	filter_entry_t * fe = &f->ruleset[ruleidx];
	
	//	convert and validate parameter to modify
	//
	int columnidx = INTEGER(columnidx_sexp)[0];
	if( columnidx >= 0 && columnidx <= 8 )
	{
		fe->filtered_column = columnidx;
		return( RBool::True() );
	}

	return RBool::False();
}







/*!	Set a filtering rule's reference values
**	- infer type (float,int) from type of comparison
**
*/
EXPORT SEXP VCF_setRuleRefValues( SEXP vcfptr, SEXP ruleidx_sexp, SEXP ref1, SEXP ref2 )
{

	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("(!!) Error : Parameter 1 not a VCFhandle EXTPTR!\n");
		return RBool::False();
	}
	
	//	if rule index is not an integer, return false or handle it somehow
	//
	if( RNumeric::isInt( ruleidx_sexp ) == false )
	{
		Rprintf("Expecting an integer as parameter 2 <ruleidx> in VCF_setRuleRefValues!\n");
		return RBool::False();
	}

	//	make sure the rule# is valid ( >= 0 and <= currently set number of rules)
	//
	int ruleidx = INTEGER( ruleidx_sexp )[0];
	if( f->num_rules_used < ruleidx || 0 > ruleidx )
	{
		Rprintf("(!!) Error : VCF_setRuleColumn : invalid rule-number to change\n");
		return RBool::False();
	}
	
	//	get pointer to wanted rule-definition
	//
	filter_entry_t * fe = &f->ruleset[ruleidx];
	
	//	convert and validate parameter to modify
	//
	if( RNumeric::isInt( ref1 ) && RNumeric::isInt( ref2 ) )
	{
		//make sure cmp op is int (or does not need ref values : DOES_EXIST)
		//
		if( ((fe->compare_how >= INT_CMP) && (fe->compare_how <= INT_CMP_CC)) || (fe->compare_how == DOES_EXIST) )
		{
			fe->ref.integers.i1 = INTEGER(ref1)[0];
			fe->ref.integers.i2 = INTEGER(ref2)[0];
		}
		else
		{
			Rprintf("(!!) Error : VCF_setRuleRefValues : integer comparisons used : specify reference values with as.integer()!\n");
			return RBool::False();
		}
	}
	//
	//
	else if( RNumeric::isFloat( ref1 ) && RNumeric::isFloat( ref2 ) )
	{
		//make sure cmp op is float (or does not need ref values : DOES_EXIST)
		//
		if( ((fe->compare_how >= FLT_CMP) && (fe->compare_how <= FLT_CMP_CC)) || (fe->compare_how == DOES_EXIST)  )
		{
			fe->ref.floats.f1 = REAL(ref1)[0];
			fe->ref.floats.f2 = REAL(ref2)[0];
		}
		else
		{
			Rprintf("(!!) Error : VCF_setRuleRefValues : floating point comparisons used : specify reference values with as.double()!\n");
			return RBool::False();
		}
	}
	else
	{
		Rprintf("(!!) Error : VCF_setRuleRefValues : reference values must both be either integer or double!\n");
		return RBool::False();
	}

	return RBool::True();
}





/*!	Set which subfield should be checked
**
**
*/
EXPORT SEXP VCF_setRuleField( SEXP vcfptr, SEXP ruleidx_sexp, SEXP field )
{

	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("(!!) Error : Parameter 1 not a VCFhandle EXTPTR!\n");
		return RBool::False();
	}
	
	//	if rule index is not an integer, return false or handle it somehow
	//
	if( RNumeric::isInt( ruleidx_sexp ) == false )
	{
		Rprintf("Expecting an integer as parameter 2 <ruleidx> in VCF_setRuleRefValues!\n");
		return RBool::False();
	}

	//	make sure the rule# is valid ( >= 0 and <= currently set number of rules)
	//
	int ruleidx = INTEGER( ruleidx_sexp )[0];
	if( f->num_rules_used < ruleidx || 0 > ruleidx )
	{
		Rprintf("(!!) Error : VCF_setRuleColumn : invalid rule-number to change\n");
		return RBool::False();
	}
	
	//	get pointer to wanted rule-definition
	//
	filter_entry_t * fe = &f->ruleset[ruleidx];
	
	//	convert and validate parameter to modify
	//
	const char * fieldname = RString::get( field , 0 );
	if( 0 == fieldname )
	{
		return RBool::False();
	}
	strncpy( & fe->filtered_subfield [0] , fieldname , sizeof( fe->filtered_subfield ) );

	return RBool::True();
}





/*!	Set a filtering rule's comparison operator
**
**
*/
EXPORT SEXP VCF_setRuleCmpOp( SEXP vcfptr, SEXP ruleidx_sexp, SEXP cmpop )
{

	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("(!!) Error : Parameter 1 not a VCFhandle EXTPTR!\n");
		return RBool::False();
	}
	
	//	if rule index is not an integer, return false or handle it somehow
	//
	if( RNumeric::isInt( ruleidx_sexp ) == false )
	{
		Rprintf("(!!) Error : Expecting an integer as parameter 2 <ruleidx> in VCF_setRuleCmpOp!\n");
		return RBool::False();
	}
	
	//	if rule index is not an integer, return false or handle it somehow
	//
	if( RNumeric::isInt( cmpop ) == false )
	{
		Rprintf("(!!) Error : Expecting an integer as parameter 3 <comparison-op> in VCF_setRuleCmpOp!\n");
		return RBool::False();
	}

	//	make sure the rule# is valid ( >= 0 and <= currently set number of rules)
	//
	int ruleidx = INTEGER( ruleidx_sexp )[0];
	if( f->num_rules_used < ruleidx || 0 > ruleidx )
	{
		Rprintf("(!!) Error : VCF_setRuleCmpOp : invalid rule-number to change\n");
		return RBool::False();
	}

	//	make sure the rule# is valid ( >= 0 and <= currently set number of rules)
	//
	int cmpopint = INTEGER( cmpop )[0];
	if( cmpopint > NUM_CMP_OPS || cmpopint < 0 )
	{
		Rprintf("(!!) Error : VCF_setRuleCmpOp : invalid comparison operator ( %d not in [0,%d] )\n",cmpopint,NUM_CMP_OPS-1);
		return RBool::False();
	}
	
	//	get pointer to wanted rule-definition
	//
	filter_entry_t * fe = &f->ruleset[ruleidx];
	
	//	convert and validate parameter to modify
	//
	fe->compare_how = cmpopint;

	return RBool::True();
}


//------------------------------------^^^^^^^^^^	DONE
//------------------------------------vvvvvvvvvv	TO DO

























