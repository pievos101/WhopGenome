/*
**
**		WhopGen
**
**		WHOle genome population genetics with popGEN
**
**
**		Tabix wrapper class implementation
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




/*!
**
*/
whop_tabix::whop_tabix(){

	//
	tabix = 0;
	index = 0;
	iter = 0;

	//
	num_seqnames = 0;
	sequence_names = 0;

	//
	last_num_fields = 0;
	field_offsets_size = 0;
	field_offsets = 0;
	current_line = 0;

	//
	//currentTid = ;
	currentBegin = 0;
	currentEnd = 0;

	//
}



/*!
**
*/
whop_tabix::whop_tabix(const char * filename ){

	//
	tabix = 0;
	index = 0;
	iter = 0;

	//
	num_seqnames = 0;
	sequence_names = 0;

	//
	last_num_fields = 0;
	field_offsets = 0;
	field_offsets_size = 0;
	current_line = 0;

	//
	//currentTid = 0;
	currentBegin = 0;
	currentEnd = 0;

	//
	this->open( filename );

	//
}



/*!
**
*/
whop_tabix::~whop_tabix()
{
	ONDBG Rprintf("whop_tabix DTOR: %lx\n",(uint64_t)tabix);

	//
	if( sequence_names != 0 )
	{
		free( sequence_names );	//were allocated with calloc() in tabix
		sequence_names=0;
		num_seqnames = 0;
	}

	//
	if( tabix != 0 )
	{
		ti_close( tabix );
		tabix = 0;
	}

	//
	if( iter != 0 )
	{
		ti_iter_destroy(iter);
		iter = 0;
	}

	//
	if( field_offsets != 0 )
	{
		last_num_fields = 0;
		field_offsets_size = 0;
		free( field_offsets );
		field_offsets = 0;
	}

	//
}



/*!
**
*/
bool			whop_tabix::open( const char * filename )
{
	//
	ti_iter_t hdriter=0;
	
	ONDBG Rprintf("whoptabix_open\n");

	//
	//
	try
	{
		//	load tabix-index of file
		//
		ONDBG Rprintf("ti_open(%s,0)...\n",filename,0);
		tabix = ti_open( filename , 0 );
		ONDBG Rprintf("ti_open=%d,OK\n",tabix);
		if( tabix == 0 )
			throw "whop_tabix::open : Failed to open tabix index file";

		//		load index
		//
		ONDBG Rprintf("ti_index_load(%s)",filename);
		index = ti_index_load( filename );
		ONDBG Rprintf("OK\n");
		if( index == 0 )
			throw "whop_tabix::open : index load failed!";

		//
		//
		ONDBG Rprintf("ti_query(%x,0,1,99999999)",tabix);
		iter = ti_query( tabix, 0, 1, 99999999 );//FIXME best way to have a catch-all iterator?
		ONDBG Rprintf("OK\n");
		if( 0 == iter )
			throw "whop_tabix::open : Failed to create whole-file iterator!";
		bEOR = (iter != 0 );

		//	iterator that includes header lines
		//
		ONDBG Rprintf("ti_query(%p,0,0,0)",tabix);
		hdriter = ti_query(tabix,0,0,0);
		ONDBG Rprintf("OK\n");
		if( 0 == hdriter )
			throw "whop_tabix::open : Failed to create header-iterator!";

		//	read header lines until the ##CHROM
		//
		int len=0;
		const char	*s;
		ONDBG Rprintf("ti_read(%p,%p,%p)*x",tabix,hdriter,&len);
		while(  (s=ti_read( tabix, hdriter, &len ))!=0 && (s[0]=='#') )
		{
			//Rprintf("[%s]\n",s);
			header_lines.push_back( s );
			
			if( s[1] == 'C' && s[2] == 'H' )	//stop parsing lines once the '#CHROM	POS	ID	REF	ALT	QUAL...' line is found
				break;							//		..helps not skipping the first data-line if no region was set
		}
		ONDBG Rprintf("OK\n");

		//	get sequence names (first column, usually chromosome name)
		//
		num_seqnames=0;
		sequence_names = ti_seqname( index, (signed*)&num_seqnames);

		ONDBG Rprintf("return TRUE from whoptabix_open\n");
#if 0
		//
		//
		for( unsigned int i = 0; i < num_seqnames; i++ )
		{
			printf("Range of nucleotides for '%s':\n",sequence_names[i] );
			
			//
			int min = 0;
			int max = 100*1000*1000, okmax=0, failmax=cmax;
			const char * tid = &sequence_names[i];
			ti_iter_t newiter = 0;
			int testiters=20;
			
			//	find maximum
			//
			while( testiters-- > 0 )
			{
				newiter = ti_query( tabix , tid, max, max+1 );
				
				// memorize which value works as maximum
				//
				if( newiter == 0 )	//half of max
				{
					failmax = max;
				}
				else
					okmax = max;
				
				// if there is too much difference between known-to-fail and known-to-work maximum
				//
				if( (failmax - okmax) > 1 )
				{
					if( newiter == 0 )
						
				}
				
				
			}
		}//...for( all contig-ids )
#endif

		return true;

		//
	}
	catch( const char * errmsg )
	{
		Rprintf("Caught exception inside whop_tabix::open('%s'):\n\t'%s'\n",filename,errmsg);
	}
	catch( ... )
	{
		Rprintf("Caught generic exception inside whop_tabix::open('%s')\n",filename);
	}

	//
	//
	if( sequence_names != 0 )
	{
		free( sequence_names );	//were allocated with calloc() in tabix
		sequence_names=0;
		num_seqnames = 0;
	}
	if( tabix != 0 )
		ti_close( tabix );
		tabix=0;
	if( hdriter )
		ti_iter_destroy( hdriter );
		hdriter=0;
	if( iter )
		ti_iter_destroy( iter );
		iter=0;
	//FIXME : free index ?
	
	Rprintf("return FALSE from whoptabix_open\n");

	return false;
}



		//-
		//-
		//-
		//-



/*!	Create a Tabix index for an already compressed VCF file 
**
**
*/
EXPORT  SEXP VCF_buildindex( SEXP filename_sexp )
{
	//
	if(! CHKSTR(filename_sexp,1) )
	{
		df0("VCF_buildindex : filename is not a single string!");
		return RBool::False();
	}
	const char * filename = CHAR(STRING_ELT(filename_sexp, 0));

	//
	if( 0 > ti_index_build( filename, &ti_conf_vcf) )
		return RBool::False();
	
	return RBool::True();
}


/*!	Open a Tabix-indexed VCF file and return a VCFhandle, a whopgen-managed EXTPTR SEXP
**
**
*/
EXPORT  SEXP tabix_build( SEXP filename_sexp, SEXP sc_sexp, SEXP bc_sexp, SEXP ec_sexp, SEXP meta_sexp, SEXP lineskip_sexp )
{
	//
	if(! CHKSTR(filename_sexp,1) ){	df0("tabix_build : filename is not a single string!");		return RBool::False();	}
	if(! CHKTYPE(sc_sexp,Integer,1) ){	df0("tabix_build : start column is not a single int!");		return RBool::False();	}
	if(! CHKTYPE(bc_sexp,Integer,1) ){	df0("tabix_build : begin column is not a single int!");		return RBool::False();	}
	if(! CHKTYPE(ec_sexp,Integer,1) ){	df0("tabix_build : end column is not a single int!");		return RBool::False();	}
	if(! CHKSTR(meta_sexp,1) ){	df0("tabix_build : meta is not a single string!");		return RBool::False();	}
	if(! CHKTYPE(lineskip_sexp,Integer,1) ){	df0("tabix_build : lineskip is not a single int!");		return RBool::False();	}
	
	//
	//	
	const char * filename = CHAR(STRING_ELT(filename_sexp, 0));
	int sc = INTEGER(sc_sexp)[0];
	int bc = INTEGER(bc_sexp)[0];
	int ec = INTEGER(ec_sexp)[0];
	char meta = CHAR(STRING_ELT(meta_sexp, 0))[0];
	int lineskip = INTEGER(lineskip_sexp)[0];
	
//	Rprintf("tabix_build %d %d %d %c %d\n",sc,bc,ec,meta,lineskip);
	
	ti_conf_t	conf;
	conf.preset = 0;
	conf.sc = sc;
	conf.bc = bc;
	conf.ec = ec;
	conf.meta_char = meta;
	conf.line_skip = lineskip;
	
	int res = ti_index_build( filename, &conf);
	if( res < 0 )
		return RBool::False();

	return RBool::True();
}



/*!
**
**
*/
EXPORT  SEXP BGZF_compress( SEXP in_filename_sexp, SEXP out_filename_sexp )
{
	//
	if(! CHKSTR(in_filename_sexp,1) ){	df0("VCF_open : in-filename is not a single string!");		return RBool::False();	}
	if(! CHKSTR(out_filename_sexp,1) ){	df0("VCF_open : out-filename is not a single string!");		return RBool::False();	}

	//
	const char * in_filename = CHAR(STRING_ELT(in_filename_sexp, 0));
	const char * out_filename = CHAR(STRING_ELT(out_filename_sexp, 0));
	
	//
	df2("Compressing file '%s' into file '%s'\n",in_filename,out_filename);
	
	//
	FILE * infh = fopen(in_filename,"rb");
	if( 0 == infh )
	{
		df0("(!!) bgzf_compress : Could not open for reading: input file '%s'\n",in_filename);
		return RBool::False();
	}
	
	//	open output file
	//
	BGZF * outfh = bgzf_open( out_filename, "w8" );
	if( 0 == outfh )
	{
		df0("(!!) bgzf_compress : Could not open for compressed writing: output file '%s'\n",out_filename);
		fclose(infh);
		return RBool::False();
	}
	
	//	determine filesize
	//
	fseek( infh, 0, SEEK_END );
	off_t fsize = ftello(infh);
	df2("(II) bgzf_compress : Input file is %ld bytes long\n",fsize);
	fseek( infh, 0, SEEK_SET );
	
	//	read + bgzf-write loop
	//	
	//
	bool		failed=false;
	char		buffer[ 4096 ];
	off_t		curpos=0;
	uint64_t	times_nothing=0;		//count how many repeated calls to bgzf_write did not write any data
				
	//	how many calls are allowed to produce no compression? - 20*4096 = 80KB
#define		MAX_NONCOMPRESSING_BGZF_WRITE_CALLS		20

	//
	while( curpos < fsize )
	{
		
		//	read from input file
		//
		uint64_t rdb = fread( buffer, 1, 4096, infh );
		if( 0 >= rdb )
		{
			//df2("READ %d bytes, curpos now %d -- breaking loop \n",rdb,curpos);
			break;
		}
		curpos += rdb;
		
		//	write and compress into output file
		//
		int r = bgzf_write( outfh, buffer, rdb );
		if( 0 >= r )
		{
			times_nothing ++;	// detect if bgzf_write does not write any compressed data anymore
		}
		else
			times_nothing = 0;
			
		//	detect read and compression errors, and handle them
		//
		if( (-1 == r) || (times_nothing > MAX_NONCOMPRESSING_BGZF_WRITE_CALLS) )
		{
			Rprintf("(!!) Aborting compression : compress result %d after %d read/compress/write repeats\n",r,times_nothing);
			failed = true;
			break;
		}
		//df2("READ %d bytes, curpos now %d, wrote with result %d\n",rdb,curpos,r);
		

	}//...while( data left to compress )
	
	//
	fclose(infh);		infh=0;
	bgzf_close(outfh);	outfh=0;
	
	//
	return failed ? RBool::False() : RBool::True();
}


