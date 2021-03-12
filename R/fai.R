
#
#
#
fai_open <- function( filename ) .Call("FAI_open", filename, PACKAGE="WhopGenome" )
fai_close <- function( faifh ) .Call("FAI_close", faifh, PACKAGE="WhopGenome")
fai_reopen <- function( faifh ) .Call("FAI_reopen", faifh, PACKAGE="WhopGenome")
fai_build <- function( filename ) .Call("FAI_build", filename, PACKAGE="WhopGenome")
fai_query2 <- function( faifh, regionstring ) .Call("FAI_query2", faifh, regionstring, PACKAGE="WhopGenome")
fai_query4 <- function( faifh, sequencename, beginpos, endpos ) .Call("FAI_query4", faifh, sequencename, beginpos, endpos, PACKAGE="WhopGenome")

fai_query3 <- function( faifh, regionstring, resultstring ){
	warning("DEPRECATED : fai_query3 : Please use fai_query2 which returns the resultstring directly or FALSE")
	.Call("FAI_query3", faifh, regionstring, resultstring, PACKAGE="WhopGenome")
}
fai_query5 <- function( faifh, sequencename, beginpos, endpos, resultstring ){
	warning("DEPRECATED : fai_query5 : Please use fai_query4 which returns the resultstring directly or FALSE")
	.Call("FAI_query5", faifh, sequencename, beginpos, endpos, resultstring, PACKAGE="WhopGenome")
}
