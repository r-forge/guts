##
# Function sv_get_levels().
# soeren.vogel@posteo.ch
# 2014-03-23
# Get levels for a variable name from codebook.
# License: GPL >= 2
#
sv_get_levels <- function( varname, codebook ) {
	levs <- unlist( strsplit( codebook[varname, "Levels"], "; " ) )
	for ( i in 1:length(levs) ) {
		if ( levs[i] == "" ) {
			levs[i] <- i
		}
	}
	return( levs )
}
