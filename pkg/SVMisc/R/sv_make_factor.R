##
# Function sv_make_factor().
# soeren.vogel@posteo.ch
# 2014-03-23
# Make a factor.
# License: GPL >= 2
#
sv_make_factor <- function( varname, numerics, codebook=NULL, add=0, ordered=FALSE, levels=NULL ) {
	if ( !is.null(codebook) ) {
		levs <- sv_get_levels( varname, codebook )
	}
	else if ( !is.null(levels) ) {
		levs <- levels
	}
	else {
		levs <- unique(as.character(sort(numerics)))
	}
	nlevs <- length(levs)
	vals <- c( 1:length(levs), numerics + add )
	fact <- factor( vals, ordered=ordered, labels=levs )
	return( fact[(nlevs+1):(length(fact))] )
}
