##
# Function sv_make_numeric().
# soeren.vogel@posteo.ch
# 2012-10-12, 2015-07-09
# Make numeric data frame from non-numeric (R-style) data frame.
# License: GPL >= 2
#
sv_make_numeric <- function( x ) {
	y <- NULL
	var.class <- class( x )[1]
	if ( var.class == "numeric" ) {
		y <- x
	}
	else if ( var.class == "logical" ) {
		y <- as.numeric( x )
	}
	else if ( var.class == "ordered" ) {
		y <- as.numeric( x )
	}
	else if ( var.class == "factor" ) {
		y <- as.numeric( x )
	}
	else {
		y <- as.character( x )
		y <- gsub( "\n", '<br>',  y, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE )
		y <- gsub( "\t", '<tab>', y, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE )
	}
	return( y )
}
