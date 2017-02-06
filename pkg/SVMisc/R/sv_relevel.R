##
# Function sv_relevel().
# soeren.vogel@posteo.ch
# 2014-03-23
# Relevel a factor.
# License: GPL >= 2
#
sv_relevel <- function( factor, levels ) {
	levs <- rev(levels)
	for ( i in levs ) {
		out <- relevel( factor, ref=i )
	}
	return( out )
}
