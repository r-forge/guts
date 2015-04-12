##
# GUTS R Definitions.
# soeren.vogel@uzh.ch
# carlo.albert@eawag.ch
# 2015-04-01
#





##
# A Little Helper For Formatted Printing.
#
.g_print_help <- function( x, width, digits ) {
	y <- paste( round(x, digits=digits), sep="", collapse=", " )
	z <- strwrap( y, width=width, indent=4, exdent=4 )
	return( z )
}





##
# Basic Hidden Print Functionality.
#
.guts.print <- function( object, width=getOption('width'), digits=getOption('digits') ) {

	# Header
	cat(
		"\n",
		"GUTS object with the following attributes:\n",
		"==========================================\n\n",
		sep=""
	)

	# Concentrations
	cat( "Concentrations (C, ", length(object$C), " elements)", sep="" )
	if ( length(object$C) > 0 ) {
		cat( ":", .g_print_help(object$C, width, digits), sep="\n" )
	} else {
		cat( "\n", sep="" )
	}

	# Concentration time points
	cat( "Concentration time points (Ct, ", length(object$Ct), " elements)", sep="" )
	if ( length(object$Ct) > 0 ) {
		cat( ":", .g_print_help(object$Ct, width, digits), sep="\n" )
	} else {
		cat( "\n", sep="" )
	}

	# Survivors
	cat( "Survivors (y, ", length(object$y), " elements)", sep="" )
	if ( length(object$y) > 0 ) {
		cat( ":", .g_print_help(object$y, width, digits), sep="\n" )
	} else {
		cat( "\n", sep="" )
	}

	# Survivor time points
	cat( "Survivor time points (yt, ", length(object$yt), " elements)", sep="" )
	if ( length(object$yt) > 0 ) {
		cat( ":", .g_print_help(object$yt, width, digits), sep="\n" )
	} else {
		cat( "\n", sep="" )
	}

	# Distribution, Model
	cat( "Distribution (dist)  : ", object$dist, "\n", sep="" )
	cat( "Model (model)        : ", object$model, "\n", sep="" )

	# Parameters
	cat( "Parameters (par, ", length(object$par), " elements)", sep="" )
	if ( length(object$par) > 0 ) {
		cat( ":", .g_print_help(object$par, width, digits), sep="\n" )
	} else {
		cat( "\n", sep="" )
	}

	# Time grid points
	cat( "Time grid points (M) : ", object$M, "\n", sep="" )
	#  cat( "Distribution (dist)  : ", dQuote(object$dist), "\n", sep="" )
	cat( "Sample length (N)    : ", object$N, "\n", sep="" )

	# Survival probabilities
	cat( "Survival probabilities (S, ", length(object$S), " elements)", sep="" )
	if ( length(object$S) > 0 ) {
		cat( ":", .g_print_help(object$S, width, digits), sep="\n" )
	} else {
		cat( "\n", sep="" )
	}

	# Loglikelihood
	cat( "Loglikelihood (LL): ", object$LL, "\n", sep="" )

	# Messages/warnings
	for ( msg in object$Messages ) {
		warning( msg )
	}
	cat( "\n", sep="" )

	# Footer
	cat( "\n", sep="" )
}

#show <- function( object ) {
#  .guts.print( object )
#}

# S3:
print.Rcpp_GUTS <- function(x, ...) {
  out <- .guts.print(x)
  return(invisible(out))
}

#print.GUTS <- function( x, ... ) {
#  out <- .guts.print( x )
#  return( invisible( out ) )
#}

# S4:
setMethod("show", "Rcpp_GUTS", function(object) .guts.print(object))



##
# Generic Function ‘logLik’.
#
logLik.Rcpp_GUTS <- function(object, ...) {
	eval.parent( substitute( object$calcLoglikelihood() ) )
	len <- length(object$y)
	if ( len < 1 | is.na(object$LL) | is.null(object$LL) ) {
		ret <- NA
		df <- NA
	} else {
		ret <- object$LL
		df <- length(object$y)
		df <- len
	}
	attr(ret, "df") <- df
	class(ret) <- "logLik"
	ret
}
