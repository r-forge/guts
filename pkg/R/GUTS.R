# GUTS R definitions.
# carlo.albert@eawag.ch, soeren.vogel@uzh.ch
# Version 0.2.2
# GPL-2

#
# Method show, print
#
.guts.print <- function( x ) {
  cat(
    "\n",
    "GUTS object with the following attributes:\n\n",

    "Vector of concentrations (C, ", length(x$C), " elements):\n",
    "\t", paste(x$C, sep="", collapse=", "), "\n",
    
    "Vector of concentration time points (Ct, ", length(x$Ct), " elements):\n",
    "\t", paste(x$Ct, sep="", collapse=", "), "\n",

    "Vector of survivors (y, ", length(x$y), " elements):\n",
    "\t", paste(x$y, sep="", collapse=", "), "\n",

    "Vector of survivor time points (yt, ", length(x$yt), " elements):\n",
    "\t", paste(x$yt, sep="", collapse=", "), "\n",

    "Parameters (par, ", length(x$par), " elements):\n",
    "\t", paste(x$par, sep="", collapse=", "), "\n",

    "Time grid points (M) : ", x$M, "\n",
    "Distribution (dist)  : ", dQuote(x$dist), "\n",
    "Sample length (N)    : ", x$N, "\n",
    
    "Sample vector (z, ", length(x$z), " elements):\n",
    "\t", "Min=", min(x$z), ", Max=", max(x$z), ", Mean=", mean(x$z), ", Sigma=", sd(x$z), "\n",
    
    "\n\n",
    sep=""
  )
}
# S3:
print.GUTS <- function( x, ... ) {
  out <- .guts.print( x )
  return( invisible( out ) )
}
# S4:
#setMethod("show", "Rcpp_GUTS", function( object ) .guts.print( object ) )
setMethod("show", "missing", function( object ) .guts.print( object ) )

#
# Generic function logLik
#
logLik.GUTS <- function (object, ...) {
  len <- length(object$y)
  if ( len < 1 | is.na(object$LL) | is.null(object$LL) ) {
    ret <- NA
    df <- NA
  }
  else {
    ret <- object$LL
    df <- len
  }
  attr(ret, "df") <- df
  class(ret) <- "logLik"
  ret
}
