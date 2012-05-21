# GUTS R definitions.
# carlo.albert@eawag.ch, soeren.vogel@uzh.ch
# 2012-05-21
# GPL-2

#
# basic hidden print functionality
#
.guts.print <- function(object, digits=6) {
  cat(
    "\n",
    "GUTS object with the following attributes:\n",
    "==========================================\n\n",
    sep=""
  )

  cat( "Vector of concentrations (C, ", length(object$C), " elements)", sep="" )
  if ( length(object$C) > 0 ) {
    cat( ":\n\t", paste(object$C, sep="", collapse=", "), sep="" )
  }
  cat( "\n", sep="" )

  cat( "Vector of concentration time points (Ct, ", length(object$Ct), " elements)", sep="" )
  if ( length(object$Ct) > 0 ) {
    cat( ":\n\t", paste(object$Ct, sep="", collapse=", "), sep="" )
  }
  cat( "\n", sep="" )

  cat( "Vector of survivors (y, ", length(object$y), " elements)", sep="" )
  if ( length(object$y) > 0 ) {
    cat( ":\n\t", paste(object$y, sep="", collapse=", "), sep="" )
  }
  cat( "\n", sep="" )

  cat( "Vector of survivor time points (yt, ", length(object$yt), " elements)", sep="" )
  if ( length(object$yt) > 0 ) {
    cat( ":\n\t", paste(object$yt, sep="", collapse=", "), sep="" )
  }
  cat( "\n", sep="" )

  cat( "Parameters (par, ", length(object$par), " elements)", sep="" )
  if ( length(object$par) > 0 ) {
    cat( ":\n\t", paste(object$par, sep="", collapse=", "), sep="" )
  }
  cat( "\n", sep="" )

  cat( "Time grid points (M) : ", object$M, "\n", sep="" )
#  cat( "Distribution (dist)  : ", dQuote(object$dist), "\n", sep="" )
  cat( "Distribution (dist)  : ", object$dist, "\n", sep="" )
  cat( "Sample length (N)    : ", object$N, "\n", sep="" )

  cat( "Sample vector (z, ", length(object$z), " elements)", sep="" )
  if ( length(object$z) > 0 ) {
    cat( ":\n\t", "Min=", min(object$z), ", Max=", max(object$z), ", Mean=", mean(object$z), ", Sigma=", sd(object$z), sep="" )
  }
  cat( "\n", sep="" )

  cat( "Vector of damages (D, ", length(object$D), " elements)", sep="" )
  if ( length(object$D) > 0 ) {
    cat( ":\n\t", "Min=", min(object$D), ", Max=", max(object$D), ", Mean=", mean(object$D), ", Sigma=", sd(object$D), sep="" )
  }
  cat( "\n", sep="" )

  cat( "Vector of survival probabilities (S, ", length(object$S), " elements)", sep="" )
  if ( length(object$S) > 0 ) {
    cat( ":\n\t", paste( round(object$S, digits=digits), sep=", ", collapse=", " ), sep="" )
  }
  cat( "\n", sep="" )

  cat( "Loglikelihood (LL) : ", object$LL, "\n", sep="" )

  cat( "Messages/warnings: ", sep="" )
  if ( sum(object$errors!="") == 0 ) {
    cat( "none", "\n", sep="" )
  }
  else {
    cat( "\n", sep="" )
    print( object$errors[object$errors!=""] )
    cat( "\n", sep="" )
  }

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
#
setMethod("show", "Rcpp_GUTS", function(object) .guts.print(object))

#
# Generic function logLik
#
logLik.Rcpp_GUTS <- function(object, ...) {
  eval.parent( substitute( object$calcLoglikelihood() ) )
  len <- length(object$y)
  if ( len < 1 | is.na(object$LL) | is.null(object$LL) ) {
    ret <- NA
    df <- NA
  }
  else {
    ret <- object$LL
    df <- length(object$y)
    df <- len
  }
  attr(ret, "df") <- df
  class(ret) <- "logLik"
  ret
}
