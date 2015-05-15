##
# GUTS R Definitions.
# soeren.vogel@uzh.ch, carlo.albert@eawag.ch
# License GPL-2
# 2015-05-15
#





##
# Function guts_setup(...).
#
guts_setup <- function(C, Ct, y, yt, par = 0, dist = 'lognormal', model = 'Proper', N = 1000, M = 10000) {

	#
	# Defaults for S and LL.
	#
	S   <-  NA
	LL  <-  NA

	#
	# Check missing arguments and arguments types (numeric, character).
	#
	args_num_names  <- c('C', 'Ct', 'y', 'yt', 'par', 'N', 'M')
	args_char_names <- c('dist', 'model')
	args_num_type   <- c(is.numeric(C), is.numeric(Ct), is.numeric(y), is.numeric(yt), is.numeric(par), is.numeric(N), is.numeric(M))
	args_char_type  <- c(is.character(dist), is.character(model))
	if ( any( !args_num_type ) ) {
		i <- which(!args_num_type)[1]
		stop( paste( "Argument ", args_num_names[i], " not numeric.", sep='' ) )
	} else if ( any( !args_char_type ) ) {
		i <- which(!args_num_type)[1]
		stop( paste( "Argument ", args_char_names[i], " not character.", sep='' ) )
	}

	#
	# Check length of single value arguments.
	#
	args_sin_names  <- c('dist', 'model', 'N', 'M')
	args_sin_len    <- c(length(dist), length(model), length(N), length(M))
	for ( i in seq_along(args_sin_len) ) {
		if ( args_sin_len[i] > 1 ) {
			warning( paste( "Argument ", args_sin_names[i], " must be of length 1, only first element used.", sep='' ) )
			assign( args_sin_names[i], get(args_sin_names[i])[1] )
		}
	}

	#
	# Check concentrations and survivors.
	#
	if ( length(C) < 2 ) {
		stop( 'Vector C must be longer than 1.' )
	} else if ( length(Ct) < 2 ) {
		stop( 'Vector Ct must be longer than 1.' )
	} else if ( length(C) != length(Ct) ) {
		stop( 'Vectors C and Ct must have the same length.' )
	} else if ( Ct[1] != 0 ) {
		stop( 'Vector Ct must start at 0.0.' )
	} else if ( any(diff(Ct) <= 0) ) {
		stop( 'Vector Ct must contain unique values in ascending order.' )
	} else if ( length(y) < 2 ) {
		stop( 'Vector y must be longer than 1.' )
	} else if ( length(yt) < 2 ) {
		stop( 'Vector yt must be longer than 1.' )
	} else if ( length(y) != length(yt) ) {
		stop( 'Vectors y and yt must have the same length.' )
	} else if ( yt[1] != 0 ) {
		stop( 'Vector yt must start at 0.0.' )
	} else if ( any(diff(yt) <= 0) ) {
		stop( 'Vector yt must contain unique values in ascending order.' )
	} else if ( any(diff(y) > 0) ) {
		stop( 'Values in vector y must not ascend.' )
	} else if ( Ct[length(Ct)] < yt[length(yt)] ) {
		stop( 'Last value of vector Ct must not be smaller than last value of vector yt.' )
	} else if ( min(c(C, Ct, y, yt)) < 0 ) {
		stop( 'Vectors C, Ct, y, yt must contain non-negative values.' )
	}

	#
	# Check dist and model, set experiment code, set par_pos.
	#
	experiment <- 11 # Default
	par_pos <- 1:5
	mdist  <- tolower(dist)
	mmodel <- tolower(model)
	if ( mdist == 'lognormal' ) {
		if ( mmodel == 'proper' ) {
			experiment <- 11
			par_pos <- 1:5
		} else if ( mmodel == 'it' ) {
			experiment <- 12
			par_pos <- c(1:2, 4:5)
		} else {
			stop( 'Model can be either "Proper" or "IT".' )
		}
	} else if ( mdist == "delta" ) {
		if ( mmodel == 'proper' ) {
			experiment <- 21
			par_pos <- 1:4
		} else if ( mmodel == 'it' ) {
			experiment <- 22
			par_pos <- c(1:2, 4)
		} else {
			stop( 'Model can be either "Proper" or "IT".' )
		}
	} else {
		stop( 'Distribution can be either "lognormal" or "delta".' )
	}

	#
	# Check parameters.
	# This is the general setup check. Function guts_calc_loglikelihood() will only check values.
	# Set wpar in attributes.
	#
	wpar <- c(0, 0, .Machine$double.xmax, 0, 0) # Default
	if ( length(par) == 1 && par == 0 ) {
		par <- rep( 0.0, length(par_pos) )
	} else if ( length(par) != length(par_pos) ) {
		stop( 'Wrong parameter vector length.' )
	} else if ( min(par) < 0 ) {
		# Here it is an error, in update it is a warning!
		stop( 'Parameters must be non-negative values.' )
	} else {
		# Assign to wpar at positions in par_pos
		wpar[par_pos] <- par
		if ( wpar[4] == 0 && wpar[5] != 0 ) {
			# Here it is an error, in update it is a warning!
			stop( 'Sample sd = 0 only if m = 0.' )
		} else if ( is.infinite(wpar[3]) ) {
			wpar[3] = .Machine$double.xmax;
		}
	}

	#
	# Check sample length and time grid points.
	#
	if ( N < 3 ) {
		stop( 'N must be greater than 2.' )
	} else if ( M < 2 ) {
		stop( 'M must be greater than 1.' )
	}

	#
	# Build GUTS object for return.
	#
	ret <- structure(
		list(
			'C'     = C,
			'Ct'    = Ct,
			'y'     = y,
			'yt'    = yt,
			'par'   = par,
			'dist'  = dist,
			'model' = model,
			'N'     = N,
			'M'     = M,
			'results' = list(
				'S'  =  S,
				'LL' =  LL
			)
		),
		class      = "GUTS",
		experiment = experiment,
		wpar       = wpar,
		par_pos    = par_pos,
		checked    = TRUE # Refers to whether the object was modified only with setup function.
	)
	invisible( return( ret ) )
}





##
# Calc Loglik
#
#guts_calc_loglikelihood <- function (gobj, pars) invisible(.Primitive(".Call")(<pointer: 0x110bc95e0>, gobj, pars))





##
# Function guts_calc_survivalprobs(...).
#
guts_calc_survivalprobs <- function( gobj, yts = 0 ) {

	#
	# Checks.
	#
	if ( !inherits( gobj, "GUTS", which=FALSE ) ) {
		stop( "No GUTS object. Use `guts_setup()` to create or modify objects." )
	} else if ( !is.numeric(yts) || length(yts) < 1 ) {
		stop( "Wrong value for argument yts." )
	}

	#
	# Create or use y and yt.
	#
	if ( length(yts) == 1 ) {
		# existing (0) or length.
		if ( yts < 0 ) {
			stop( "Wrong value for argument yts." )
		} else if ( yts == 0 ) {
			# existing y and yt
			y  <- gobj$y
			yt <- gobj$yt
		} else {
			# number of yts.
			y  <- rep(0.0, yts)
			yt <- 0:(yts-1)
		}
	} else {
		# time points.
		y  <- rep(0.0, length(yts))
		yt <- yts
	}

	#
	# Use guts_setup to check for correct vector.
	#
	# tryCatch(print(a), error = function(e) FALSE)
	tmp <- guts_setup(gobj$C, gobj$Ct, y, yt, par = gobj$par, dist = gobj$dist, model = gobj$model, N = gobj$N, M = gobj$M)
	if ( is.list(tmp) && class(tmp)=="GUTS" ) {
		# No error.
		# a <- 1:3; a
		# f <- function(obj) {
		# 	objname <- deparse(substitute(obj)); print(objname)
		# 	assign(objname, 7:9, parent.frame()); print(obj)
		# }
		# f(a); a

		# Calculate everything on tmp.
		guts_calc_loglikelihood( tmp, gobj$par )

		# Create a random hidden temp file in parent.frame() and assign tmp to it.
		# Because assign does not work with lists.
		fn <- paste( '.', paste(sample(letters, 8), sep='', collapse=''), as.numeric(Sys.time()), sep='' )
		assign(fn, tmp, parent.frame())

		# Now do the actual assignment of fn to the original object.
		gobjname <- deparse(substitute(gobj))
		txt <- paste( gobjname, '<-', fn, sep='' )
		eval(parse(text=txt), parent.frame())

		# Remove temporary file from parent.frame().
		rm( list=c(fn), envir=parent.frame() )
	} else {
		stop( "" )
	}
}





##
# Printing.
#

#
# A small helper for printing and wrapping.
#
.g_print_help <- function( x, width, digits, prefix=NULL ) {
	y <- paste( round(x, digits=digits), sep="", collapse=", " )
	z <- strwrap( y, width=width, indent=0, exdent=6, initial=prefix )
	return( z )
}

#
# The actual print function.
#
.g_print <- function( object, width=getOption('width'), digits=getOption('digits') ) {

	# Header
	cat(
		"\n",
		"GUTS object:\n",
		"============\n",
		sep=""
	)

	# Distribution, Model
	cat( "Distribution: ", object$dist, ", model: ", object$model, ".\n", sep="" )

	# Concentrations, Survivors
	cat( "Concentrations (n=", length(object$C), "), survivors (n=", length(object$y), ")", sep="" )
	if ( length(object$C) > 0 ) {
		cat( ":", sep="\n" )
		cat( .g_print_help(object$Ct, width, digits, prefix="  Ct: "), sep="\n" )
		cat( .g_print_help(object$C,  width, digits, prefix="   C: "), sep="\n" )
	}
	if ( length(object$y) > 0 ) {
		cat( .g_print_help(object$yt, width, digits, prefix="  yt: "), sep="\n" )
		cat( .g_print_help(object$y,  width, digits, prefix="   y: "), sep="\n" )
	} else {
		cat( "\n", sep="" )
	}

	# Sample length, Time grid points
	cat( "Sample length: ", object$N, ", Time grid points: ", object$M, ".\n", sep="" )

	# Parameters
	cat( "Parameters (n=", length(object$par), ")", sep="" )
	if ( length(object$par) > 0 ) {
		cat( .g_print_help(object$par, width, digits, prefix=": "), sep="\n" )
	} else {
		cat( "\n", sep="" )
	}

	# Survival probabilities
	cat( "Survival probabilities (n=", length(object$results$S), ")", sep="" )
	if ( length(object$results$S) > 0 ) {
		cat( .g_print_help(object$results$S, width, digits, prefix=": "), sep="\n" )
	} else {
		cat( "\n", sep="" )
	}

	# Loglikelihood
	cat( "Loglikelihood: ", object$results$LL, "\n", sep="" )

	# Footer
	cat( "\n", sep="" )
}

#
# Do print.
#
print.GUTS <- function(x, ...) {
	out <- .g_print(x)
	return(invisible(out))
}





##
# Setters of Fields.
#
#"[[<-" <- function(x, field, value) {
#	UseMethod("[[<-", x)
#}
"[[<-.GUTS" <- function(x, field, value) {
	stop( "Use function `guts_setup()` for changing fields." )
}
#"$<-" <- function(x, field, value) {
#	UseMethod("$<-", x)
#}
"$<-.GUTS" <- function(x, field, value) {
	stop( "Use function `guts_setup()` for changing fields." )
}





##
# Attributes.
#
#'attr<-' <- function(x, which, value) {
#	UseMethod('attr<-', x)
#}
"attr<-.GUTS" <- function(x, which, value) {
	stop( "Use function `guts_setup()` for changing GUTS object fields." )
}
#'attributes<-' <- function(x, which, value) {
#	UseMethod('attributes<-', x)
#}
"attributes<-.GUTS" <- function(x, which, value) {
	stop( "Use function `guts_setup()` for changing attributes of a GUTS object." )
}
#'mostattributes<-' <- function(x, which, value) {
#	UseMethod('mostattributes<-', x)
#}
"mostattributes<-.GUTS" <- function(x, which, value) {
	stop( "Use function `guts_setup()` for changing attributes of a GUTS object." )
}




































