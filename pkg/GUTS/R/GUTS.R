##
# GUTS R Definitions.
# soeren.vogel@uzh.ch, carlo.albert@eawag.ch, oliver.jakoby@rifcon.de, alexander.singer@rifcon.de, dirk.nickisch@rifcon.de
# License GPL-2
# 2019-05-24
# updated: 2021-11-30


##
# Function guts_setup(...).
#
guts_setup <- function(
	C, Ct, y, yt,
	dist = 'lognormal', model = 'Proper',
	N = 1000L, 
	MF = 100L, 
	M = max(
		5000L, 
		as.integer(ceiling(MF * length(union(Ct, yt)))),
		as.integer(ceiling(MF * max(union(Ct, yt))))
		),
	SVR = 1L,
	study = "", Clevel = ""
) {
	
	#
	# Check missing arguments and arguments types (numeric, character).
	#
	args_num_names  <- c('C', 'Ct', 'y', 'yt', 'N', 'M', 'SVR')
	args_char_names <- c('dist', 'model', 'study', 'Clevel')
	if (length(y) == 1) {
		if (is.na(y) | is.null(y)) y <- numeric()
	}
	if (is.na(M) | is.null(M)) M <- as.numeric(NA)
	if (is.na(N) | is.null(N)) N <- as.numeric(NA)
	if (any(is.na(SVR), is.nan(SVR), is.null(SVR), is.infinite(SVR))) SVR <- 1L
	args_num_type   <- c(is.numeric(C), is.numeric(Ct), is.numeric(y), is.numeric(yt), is.numeric(N), is.numeric(M), is.numeric(SVR))
	args_char_type  <- c(is.character(dist), is.character(model), is.character(study), is.character(Clevel))
	if ( any( !args_num_type ) ) {
		i <- which(!args_num_type)
		stop( paste( "Argument ", paste0(args_num_names[i], collapse = ", "), " must be numeric.", sep='' ) )
	} else if ( any( !args_char_type ) ) {
		i <- which(!args_char_type)
		stop( paste( "Arguments ", paste0(args_char_names[i], collapse = ", "), " must be character.", sep='' ) )
	}
	
	check_finite(Ct)
	check_finite(C)
	check_finite(yt)
	
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
	
	
	# Check concentrations and survivors.
	
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
	} else if ( min(c(C, Ct, y, yt)) < 0 ) {
		stop( 'Vectors C, Ct, y, yt must contain non-negative values.' )
	}
	
	
	#
	# Check Ct and yt length and, if needed, truncate y and yt.
	#
	Ct.last <- Ct[length(Ct)]
	if ( Ct.last < yt[length(yt)] ) {
		i <- which(yt <= Ct.last)
		y <- y[i]
		yt <- yt[i]
		warning( 'Survivor information at time points later than the latest concentration time point are disregarded.' )
	}
	
	#list reflects enums in C++
	TD_types <- list(PROPER = 0L, IT = 1L, SD = 2L )
	dist_types <- list(LOGLOGISTIC = 0L, LOGNORMAL = 1L, DELTA = 2L, EXTERNAL = 3L)
	
	TD <- toupper(model)
	dist_type <- toupper(dist)
	
	#organise parameters
	
	par_len <- switch (TD,
										 PROPER =
										 	switch(dist_type,
										 				 LOGLOGISTIC =, LOGNORMAL = 5L,
										 				 DELTA = 4L,
										 				 EXTERNAL = 3L
										 	),
										 IT =
										 	switch(dist_type,
										 				 LOGLOGISTIC =, LOGNORMAL = 4L,
										 				 EXTERNAL = 2L
										 	),
										 SD = 4L
	)
	if (is.null(par_len)) {
		stop("Cannot construct GUTS from model = '", model, "' and dist = '", dist, "'")
	}
	
	if (any(is.na(M), is.nan(M), is.infinite(M), is.null(M))) {
		D <- NA
		Dt <- NA
	} else {
		D <- rep(NA, M)
		Dt <- rep(NA, M)
	}
	
	#
	# Check correctness of model specific parameters
	#
	if (TD %in% c("PROPER", "SD")) {
		if (any(is.na(M), is.nan(M), is.infinite(M), is.null(M), M<2)) {
			stop(
				paste0(
					"For models 'PROPER' and 'SD':\n",
					"  The number of discretization time steps M must be an integer >= 2.\n",
					"  Current value: ", M
				)
			)
		}
	}
	if (TD == "PROPER") {
		if (any(is.na(N), is.nan(N), is.infinite(N), is.null(N), N<3)) {
			stop(
				paste0(
					"For model 'PROPER':\n",
					"  The distribution sample size N must be an integer >= 3.",
					"  Current value: ", N
				)
			)
		}
	}
	if (SVR<=0) {
		stop(
			paste0(
				"The surface volume ratio SVR must be a positive value.",
				"  Current value: ", SVR
			)
		)
	}
	
	#
	# Build GUTS object for return.
	#
	
	ret <- structure(
		list(
			'study' = study,
			'Clevel'= Clevel,
			'C'     = C,
			'Ct'    = Ct,
			'y'     = y,
			'yt'    = yt,
			'dist'  = dist,
			'model' = model,
			'N'     = N,
			'M'     = M,
			'par'   = rep(NA, par_len),
			'external_dist' = NULL,
			'S'     = rep(NA, length(yt)),
			'D'     = D,
			'Dt'    = Dt,
			'LL'    = NA,
			'SPPE'  = NA,
			'squares' = NA,
			'SVR'   = SVR
		),
		class      = "GUTS",
		TD_type    = TD_types[[TD]],
		dist_type  = dist_types[[dist_type]],
		par_len    = par_len,
		update_ID  = c(S = 0, SPPE = -1, squares = -1)
	)
	invisible( return( ret ) )
} # End of guts_setup()

check_finite <- function(x) {
	if (any(is.na(x), is.nan(x), is.infinite(x), is.null(x))) {
		stop(
			paste0(
				deparse(substitute(x)),
				" needs to be finite."
			)
		)
	}
}

##
# Function guts_calc_loglikelihood(...).
guts_calc_loglikelihood <- function(gobj, par, external_dist = NULL, use_multinomial_coefficient = FALSE) {
	invisible(.Call('_GUTS_guts_engine', PACKAGE = 'GUTS', gobj, par, z_dist = external_dist))
	if (use_multinomial_coefficient) {
		return(gobj[['LL']] + log_multinomial_coefficient(gobj))
	} else {
		return(gobj[['LL']])
	}
}

##
# Function guts_calc_survivalprobs(...).
guts_calc_survivalprobs <- function(gobj, par, external_dist = NULL) {
	invisible(.Call('_GUTS_guts_engine', PACKAGE = 'GUTS', gobj, par, z_dist = external_dist))
	return(gobj[['S']])
}

##
# Function guts_report_damage(...).
guts_report_damage <- function(gobj) {
	if (gobj$model == "IT") {
		# To get a more complete damage profile, the fast projector internally calculates damage at some relevant time points.
		# These additional calculations are appended and might be duplicates of previous calculations
		# Therefore, duplicate time points are removed and the order is corrected.
		dam_time <- gobj[['Dt']]
		dam_time_sorted <- sort.int(unique(dam_time))
		return(
			data.frame(
				time = dam_time_sorted,
				damage = gobj[['D']][match(dam_time_sorted, dam_time)]
			)
		)
	} else {
		return(
			data.frame(
				time = gobj[['Dt']],
				damage = gobj[['D']]
			)
		)
	}
}

##
# Function guts_report_sppe(...).
guts_report_sppe <- function(gobj) {
	return(gobj[['SPPE']])
}

##
# Function guts_report_squares(...).
guts_report_squares <- function(gobj) {
	return(gobj[['squares']])
}

###
# multinomial coefficients
faculty <- function(x) sapply(x, function(y) prod(seq_len(y)))
multinomial_coefficient <- function(gts) faculty(sum(-diff(c(gts$y,0))))/prod(faculty(-diff(c(gts$y,0))))
log_multinomial_coefficient <- function(gts) log(faculty(sum(-diff(c(gts$y,0))))) - log(prod(faculty(-diff(c(gts$y,0)))))


##
# Printing.
#

# A small helper for printing and wrapping.
.g_print_help <- function( x, width, digits, exdent = 6, prefix = NULL ) {
	y <- paste( round(x, digits = digits), sep = "", collapse = ", " )
	z <- strwrap( y, width = (width-6), indent = 0, exdent = exdent, initial = prefix )
	return( z )
}

# The actual print function.
.g_print <- function( object, width = getOption('width'), digits = getOption('digits') ) {
	
	# Header
	cat(
		"\n",
		"GUTS object:\n",
		"============\n",
		sep=""
	)
	
	# Distribution, Model
	cat( "Name: ", object$study, ", CLevel: ", object$Clevel, ".\n", sep="" )
	
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
	prf <- paste("Parameters (n=", length(object$par), ")", sep="")
	if ( length(object$par) > 0 ) {
		prf <- paste(prf, ": ", sep="")
		cat( .g_print_help(object$par, width, digits, prefix=prf), sep="\n" )
	} else {
		cat( "\n", sep="" )
	}
	
	# External distribution
	if (is.null(object$external_dist)) {
		cat("External distribution: NULL\n", sep = "")
	} else {
		prf <- paste0("External distribution (n=", length(object$external_dist), "; subset: first 6 elements)")
		if ( length(object$external_dist) > 0 ) {
			prf <- paste(prf, ": ", sep="")
			cat( .g_print_help(head(object$external_dist), width, digits, prefix=prf), sep="\n" )
		} else {
			cat( "\n", sep="" )
		}
	}
	
	# Survival probabilities
	prf <- paste("Survival probabilities (n=", length(object$S), ")", sep="")
	if ( length(object$S) > 0 ) {
		prf <- paste(prf, ": ", sep="")
		cat( .g_print_help(object$S, width, digits, exdent = 2, prefix = prf), sep="\n" )
	} else {
		cat( "\n", sep="" )
	}
	
	# Damage
	if (any(is.na(object$D))) {
		damage <- object$D
	} else {
		damage <- guts_report_damage(object)$damage
	}
	if (length(damage) == 0) {
		cat("\n", sep="")
	} else {
		prf <- paste0("Damage (n=", length(damage))
		if (length(damage) > 30) {
			prf <- paste0(prf, "; subset: 30 regularly spaced values)")
			damage <- damage[seq(1, length(damage), length.out = 30)]
		} else {
			prf <- paste0(prf, ")")
		}
		prf <- paste(prf, ": ", sep="")
		cat( .g_print_help(damage, width, digits, exdent = 2, prefix = prf), sep="\n" )
	}
	
	#Sum of squares
	cat( "Sum of squares: ", object$squares, "\n", sep="" )
	
	# Loglikelihood
	cat( "Loglikelihood (ignoring multinomial_coefficient): ", object$LL, "\n", sep="" )
	
	#SPPE
	cat( "SPPE: ", object$SPPE, "\n", sep="" )
	
	#SVR
	cat( "SVR: ", object$SVR, "\n", sep="" )
	
	# Footer
	cat( "\n", sep="" )
}

# Do print.
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