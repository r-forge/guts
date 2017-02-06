##
# Function sv_export_to_SPSS().
# soeren.vogel@posteo.ch
# 2012-10-12, 2015-07-09
# Export an R data frame to SPSS with labels, levels and attributes.
# License: GPL >= 2
#
sv_export_to_SPSS <- function( data, info, filebasename, labelcolname = "label", levelscolname = "levels", levelsseparator = "; ", stringsAsFactors = FALSE ) {

	# Error checks
	if ( !is.data.frame( data ) ) stop ( "data must be a data.frame!" )
	if ( !is.data.frame( info ) ) stop ( "info must be a data.frame!" )
	if ( ncol(data) != nrow(info) ) stop( "ncol(data) != nrow(info)" )
	# check whether cols exist

	# Levels
	takelevelsfromdata <- FALSE
	if ( is.logical( levelscolname ) ) {
		if ( !levelscolname ) {
			takelevelsfromdata <- TRUE
		}
		else {
			stop( "levelscolname can be either false or a character!" )
		}
	}
	else if ( is.character( levelscolname ) ) {
		if ( !exists( levelscolname, info ) )
			stop( "levelscolname does not exist in info!" )
		}
	else {
		stop( "levelscolname can be either false or a character!" )
	}

	# File names
	filecsv <- paste( filebasename, '.csv', sep = '', collapse = '' )
	filesps <- paste( filebasename, '.sps', sep = '', collapse = '' )
	filesav <- paste( filebasename, '.sav', sep = '', collapse = '' )
	filespv <- paste( filebasename, '.spv', sep = '', collapse = '' )
	cat( "", file = filecsv, append = FALSE )
	cat( "", file = filesps, append = FALSE )

	# Numeric version of the data frame
	data.num <- as.data.frame( sapply( data, sv_make_numeric ), stringsAsFactors = stringsAsFactors )
	write.csv( data.num, file = filecsv, row.names = FALSE, na = '' )
#	write.table( data.num, file = filetxt, append = FALSE, quote = TRUE, sep = ",", eol = ",", na = "\"\"", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = "double", fileEncoding = "UTF-8")

	# SPS syntax
	cat( "/* „–Ä–“: Syntax Script Auto Generated in UTF 8 „–Ü–“ */\n\n", file = filesps, sep = "", append = TRUE )
	
	cat( "SET UNICODE=YES.\n\n", file = filesps, sep = "", append = TRUE )
	
	# /ENCODING=\'UTF8\'
	# /FILE=\'', getwd(), '/', filecsv, '\'
	cat( 'GET DATA
	/TYPE=TXT
	/FILE=\'', filecsv, '\'
	/DELCASE=LINE
	/DELIMITERS=","
	/QUALIFIER=\'"\'
	/ARRANGEMENT=DELIMITED
	/FIRSTCASE=2
	/IMPORTCASE=ALL
	/VARIABLES=', file = filesps, sep = "", append = TRUE )
	for ( i in names(data) ) {
		## Format according to type
		
		## Logicals, (ordered) factors. No number of levels should exceed 999.
		if ( is.logical(data[, i]) | is.factor(data[, i]) ) {
			y <- 'F3.0'
		} else if ( is.numeric(data[, i]) ) {
			if ( sum(as.numeric(data[, i]), na.rm=TRUE) == 0 ) {
				# prevent integer overflow
				nc <- 0
			} else {
				nc <- max(nchar(as.numeric(na.exclude(round(data[, i])))))
			}
			
			xd <- na.exclude(data[, i] - trunc(data[, i]))
			if ( any(xd!=0) ) {
				nd <- max(nchar(as.character(xd)))-2 # remove leading digit and point, if present
			} else {
				nd <- 0
			}
			if ( nd > 16 ) nd <- 16
			if ( is.na(nd) ) nd <- 0
			
			nc <- nc + nd
			if ( nc==0 | is.na(nc) ) {
				# at least length 1
				nc <- 1
				nd <- 0
			}
			
			# Longer digits can only go to character
			if ( nc > 40 ) {
				y <- paste( 'A', nc, sep='' )
			} else {
				y <- paste( 'F', nc, '.', nd, sep='' )
			}
		} else {
			nc <- max(nchar(as.character(data[, i])))
			if ( nc==0 | is.na(nc) ) nc <- 1
			if ( nc>32000 ) nc <- 32000
			y <- paste( 'A', nc, sep='' )
		}
		
		cat( "  ", i, " ", y, "\n", file = filesps, sep = "", append = TRUE )
	}
	cat('.
CACHE.
EXECUTE.
DATASET NAME DataSet1 WINDOW=FRONT.', "\n\n", file = filesps, sep = "", append = TRUE )

	cat( "VARIABLE LABELS\n", file = filesps, sep = "", append = TRUE )
	for ( i in names(data) ) {
		cat( "  ", i, " \"", substr(info[i, labelcolname], 1, 100), "\"\n", file = filesps, sep = "", append = TRUE ) # Labels
	}
	cat( ".\nEXECUTE.\n\n", file = filesps, sep = "", append = TRUE )


	cat( "VALUE LABELS\n", file = filesps, sep = "", append = TRUE )
	options( useFancyQuotes = FALSE )
	for ( i in names(data) ) {
		if ( is.logical( data[, i] ) | is.factor( data[, i] ) ) {
			if ( is.logical( data[, i] ) ) {
				y <- '0 "FALSE" 1 "TRUE"'
			}
			else {
				if ( takelevelsfromdata ) {
					x <- levels( data[, i] )
				}
				else {
					x <- strsplit( info[i, levelscolname], levelsseparator, fixed = TRUE, perl = FALSE )[[1]]
				}
				x <- gsub( " *$", "", x )
				x <- gsub( "^ *", "", x )
				y <- paste( ' ', 1:length(x), dQuote(substr(x, 1, 100)), sep = ' ') # Levels
			}
			cat( "  /", i, " ", y, file = filesps, sep = '', append = TRUE, fill = 80 )
		}
	}
	options( useFancyQuotes = TRUE )
	cat( ".\nEXECUTE.\n\n", file = filesps, sep = '', append = TRUE )

	cat( "VARIABLE LEVEL\n", file = filesps, sep = '', append = TRUE )
	for ( i in names(data) ) {
		if ( is.logical(data[, i]) | is.factor(data[, i]) | is.numeric(data[, i]) ) {
			cat( "  ", i, "\n", file = filesps, sep = "", append = TRUE )
		}
	}
	cat( "(SCALE)\n.\nEXECUTE.\n\n", file = filesps, sep = '', append = TRUE )

	cat( 'SAVE OUTFILE=\'', filesav, '\'', "\n", '  /COMPRESSED', "\n.\n", file = filesps, sep = '', append = TRUE )
	cat( 'OUTPUT SAVE OUTFILE=\'', filespv, '\'', "\n.", file = filesps, sep = '', append = TRUE )
}
