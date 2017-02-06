\encoding{UTF-8}

\name{sv_export_to_SPSS}
\alias{sv_export_to_SPSS}

\title{Export Data to IBM SPSS}

\description{Export R data to CSV and SAV, for import in SPSS.}

\usage{
sv_export_to_SPSS(data, info, filebasename, labelcolname = "label", levelscolname = "levels", levelsseparator = "; ", stringsAsFactors = FALSE)
}

\arguments{
	\item{data}{Data frame to export.}
	\item{info}{Data frame with variable names, labels, and levels. Will be used to generate labels in SPSS.}
	\item{filebasename}{Name of the file to export to, without extension.}
	\item{labelcolname}{Name of the column in \code{info} that holds the labels of the variables.}
	\item{levelscolname}{Name of the column in \code{info} that holds the levels of the variables.}
	\item{levelsseparator}{Separator used in the levels’ column to separat the levels.}
	\item{stringsAsFactors}{Should strings be exported as factors.}
}

%\details{}

\value{Writes the CSV- and the SAV-files to the file system.}

%\references{}

\author{\packageAuthor{SVMisc}}

%\note{}

% \section{Warning}{}

\seealso{\code{\link{SVMisc-package}}}

\examples{
## Currently no example.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }
