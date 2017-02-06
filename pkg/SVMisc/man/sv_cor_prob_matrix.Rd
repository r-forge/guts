\encoding{UTF-8}

\name{sv_cor_prob_matrix}
\alias{sv_cor_prob_matrix}

\title{Matrix of correlation coefficients and p-values}

\description{The function produces a matrix of correlation coefficients (upper tri) and p-values (lower tri) together.}

\usage{
sv_cor_prob_matrix(X, df = nrow(X) - 2)
}

\arguments{
	\item{X}{Data frame from which the correlations are calculated.}
	\item{df}{Integer. Degrees of freedom.}
}

%\details{}

\value{A matrix of correlation coefficients (upper tri) and p-values (lower tri).}

%\references{}

\author{\packageAuthor{SVMisc}}

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{SVMisc-package}}}

\examples{
sv_cor_prob_matrix(data.frame("x"=1:10, "y"=rnorm(10), "z"=runif(10)))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
