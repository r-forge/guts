\encoding{UTF-8}

\name{GUTS-package}
\alias{GUTS-package}
%\alias{attr<-}
\alias{attr<-.GUTS}
%\alias{attributes<-}
\alias{attributes<-.GUTS}
%\alias{mostattributes<-}
\alias{mostattributes<-.GUTS}
\alias{print.GUTS}
\alias{modguts}
\alias{[[<-.GUTS}
\alias{$<-.GUTS}

\docType{package}

\title{Fast Calculation of the Likelihood of a Stochastic Survival Model}

\description{GUTS (General Unified Theory of Survival) is a stochastic survival model for ecotoxicology. The package allows for the definition of exposure and survival time series as well as parameter values, and the fast calculation of the survival probabilities as well as the logarithm of the corresponding likelihood.}

\details{
\tabular{ll}{
Package: \tab GUTS\cr
Type: \tab Package\cr
%Version: \tab 1.0\cr
%Date: \tab 2015-05-03\cr
License: \tab GPL (>= 2)\cr
}
A GUTS object is a special list of class \dQuote{GUTS}. Functions \code{\link{guts_setup}}, \code{\link{guts_calc_loglikelihood}} and \code{\link{guts_calc_survivalprobs}} are available to create and work with GUTS objects. A data set \link{diazinon} is also included. See links for more details.
}

\author{Carlo Albert \email{carlo.albert@eawag.ch} and Sören Vogel \email{soeren.vogel@posteo.ch}

Maintainer: Sören Vogel \email{soeren.vogel@posteo.ch}}

\note{%
\itemize{
	\item Caution: This is a complete reimplementation of the former package. Please update your code.%
	\item The package is under current development. Parts of the software are likely to change with future releases. Refer to the manual pages after updating your copy. Save your work before using this software.%
	\item Extended help and an in-depth documentation can be found at the project web site \href{http://guts.r-forge.r-project.org}{http://guts.r-forge.r-project.org}.%
}
} % End of \note.

\references{Jager, T., Albert, C., Preuss T. and Ashauer R. (2011) General unified theory of survival -- a toxicokinetic toxicodynamic framework for ecotoxicology. Published online in \emph{Env. Sci. and Tech.} 
(\href{http://www.mendeley.com/research/general-unified-threshold-model-survival-toxicokinetictoxicodynamic-framework-ecotoxicology/}{http://www.mendeley.com/research/general-unified-threshold-model-survival-toxicokinetictoxicodynamic-framework-ecotoxicology/}).}

\seealso{\code{\link{guts_setup}}, \code{\link{guts_calc_loglikelihood}}, \code{\link{guts_calc_survivalprobs}}, \code{\link{diazinon}}, \code{\link{Rcpp}}}
