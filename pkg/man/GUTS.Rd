\encoding{UTF-8}

\name{GUTS}

\alias{GUTS}
\alias{guts}
\alias{guts_setup}
\alias{guts_calc_loglikelihood}
\alias{guts_calc_survivalprobs}



\title{Fast Calculation of the Likelihood of a Stochastic Survival Model}



\description{GUTS (General Unified Theory of Survival) is a stochastic survival model for ecotoxicology. Use \code{guts_setup} to define exposure and survival time series as well as parameter values. Use \code{guts_calc_loglikelihood} to calculate the survival probabilities and the corresponding loglikelihood. \code{guts_calc_survivalprobs} is a convenient wrapper for calculate survival probabilities.}



\usage{
guts_setup(C, Ct, y, yt, par = 0, dist = "lognormal",
	model = "Proper", N = 1000, M = 10000)

guts_calc_loglikelihood(gobj, pars)

guts_calc_survivalprobs(gobj, yts = 0)
}



\arguments{%
	\item{C}{Numeric vector of concentrations. Vector must contain at least 2 values and be of the same length as Ct.%
	}
	\item{Ct}{Numeric vector of concentration time points. Vector must contain at leas 2 values and be of the same length as C. Numerics must start at 0, and contain unique values in ascending order. The last value of \code{Ct} must not be smaller than the last value of \code{yt}, i.e., there must be information about concentrations when no survivor remains.%
	}
	\item{y}{Integer vector (counts) of survivors. Vector must contain at least 2 values and be of the same length as yt. y must not be ascending.%
	}
	\item{yt}{Numeric vector of survivor time points. Vector must contain at least 2 values and be of the same length as y. Numerics must start at 0, and contain unique values in ascending order. The last value of \code{yt} must not be larger than the last value of \code{Ct}, i.e., there must be information about concentrations when no survivor remains.%
	}
	\item{par}{Numeric vector of parameters. If \code{0}, all parameters are set to 0. See details below.%
	}
	\item{dist}{Distribution as character, either \dQuote{lornormal} (default) or \dQuote{delta}.%
	}
	\item{model}{Model as character, either \dQuote{Proper} (default) or \dQuote{IT}.%
	}
	\item{N}{Integer. Sample length. Must be greater than 2.%
	}
	\item{M}{Integer. Time grid points. Must be greater than 1.%
	}
	\item{gobj}{GUTS object. The object to be updated (and used for the calculation).%
	}
	\item{pars}{Numeric vector of new parameters. See details below.%
	}
	\item{yts}{Either a single positive integer or a vector of survivor time points. If \code{0}, survivor information in \code{gobj} are used.%
	}
} % End of \arguments



\details{%

\subsection{Field and Attribute Access}{%
Fields and attributes of an object of class \dQuote{GUTS} are read-only. To prevent accidental change of fields or attributes, replacement functions have be rewritten throwing an error. Use function \code{guts_setup} to modify fields.

However, parameters can be updated through the use of the function \code{guts_calc_loglikelihood}, and survivors and survivor time points are set when using the function \code{guts_calc_survivalprobs}.
} % End of \subsection{Field and Attribute Access}.

\subsection{Parameters}{%
For distribution \dQuote{lognormal} and model \dQuote{Proper}, 5 parameters are required:
\itemize{%
	\item hb: background mortality rate
	\item kr: recovery rate
	\item kk: killing rate
	\item m:  sample mean
	\item sd: sample standard deviation
}

The number of parameters is checked according to \code{dist} and \code{type}.
For type \dQuote{Proper}…
\itemize{%
	\item …if distribution is \dQuote{lognormal}, 5 non-negative parameters are required. If \code{sd} is not 0, \code{m} must not be 0.%
	\item …if distribution is \dQuote{delta}, 4 non-negative parameters are required.%
}
If dist is \dQuote{delta}, \code{sd} is set internally and must not be provided. If model is \dQuote{IT}, the number of required parameters is reduced by 1; the parameter \dQuote{kk} is set internally and must not be provided.

Wrong number of parameters will give an error. Wrong parameter values (e.g., negative values) give a warning to the object, setting the loglikelihood to \code{-Inf}.
} % End of \subsection{Parameters}.

\subsection{Distributions}{%
Samples are created internally according to \code{dist}, \code{model}, \code{N} (and \code{pars}).

\itemize{%
	\item If \code{dist} is set to \dQuote{lognormal} (default), a sample of size \code{N} is drawn from a lognormal distribution with \code{m, sd}. Note that \code{m} and \code{sd} are parameters for the \bold{lognormal} distribution, not its logarithm (differs from \code{rlnorm} in \R)!%
	\item If \code{dist} is set to \dQuote{delta}, a sample vector of size \code{N} from the delta distribution is created (with one parameter).%
	\item Model type \dQuote{SD} is implemented through distribution \dQuote{delta}.
}
} % End of \subsection{Distributions}.

} % End of \details



\value{
\code{guts_setup} returns a list of class \dQuote{GUTS} with the following fields:

\item{C}{Concentrations.}
\item{Ct}{Concentration time points.}
\item{y}{Survivors.}
\item{yt}{Survivor time points.}
\item{par}{Parameters.}
\item{dist}{Distribution.}
\item{model}{Model.}
\item{N}{Sample length.}
\item{M}{Time grid points.}
\item{results}{List of results.}

The list of results contain the following fields:

\item{S}{Vector of Survivor probabilities.}
\item{LL}{The loglikelihood.}

The object has the following attributes (for internal use):

\item{class}{\code{GUTS}}
\item{experiment}{Internal representation of distribution and model.}
\item{wpar}{Internal representation of parameters.}
\item{par_pos}{Internal representation of parameter positions.}
\item{checked}{Force creation of GUTS objects through function use.}

\code{guts_calc_loglikelihood} and \code{guts_calc_survivalprobs} update objects directly and do not return a value. Use the fields \code{object$results$LL} and \code{object$results$S} to access the loglikelihood and survival probabilities.
} % End of \value.



\references{Jager, T., Albert, C., Preuss T. and Ashauer R. (2011) General unified theory of survival -- a toxicokinetic toxicodynamic framework for ecotoxicology. Published online in \emph{Env. Sci. and Tech.} 
(\href{http://www.mendeley.com/research/general-unified-threshold-model-survival-toxicokinetictoxicodynamic-framework-ecotoxicology/}{http://www.mendeley.com/research/general-unified-threshold-model-survival-toxicokinetictoxicodynamic-framework-ecotoxicology/}).}



\author{Carlo Albert \email{carlo.albert@eawag.ch} and Sören Vogel \email{soeren.vogel@posteo.ch}

Maintainer: Sören Vogel \email{soeren.vogel@posteo.ch}}



\note{%
\itemize{
	\item Caution: This is a complete reimplementation of the former package. Please update your code.%
	\item The package is under current development. Parts of the software are likely to change with future releases. Refer to this manual page after updating your copy. \bold{Save your work before using this software.}%
	\item Extended help and documentation}{Extended help and an in-depth documentation can be found at the project web site \href{http://guts.r-forge.r-project.org}{http://guts.r-forge.r-project.org}.%
}
} % End of \note.



\seealso{\code{\link{diazinon}}, \code{\link[adaptMCMC]{adaptMCMC}}, \code{\link{GUTS-package}}}



\examples{
data(diazinon)

gts <- guts_setup(C = diazinon$C1, Ct = diazinon$Ct1,
	y = diazinon$y1, yt = diazinon$yt1)
guts_calc_loglikelihood(gts, c(0.05084761, 0.12641525,
	1.61840054, 19.09911, 6.495246))
gts
# LL: -183.566.

guts_calc_survivalprobs( gts, 5 )
gts
# LL: 0, S: 1, 0.9293491, 0.8409267, 0.7987192, 0.4578539.

guts_calc_survivalprobs( gts, 0:6 )
gts
# LL: 0, S: 1, 0.9293336, 0.8409121, 0.7987059,
#           0.4577675, 0.3383611, 0.3213131.

\dontrun{gts[["C"]] <- 1:3} # Error.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line




















