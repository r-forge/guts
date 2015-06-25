\encoding{UTF-8}

\name{GUTS}

\alias{GUTS}
\alias{guts}
\alias{guts_setup}
\alias{guts_calc_loglikelihood}
\alias{guts_calc_survivalprobs}



\title{Fast Calculation of the Likelihood of a Stochastic Survival Model}



\description{GUTS (General Unified Threshold model of Survival) is a stochastic survival model for ecotoxicology.  The package allows for the definition of exposure and survival time series as well as parameter values, and the fast calculation of the survival probabilities as well as the logarithm of the corresponding likelihood.}



\usage{
guts_setup(C, Ct, y, yt, dist = "lognormal",
	model = "Proper", N = 1000, M = 10000)

guts_calc_loglikelihood(gobj, par)

guts_calc_survivalprobs(gobj, par)
}



\arguments{%
	\item{C}{Numeric vector of concentrations.  Vector must contain at least 2 values and be of the same length as \code{Ct}.%
	}
	\item{Ct}{Numeric vector of concentration time points.  Vector must contain at least 2 values and be of the same length as \code{C}.  Numerics must start at 0, and contain unique values in ascending order.%
	}
	\item{y}{Integer vector (counts) of survivors.  Vector must contain at least 2 values and be of the same length as \code{yt}.  \code{y} must not be ascending.%
	}
	\item{yt}{Numeric vector of survivor time points.  Vector must contain at least 2 values and be of the same length as \code{y}.  Numerics must start at 0, and contain unique values in ascending order.  Survivor information at time points later than the latest concentration time point will be disregarded (with a warning).%
	}
	\item{dist}{Distribution as character, either \dQuote{lornormal} (default) or \dQuote{delta}.%
	}
	\item{model}{Model as character, either \dQuote{Proper} (for full model, the default) or \dQuote{IT} (for individual tolerance).%
	}
	\item{N}{Integer.  Thresholds sample length.  Must be greater than 2.%
	}
	\item{M}{Integer.  Number of time grid points.  Must be greater than 1.%
	}
	\item{gobj}{GUTS object.  The object to be updated (and used for the calculation).%
	}
	\item{par}{Numeric vector of parameters.  See details below.%
	}
} % End of \arguments



\details{%


\subsection{Functions}{%
Use \code{guts_setup} to define (or alter) a GUTS object.  Various checks are applied to the data.  On success, a GUTS object will be created.

Use \code{guts_calc_loglikelihood} to calculate the survival probabilities and the corresponding loglikelihood for a given set of parameters.  The function is very fast and can be used in routines for parameter estimation.  The function returns the loglikelihood, however it also updates the fields \code{par}, \code{S} and \code{LL} of the object.

\code{guts_calc_survivalprobs} is a convenience wrapper (calling \code{cuts_calc_loglikelihood}) that can be used for predictions; it returns the survival probabilities.
}


\subsection{Parameters, Models, and Distributions}{%

GUTS uses the following parameters:
\itemize{%
	\item hb: background mortality rate
	\item ke: dominant rate constant% was kr: recovery rate
	\item kk: killing rate
	\item mn: mean of thresholds z from distribution \code{dist}
	\item sd: standard deviation of thresholds z from distribution \code{dist}
}

For model \dQuote{Proper}, the first 3 parameters are needed as well as the parameters for the distribution (\code{mn} and \code{sd} for \dQuote{lognormal}, or \code{mn} for \dQuote{delta}).  For model \dQuote{IT} (individual tolerance), the third parameter (\code{kk}) is set internally to infinity and must not be provided.  For model type \dQuote{SD} (stochastic death), use distribution \dQuote{delta} and model \dQuote{Proper} with the first 4 parameters.

The number of parameters is checked according to \code{dist} and \code{model}.  Wrong number of parameters invoke an error, wrong parameter values (e.g., negative values) invoke a warning, and the loglikelihood is set to -Inf.

Thresholds are created internally according to \code{dist}, \code{model}, \code{N} (and \code{par}).  If \code{dist} is \dQuote{lognormal} (default), \code{N} ordered thresholds are created from the lognormal distribution with \code{mn} and \code{sd}.  Note that \code{mn} and \code{sd} are parameters for the \bold{lognormal} distribution, not its logarithm (differs from \code{rlnorm} in \R)!  If \code{dist} is \dQuote{delta}, all \code{N} thresholds are equal to \code{mn}.


} % End of \subsection{Parameters, Models, and Distributions}.


\subsection{Field and Attribute Access}{%
Fields and attributes of an object of class \dQuote{GUTS} are read-only.  To prevent accidental change of fields or attributes, replacement functions were rewritten throwing an error when used.  Always use function \code{guts_setup} to create objects or modify fields on existing objects.  However, \code{guts_calc_loglikelihood} and \code{guts_calc_survivalprobs} update an object’s fields \code{par} (parameters), \code{S} (survival probabilities) and \code{LL} (the loglikelihood).
} % End of \subsection{Field and Attribute Access}.
} % End of \details



\value{
\code{guts_setup} returns a list of class \dQuote{GUTS} with the following fields:

\item{C}{Concentrations.}
\item{Ct}{Concentration time points.}
\item{y}{Survivors.}
\item{yt}{Survivor time points.}
\item{dist}{Distribution.}
\item{model}{Model.}
\item{N}{Sample length.}
\item{M}{Time grid points.}
\item{par}{Parameters.}
\item{S}{Vector of survivor probabilities.}
\item{LL}{The loglikelihood.}

%The object has the following attributes (for internal use):
%
%\item{class}{\code{GUTS}}
%\item{experiment}{Internal representation of distribution and model.}
%\item{wpar}{Internal representation of parameters.}
%\item{par_pos}{Internal representation of parameter positions.}

\code{guts_calc_loglikelihood} returns the loglikelihood.

\code{guts_calc_survivalprobs} returns the survival probabilities.
} % End of \value.



\note{%
The GUTS project web site can be found here: \href{http://guts.r-forge.r-project.org}{http://guts.r-forge.r-project.org}.
} % End of \note.



\references{Jager, T., Albert, C., Preuss T. and Ashauer R. (2011) General unified threshold model of survival -- a toxicokinetic toxicodynamic framework for ecotoxicology. Environmental Science \& Technology, 45(7), 2529--2540, doi: 10.1021/es103092a.}



\author{Carlo Albert \email{carlo.albert@eawag.ch} and Sören Vogel \email{soeren.vogel@posteo.ch}

Maintainer: Sören Vogel \email{soeren.vogel@posteo.ch}}



\seealso{\code{\link{diazinon}}, \code{\link[adaptMCMC]{adaptMCMC}}, \code{\link{GUTS-package}}}



\examples{
data(diazinon)

gts <- guts_setup(C = diazinon$C1, Ct = diazinon$Ct1,
	y = diazinon$y1, yt = diazinon$yt1)

guts_calc_loglikelihood(gts, c(0.05084761, 0.12641525,
	1.61840054, 19.09911, 6.495246)) # -183.566.
guts_calc_survivalprobs(gts, rep(.5, 5)) # 1.000000e+00 1.737352e-05 ... 1.674316e-33

guts_calc_loglikelihood(gts, -1:3) # Warning.
gts

\dontrun{guts_calc_survivalprobs( gts, 1:4 )} # Error.

\dontrun{gts[["C"]] <- 1:3} # Error.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line




















