\encoding{UTF-8}
\name{Rcpp_GUTS-class}
\Rdversion{1.1}
\docType{class}

\alias{GUTS}
\alias{GUTS-class}
\alias{GUTS-method}
\alias{GUTS-package}

\alias{Rcpp_GUTS}
\alias{Rcpp_GUTS-class}
\alias{Rcpp_GUTS-method}

\alias{modguts}

\alias{print.Rcpp_GUTS}
\alias{logLik.Rcpp_GUTS}



\title{Class \code{Rcpp_GUTS}}

\description{%
\code{GUTS} (General Unified Theory of Survival) is a stochastic survival model for ecotoxicology.
%The class \code{Rcpp_GUTS} is the reference class for the C++ class \code{GUTS} exposed to \R through the \link{Rcpp} module \code{modguts}.
The class allows for the definition of exposure and survival time series as well as parameter values, and the calculation of the survival probabilities and the logarithm of the corresponding likelihood.
%\code{GUTS} is the constructor object of \code{Rcpp_GUTS}.%
}



\section{Objects from the Class}{%
New objects can be created by calls to the method \code{new()} (see \link{ReferenceClasses}), or the function \code{new("Rcpp_GUTS")}.
A \dQuote{factory fresh} \code{GUTS} object contains default (or empty) values.
Fields are read only.
Setter methods must be used to set values on fields.
} % End of \section{Objects from the Class}.



%\section{Slots}{The class has no slots. Private fields can be set using methods (see below).}



\section{Methods}{
\describe{
	\item{new}{\code{signature()}: Create a new \code{GUTS} object.}

	\item{setConcentrations}{\code{signature(C="double", Ct="double")}:
		Set the vectors of concentrations (\code{C}) and concentration time points (\code{Ct}).
		\code{C} and \code{Ct} must be unnamed numerics and of the same length of at least 2 values each.
		\code{Ct} must start at 0, and contain unique values in ascending order.
		The last value of \code{Ct} must not be smaller than the last value of \code{yt}, i.e., there
			must be information about concentrations when no survivor remains.
	}

	\item{setSurvivors}{\code{signature(y="integer", yt="double")}:
		Set the vectors of survivors (\code{y}) and survivor time points (\code{yt}).
		\code{y} and \code{yt} must be unnamed numerics and of the same length of at least 2 values each.
		\code{y} must be a vector of integers (counts).
		\code{yt} must start at 0, and contain unique values in ascending order.
		The last value of \code{yt} must not be larger than the last value of \code{Ct}, i.e., there
			must be information about concentrations when no survivor remains.
	}

	\item{setDistribution}{\code{signature(dist="character")}:
		Set the name of the distribution to create (currently either \dQuote{lognormal} (default),
			\dQuote{delta} or \dQuote{empirical}).
		Distribution name is case-insensitive. See details below.
		Do not append distribution parameters here!
	}

	\item{setModel}{\code{signature(type="character")}:
		Set the model type (currently either \dQuote{Proper} (default) or \dQuote{IT}).
		Model name is case-insensitive.
		Model \dQuote{SD} is implemented through distribution \dQuote{delta}.
		See details below.
	}

	\item{setParameters}{\code{signature(par="double")}:
		Set parameters.
		\code{par} must be unnamed numerics.
		See details below.
	}

	\item{setTimeGridPoints}{\code{signature(M="integer")}:
		Set the number of grid points on the time axis (for numerical exactness).
		\code{M} must be an unnamed integer greater than 1.
	}

	\item{setSampleLength}{\code{signature(N="integer")}:
		Set the length of the sample that is drawn from the distribution (numerical exactness).
		\code{N} must be an unnamed integer greater than 1.
		\code{setSample} calls this method to set \code{N} to the length of the vector passed to \code{setSample}, see below.
	}

	\item{calcLoglikelihood}{\code{signature()}:
		Calculate the logarithm of the likelihood of a properly set-up \code{GUTS} object.
		Sets the logliklihood field (\code{LL}), if calculated.
	}

	\item{show}{\code{signature(object="Rcpp_GUTS")}:
		Display current values of a GUTS object.
		S4 wrapper for field access with a formatted output.
	}

	%\item{setInterval}{\code{signature(Interval="double")}: Experimental method. Do not use.%
	%}
}
} % End of \section{Methods}.





\section{Optional Methods}{
\describe{
	\item{calcSurvivalProbabilities}{\code{signature(...)}:
		Calculate the survival probabilities (\code{S}).
		\code{calcLoglikelihood} calls \code{calcSurvivalProbabilities} to set the vector of survival probabilities \code{S}.
		However, \code{calcSurvivalProbabilities} can be used directly to calculate vector \code{S}.
		See details below.
	}
	\item{setSample}{\code{signature(z="double")}:
		Set an empirical sample instead of a theoretical distribution to sample from.
		\code{z} must contain at least 2 non-negative unnamed numerics.
		\code{z} can be unordered.
		\code{dist} and \code{N} are adjusted accordingly.
		The sample values are not stored in \code{GUTS}, and hence cannot be accessed from a \code{GUTS} object.
	}
}
} % End of \section{Optional Methods}.





\section{Fields}{%
To access the survival probabilities or the loglikelihood, use the fields \code{S} and \code{LL}.

\code{C}, \code{Ct}, \code{y}, \code{yt}, \code{dist}, \code{model}, \code{par}, \code{M}, \code{N},
\code{S}, \code{LL} and \code{Messages} are read-only fields of a GUTS object.
} % End of \section{Fields}.





\section{Details}{%



\subsection{Distributions}{%
\code{GUTS} creates samples internally according to \code{dist}, \code{pars} and \code{N}.

\itemize{%
	\item If \code{dist} is set to \code{"lognormal"} (default), a sample of size \code{N} is drawn from a
			lognormal distribution with \code{m, sd}.
		Note that m and sd are parameters for the \bold{lognormal} distribution, not its logarithm (differs from \code{rlnorm} in \R)!
	\item If \code{dist} is set to \code{"delta"}, a sample vector of size \code{N} from the
			delta distribution is created (with one parameter).
		Model type \dQuote{SD} is implemented through distribution \dQuote{delta}.
	\item If \code{dist} is set to \code{"empirical"}, \code{setSample} must be invoked with a valid sample vector.
}
} % End of \subsection{Distributions}.



\subsection{Parameters}{%
For distribution \dQuote{lognormal} and model \dQuote{Proper}, \code{GUTS} requires 5 parameters:
\itemize{%
	\item hb: background mortality rate
	\item kr: recovery rate
	\item kk: killing rate
	\item m:  sample mean
	\item sd: sample standard deviation
}

\code{GUTS} checks the number of parameters according to \code{dist} and \code{type}.
For type \dQuote{Proper}…
\itemize{%
	\item …if distribution is \dQuote{lognormal}, 5 non-negative parameters are required
			and \code{m} must be greater than 0.
	\item …if distribution is \dQuote{delta}, 4 non-negative parameters are required.
	\item …if distribution is \dQuote{empirical}, 3 non-negative parameters are required.
}
If dist is \dQuote{delta}, the standard deviation is set internally and must not be provided. If model is \dQuote{IT}, the number of required parameters is reduced by 1; the parameter \dQuote{kk} is set internally and must not be provided.

\code{m} can be 0 if and only if sd is 0.

Wrong number of parameters will set an error to the object, setting the loglikelihood to NaN.
Wrong parameter values (e.g., negative) set a warning to the object, setting the loglikelihood to -Inf.
} % End of \subsection{Parameters}.



\subsection{Method Calculate Survival Probabilities}{%
This can be called in either of three ways:
\describe{%
	\item{with no argument}{It checks for all fields properly set up, except the survivors field \code{y}.}

	\item{\code{signature(Stlength="integer")}}{\code{Stlength} must be a single unnamed integer.
		The method creates a vector of zeros of length \code{Stlength} for \code{y}, and
			a vector \code{0:Stlength} for \code{yt}, and passes both to \code{setSurvivors()}.
		Note that this will render the vector of survivors \code{y} \emph{invalid} for the calculation of the loglikelihood!}

	\item{\code{signature(St="double")}}{\code{St} must be an unnamed numeric vector of survivor time points.
		The method creates a vector of zeros of length of St, and passes both vectors to \code{setSurvivors()}.
		Note that this will render the vector of survivors \code{y} \emph{invalid} for the calculation of the loglikelihood!}
}
} % End of \subsection{Method Calculate Survival Probabilities}.



\subsection{Errors, Messages}{%
\code{GUTS} objects store errors and messages internally (but do not throw them to the \R console).
Messages can be displayed by printing \code{GUTS} objects.
} % End of \subsection{Errors, Messages}.
} % End of \section{Details}.





\section{Functions}{%
\describe{%
	\item{\code{print(object)}}{Simple wrapper for \code{show}.}%
	\item{\code{logLik(object, ...)}}{Wrapper for \code{calcLoglikelihood()}.
		Returns an S3 object of class \dQuote{logLik} with the loglikelihood as value and the
			length of the survivor vector (\code{y}) as the degrees of freedom.
		Additional arguments (\code{...}) may be supplied but are ignored.}
}
} % End of \section{Functions}





\section{Notes and Warnings}{%
\describe{%
	\item{Copy objects}{New objects \bold{must} be created using the \code{new()} method (see above).
		Assigning \R objects to other \R objects would \bold{not} create a \bold{copy of the object}
			but a \bold{copy of the reference} (pointing to the same object).
		See examples below for details.}

	\item{Serialisation}{Serialisation is currently not implemented.
		Objects of class \code{GUTS} cannot be saved across sessions.
		Clean out every \code{GUTS} object from your workspace before leaving \R.
		On restart, rebuild \code{GUTS} objects if necessary.}

	\item{Development status}{The package is under current development.
		Parts of the software are likely to change with future releases.
		Refer to this manual page after updating your copy.
		\bold{Save your work before using this software.}}

	\item{Extended help and documentation}{Extended help and an in-depth documentation can be found
			at the project web site
			\href{http://guts.r-forge.r-project.org}{http://guts.r-forge.r-project.org}.}
	}
} % \section{Notes and Warnings}.





\references{Jager, T., Albert, C., Preuss T. and Ashauer R. (2011) General unified theory of survival -- a toxicokinetic toxicodynamic framework for ecotoxicology. Published online in \emph{Env. Sci. and Tech.} 
(\href{http://www.mendeley.com/research/general-unified-threshold-model-survival-toxicokinetictoxicodynamic-framework-ecotoxicology/}{http://www.mendeley.com/research/general-unified-threshold-model-survival-toxicokinetictoxicodynamic-framework-ecotoxicology/}).}

\author{Carlo Albert \email{carlo.albert@eawag.ch} and Sören Vogel \email{soeren.vogel@posteo.ch}

Maintainer: Sören Vogel \email{soeren.vogel@posteo.ch}}

\seealso{\code{\link{ReferenceClasses}}, \code{\link{Rcpp}}, \code{\link{new}}, \code{\link{diazinon}}, \code{\link{logLik}}, \code{\link{print}}, \code{\link[adaptMCMC]{adaptMCMC}}}





\examples{
library("GUTS")

## what GUTS is
\dontrun{getRefClass("Rcpp_GUTS")}

## defaults
obj0 <- new("Rcpp_GUTS") # also: obj0 <- GUTS$new()
obj0$LL # NaN
logLik( obj0 )

## one experiment
data(diazinon)
obj0$setConcentrations(C = diazinon$C1, Ct = diazinon$Ct1)
obj0$setSurvivors(y = diazinon$y1, yt = diazinon$yt1)
obj0$setDistribution(dist = "lognormal")
obj0$setModel(model = "Proper")
obj0$setParameters(par = c(0.05084761, 0.12641525, 1.61840054, 19.09911, 6.495246))
obj0$setTimeGridPoints(M = 10000)
obj0$setSampleLength(N = 1000)
obj0$calcLoglikelihood()
obj0$LL        # abt -183
logLik(obj0)   # abt -183, df 23
obj0$S         # vector of survival probabilities


## Print:
obj0

## Copy copies the reference, not the object
obj1 <- GUTS$new()
obj1$setSampleLength( n=2050 )
obj2 <- obj1
obj2$setSampleLength( N = 1999 ) # sets N also on obj1!
obj1$N # 1999
}

\keyword{classes}
%\keyword{GUTS-class}
