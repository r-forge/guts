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
The class \code{Rcpp_GUTS} is the reference class for the C++ class \code{GUTS} exposed to R through the \link{Rcpp} module \code{modguts}.
The class allows for the definition of exposure and survival time series as well as parameter values, and the calculation of the survival probabilities and the logarithm of the corresponding likelihood.
\code{GUTS} is the constructor object of \code{Rcpp_GUTS}.%
}

\section{Objects from the Class}{%
New objects can be created by calls to the method \code{new()} (see \link{ReferenceClasses}), or the function \code{new("Rcpp_GUTS")}.
No arguments must be provided, and a “factory fresh” \code{GUTS} object contains default (empty) values.
Fields are read only. Setter methods must be used to set values on fields.%
}

%\section{Slots}{The class has no slots. Private fields can be set using methods (see below).}

\section{Methods}{
  \describe{
    \item{new}{\code{signature()}: Create a new \code{GUTS} object with non-sense values.}

    \item{setConcentrations}{\code{signature(C = "double", Ct = "double")}:
      Set the vectors of concentrations (\code{C}) and concentration time points (\code{Ct}).
      \code{C} and \code{Ct} must have the same length of at least 2 values each.
      \code{Ct} must start at 0, and contain unique values in ascending order.
      The last value of \code{Ct} must not be smaller than the last value of \code{yt},
        i.e., there must be information about concentrations when no survivor remain.
    }

    \item{setSurvivors}{\code{signature(y = "integer", yt = "double")}:
      Set the vectors of survivors (\code{y}) and survivor time points (\code{yt}).
      \code{y} must be a vector of integers (counts).
      \code{y} and \code{yt} must have the same length of at least 2 values each.
      \code{yt} must start at 0, and contain unique values in ascending order.
      The last value of \code{yt} must not be larger than the last value of \code{Ct},
        i.e., there must be information about concentrations when no survivor remain.
    }

    \item{setParameters}{\code{signature(par = "double")}:
      Set parameters (\code{length(par) > 2}) for the experiment.
      The first three parameters (background mortality rate, recovery rate, and killing rate) are used
        for the calculation of the logarithm of the likelihood.
      Two more parameters may be set and are used to construct a sample of length \code{N}
        (see \code{setSampleLength()}) from distribution \code{dist} (see \code{setDistribution()}).
      If \code{dist} is set to “empirical”, the two parameters are ignored.
      If \code{dist} is set to “lognormal” (the default), the two parameters specify
        the \emph{mean} and \emph{standard deviation} of the \emph{lognormal} distribution.
      Note that this differs from \code{\link{rlnorm}} where the parameters denote mean and standard deviation
        of the \emph{corresponding normal} distribution.)
    }

    \item{setTimeGridPoints}{\code{signature(M = "integer")}:
      Set the number of grid points on the time axis (for numerical exactness).
      \code{M} must be an integer greater than 1.
    }
    
    \item{setDistribution}{\code{signature(dist = "character")}:
      Set the \bold{name} of the distribution to create (currently either \dQuote{lognormal} – the default – or \dQuote{empirical}).
      See section \bold{Distributions} for details.
      Do not append distribution parameters here!
      Method \code{setSample} calls this method to set \code{dist} to \dQuote{empirical}.
    }

    \item{setSampleLength}{\code{signature(N = "integer")}:
      Set the length of the sample that is drawn from the distribution (numerical exactness).
      \code{N} must be an integer greater than 1.
      \code{setSample} calls this method to set \code{N} to the length of the vector passed to \code{setSample}.
    }

    \item{setSample}{\code{signature(z = "double")}:
      Set an empirical sample \bold{instead} of a theoretical distribution to sample from.
      \code{z} must contain at least 2 non-negative values.
      \code{z} can be unordered, it is ordered internally where the original numeric vector remains as set.
      \code{setSample} calls \code{setDistribution} to set \code{dist} to \dQuote{empirical},
        and calls \code{setSampleLength} to set \code{N} to the length of \code{z}.
    }

    \item{calcSample}{\code{signature()}:
      Compute a sample from \code{par[4:5]}, \code{N}, and \code{dist}.
      Either \code{calcSurvivalProbabilities} or \code{calcLoglikelihood} will call
        \code{calcSample} to set the sample vector \code{z}.
      However, \code{calcSample} can be used directly to calculate set a new sample vector.
    }

    \item{calcSurvivalProbabilities}{\code{signature(...)}:
      Calculate the vector of damages (\code{D}) and the survival probabilities (\code{S}).
      \code{calcLoglikelihood} calls \code{calcSurvivalProbabilities} to set the vector of survival probabilities \code{S}.
      However, \code{calcSurvivalProbabilities} can be used directly to calculate and set a vectors \code{D} and \code{S}.
      This method is overloaded, which means it may be called with various arguments.
        See section \bold{Calculate Survival Probabilities} for details.
    }

    \item{calcLoglikelihood}{\code{signature()}:
      Calculate the logarithm of the likelihood of a properly set-up \code{"GUTS"} object.
      Sets the logliklihood (\code{LL}), if calculated, else to \code{NA}.
    }

    \item{show}{\code{signature(object = "Rcpp_GUTS")}:
      Display current values of a GUTS object. S4 wrapper for field access with a formatted output.
    }
  }
}

\section{Fields}{%
To access the survival probabilities or the loglikelihood, use the fields \code{S} and \code{LL}.

\code{C}, \code{Ct}, \code{y}, \code{yt}, \code{par}, \code{M}, \code{dist}, \code{N}, \code{z}, \code{D}, \code{S}, \code{LL}, and \code{errors} are fields of a GUTS object. Note that any of these fields is \bold{read-only}. Content can be modified using setters and methods (see above). However, field contents can be displayed, or assigned to other \R variables using the standard \R procedures.%
}

\section{Distributions}{%
  \code{GUTS} creates samples internally according to \code{dist}, \code{N}, and \code{pars[4:length(pars)]}.
  If no values for the creation of the distribution are given, no sample is created.
  If \code{dist} is set to \code{"empirical"}, \code{setSample} must be invoked with a valid sample vector.
  If \code{dist} is set to \code{"lognormal"}, a sample of size \code{N} is drawn from a lognormal distribution
    with \code{mean = pars[4]} and \code{sigma = pars[5]}. Note that mean and sigma are parameters for the
    \bold{lognormal} distribution, not its logarithm (differs from \code{rlnorm} in \R)!
}

\section{Method: Calculate Survival Probabilities}{%
  This method is overloaded. It can be called in either of three ways:
  \describe{
    \item{with no argument}{%
      It checks for all fields properly set up, except the survivors field \code{y}.%
    }
    \item{\code{signature(Stlength = "integer")}}{%
      \code{Stlength} must be a single integer.
      The method creates a vector of zeros of length Stlength, and a vector \code{0:Stlength},
        and passes both vectors to \code{setSurvivors()}.
        Note that this will render the vector of survivors \code{y} \emph{invalid} for the
          calculation of the loglikelihood!%
    }
    \item{\code{signature(St = "double")}}{%
      \code{St} must be a numeric vector.
      The method creates a vector of zeros of length of St,
        and passes both vectors to \code{setSurvivors()}.
        Note that this will render the vector of survivors \code{y} \emph{invalid} for the
          calculation of the loglikelihood!%
    }
  }
}

\section{Functions}{
  \describe{
    \item{\code{print(object)}}{%
      Simple wrapper for \code{show}.%
    }
    \item{\code{logLik(object, ...)}}{%
      Wrapper for \code{calcLoglikelihood()}.
      Returns an S3 object of class “logLik” with the loglikelihood as value and the length of the survivor’s vector (\code{y})
        as the degrees of freedom.
      Additional arguments (\code{...}) may be supplied but are ignored.%
    }
  }
}

\section{Notes and warnings}{
  \describe{
%    \item{Compiling from source}{Compilation from source would require an installation of the C++ boost library (\href{http://www.boost.org/}{http://www.boost.org/}). The library is expected in \code{/opt/local/include} under UNIX-likes, and in \code{/Boost} under Windows. Change the makefiles in \code{GUTS/src} to fit your installation.}
    \item{Copy objects}{%
      New objects \bold{must} be created using the \code{new{}} method (see above).
      Assigning \R objects to other \R objects would \bold{not} create a \bold{copy of the object}
        but a \bold{copy of the reference} (pointing to the same object).
      See examples below for details.%
    }
    \item{Serialisation}{%
      Serialisation is currently not implemented.
      Objects of class \code{GUTS} cannot be saved across sessions.
      Clean out every \code{GUTS} object from your workspace before leaving \R.
      On restart, rebuild \code{GUTS} objects if necessary.%
    }
    \item{Development status}{%
      The package is under current development.
      Party of the software are likely to change with future releases. Refer to this manual page after updating your copy.
      \bold{Save your work before using this software.}%
    }
  }
}

\references{Jager T., Albert C., Preuss T., Ashauer R., General unified theory of survival – a toxicokinetic toxicodynamic framework for ecotoxicology, published online in Env. Sci. and Tech. (2011) (\href{http://www.mendeley.com/research/general-unified-threshold-model-survival-toxicokinetictoxicodynamic-framework-ecotoxicology/\#}{http://www.mendeley.com/research/general-unified-threshold-model-survival-toxicokinetictoxicodynamic-framework-ecotoxicology/#}).}

\author{%
Carlo Albert \email{carlo.albert@eawag.ch} & Sören Vogel \email{soeren.vogel@uzh.ch}

Maintainer: Sören Vogel \email{soeren.vogel@uzh.ch}
}

\seealso{\code{\link{ReferenceClasses}}, \code{\link{Rcpp}}, \code{\link{new}}, \code{\link{diazinon}}, \code{\link{logLik}}, \code{\link{print}}, \code{\link[adaptMCMC]{adaptMCMC}}}

\examples{
library("GUTS")

## documentation
\dontrun{vignette("GUTS-Documentation")}

## what GUTS is
\dontrun{getRefClass("Rcpp_GUTS")}

## defaults
obj0 <- new("Rcpp_GUTS") # also: obj0 <- GUTS$new()
obj0$LL # NaN
logLik(obj0)

## one experiment
data(diazinon)

obj0$setConcentrations(
   C = diazinon$C1,
  Ct = diazinon$Ct1
)
obj0$setSurvivors(
   y = diazinon$y1,
  yt = diazinon$yt1
)
obj0$setParameters(
  par = c(0.05084761, 0.12641525, 1.61840054, 19.09911, 6.495246)
)
obj0$setTimeGridPoints(M = 10000)
obj0$setDistribution(dist = "lognormal")
obj0$setSampleLength(N = 10000)
obj0$calcLoglikelihood()
obj0$LL # abt -183
logLik(obj0) # abt -183
obj0$S # vector of survival probabilities

## R: use meanlog, meansd; GUTS: use mean, sigma
zR <- rlnorm(10000, 2.8943280, 0.33420543)
obj1 <- GUTS$new()
obj1$setParameters(c(rep(1, 3),
  mean(zR), sd(zR))) # mean and sigma of lognormal!
obj1$setSampleLength(10000)
obj1$setDistribution("lognormal")
obj1$calcSample()

## compare:
\dontrun{zDiff <- sort(zR) - sort(obj1$z)}
\dontrun{par(mfrow=c(3, 1))}
\dontrun{plot(density(obj1$z))}
\dontrun{plot(density(zR))}
\dontrun{plot(density(zDiff))}

## Copy copies the reference, not the object
obj2 <- obj1
obj2$setSampleLength(N = 1999) # sets N also on obj1!
obj1$N # 1999
}

\keyword{classes}
%\keyword{GUTS-class}
