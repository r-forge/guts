\encoding{UTF-8}

\name{sv_make_factor}
\alias{sv_make_factor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sv_make_factor(varname, numerics, codebook = NULL, add = 0, ordered = FALSE, levels = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varname}{
%%     ~~Describe \code{varname} here~~
}
  \item{numerics}{
%%     ~~Describe \code{numerics} here~~
}
  \item{codebook}{
%%     ~~Describe \code{codebook} here~~
}
  \item{add}{
%%     ~~Describe \code{add} here~~
}
  \item{ordered}{
%%     ~~Describe \code{ordered} here~~
}
  \item{levels}{
%%     ~~Describe \code{levels} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}

\author{\packageAuthor{SVMisc}}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (varname, numerics, codebook = NULL, add = 0, ordered = FALSE, 
    levels = NULL) 
{
    if (!is.null(codebook)) {
        levs <- sv_get_levels(varname, codebook)
    }
    else if (!is.null(levels)) {
        levs <- levels
    }
    else {
        levs <- unique(as.character(sort(numerics)))
    }
    nlevs <- length(levs)
    vals <- c(1:length(levs), numerics + add)
    fact <- factor(vals, ordered = ordered, labels = levs)
    return(fact[(nlevs + 1):(length(fact))])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
