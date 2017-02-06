\encoding{UTF-8}

\name{sv_get_levels}
\alias{sv_get_levels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sv_get_levels(varname, codebook)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varname}{
%%     ~~Describe \code{varname} here~~
}
  \item{codebook}{
%%     ~~Describe \code{codebook} here~~
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
function (varname, codebook) 
{
    levs <- unlist(strsplit(codebook[varname, "Levels"], "; "))
    for (i in 1:length(levs)) {
        if (levs[i] == "") {
            levs[i] <- i
        }
    }
    return(levs)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
