\encoding{UTF-8}

\name{sv_make_numeric}
\alias{sv_make_numeric}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sv_make_numeric(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (x) 
{
    y <- NULL
    var.class <- class(x)[1]
    if (var.class == "numeric") {
        y <- x
    }
    else if (var.class == "logical") {
        y <- as.numeric(x)
    }
    else if (var.class == "ordered") {
        y <- as.numeric(x)
    }
    else if (var.class == "factor") {
        y <- as.numeric(x)
    }
    else {
        y <- as.character(x)
        y <- gsub("\n", "<br>", y, ignore.case = FALSE, perl = FALSE, 
            fixed = FALSE, useBytes = FALSE)
        y <- gsub("\t", "<tab>", y, ignore.case = FALSE, perl = FALSE, 
            fixed = FALSE, useBytes = FALSE)
    }
    return(y)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
