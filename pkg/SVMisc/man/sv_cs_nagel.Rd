\encoding{UTF-8}

\name{sv_cs_nagel}
\alias{sv_cs_nagel}

\title{Nagelkerke’s R and Cox & Snell R}

\description{Compute Nagelkerke’s R, Cox-&-Snell R from a fitted binary logistic regression model.}

\usage{
sv_cs_nagel(mod)
}

\arguments{
	\item{mod}{A fitted binary logistic regression model.}
}

%\details{}

\value{A list/table with Nagelkerke’s and Cox & Snell R.}

%\references{}

\author{\packageAuthor{SVMisc}}

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{SVMisc-package}}}

\examples{
y <- sample(c(T, F), 50, repl=T)
x <- sample(1:7, 50, repl=T)
mod <- glm(y ~ x, family=binomial("logit"))
sv_cs_nagel(mod)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }

