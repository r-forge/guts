##
# Function sv_cor_prob_matrix().
# soeren.vogel@posteo.ch
# 2015-07-12
# Function produces a matrix of correlation coefficients (upper tri) and p-values (lower tri) together.
# Usage: sv_cor_prob_matrix(data.frame("x"=1:10, "y"=rnorm(10), "z"=runif(10)))
# License: GPL >= 2
#
sv_cor_prob_matrix <- function(X, df=nrow(X)-2) {
	stopifnot( is.data.frame(X) )
	M <- cor(X)
	below <- lower.tri(M) # row(M) > col(M)
	r2 <- M[below]^2
	Fstat <- r2 * df / (1 - r2)
	M[below] <- 1 - pf(Fstat, 1, df)
	return(M)
}
