##
# Function sv_cs_nagel().
# soeren.vogel@posteo.ch
# 2014-07-06
# Compute Nagelkerke’s R, Cox-&-Snell R from a fitted binary logistic regression model.
# Example:
#   y <- sample(c(T, F), 50, repl=T)
#   x <- sample(1:7, 50, repl=T)
#   mod <- glm(y ~ x, family=binomial("logit"))
#   sv_cs_nagel(mod)
# License: GPL >= 2
# 
sv_cs_nagel <- function( mod ) {
	llnull <- mod$null.deviance
	llmod <- mod$deviance
	n <- nrow(mod$model) # n <- length(mod$fitted.values)
	Rcs <- 1 - exp( (llmod - llnull) / n )
	Rnagel <- Rcs / (1 - exp(-llnull/n))
	out <- list( 'Rcs'=Rcs, 'Rnagel'=Rnagel, 'N'=n )
	class(out) <- c('list', 'table')
	return(out)
}
