loadModule("modguts", TRUE)

# this has changed since R 2.15.1 with Rcpp 0.9.13
#.onLoad <- function(libname, pkgname){
#  require("methods", character.only = TRUE, quietly = TRUE)  # check complains!
#  loadRcppModules()
#}
