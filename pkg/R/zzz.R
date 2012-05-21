.onLoad <- function(libname, pkgname){
  require("methods", character.only = TRUE, quietly = TRUE)  # check complains!
  loadRcppModules()
}
