# .onLoad <- function(lib, pkg) {
#  #require(methods)
#   cat("nothing to do here...\n")
#   cat("please run topOnto::intiONT(ontology_name) to choose ontology.")
# }

.onAttach <- function(lib, pkg) {
  # where <- match(paste("package:", pkg, sep=""), search())
  # initVar(where) 
  ## some preprocessing
  init()
}
