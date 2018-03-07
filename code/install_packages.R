# load packages and install missing packages

# devtools
if (!require(devtools)) {
  install.packages("devtools")
}

# install greta package and all dependencies
if (!require(greta)) {
  devtools::install_github("greta-dev/greta", ref = "master")
  library(greta)
}

# install modified version of gretaDynamics from jdyen repo
if (!require(gretaDynamics)) {
  devtools::install_github("jdyen/gretaDynamics", ref = "master")
  library(gretaDynamics)
}


