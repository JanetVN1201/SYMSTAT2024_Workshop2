##Libraries necessary for the course

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("graph", "Rgraphviz"), dep=TRUE)

#Install INLA
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)

#Other libraries
libs = c("DiagrammeR","car","ggpubr","spdep","RColorBrewer",
         "spatstat","sp","latticeExtra", "viridis", "cleangeo",
         "gridExtra","gstat","raster","ggplot2","ggfortify",
         "survival","joineR","BayesSurvival","icenReg","nloptr",
         "faraway","lme4","boot","sf","coda","spBayesSurv",
         "BayesX", "R2BayesX", "fields", "R.rsp", "devtools",
         "rnaturalearth", "leaflet", "remotes", "hydroGOF",
         "INLAspacetime", "ggthemes", "ggthemes","splancs",
         "fmesher")

ix <- which(!sapply(libs, require, char = TRUE))
if (length(ix) > 0) {install.packages(libs[ix], repos = "https://cloud.r-project.org/")
  sapply(libs[ix], require, char = TRUE)}

remotes::install_github("inlabru-org/inlabru", ref = "devel")
