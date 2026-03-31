## Package Loader
req_pkgs <- c(
  # Base math / matrix / sampling
  "Matrix", "MASS", "mvtnorm", "invgamma",
  # Bayesian / Polya-Gamma
  "BayesLogit", "coda", "Metrics",
  # Spatial data handling
  "sf", "spdep", "spatialreg",
  # Visualization
  "ggplot2", "RColorBrewer", "scales",
  # Data manipulation
  "tidyverse", "dplyr", "tidycensus", "foreign", "tidyr", "tigris",
  # Simulation / survey
  "sampling", "survey", "ltm", "mase", "srvyr",
  # Misc
  "rlang", "HDInterval", "readr"
)

for (pkg in req_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing missing package: ", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(Matrix)
  library(MASS)
  library(mvtnorm)
  library(invgamma)
  library(BayesLogit)
  library(coda)
  library(Metrics)
  library(sf)
  library(spdep)
  library(spatialreg)
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
  library(tidyverse)
  library(dplyr)
  library(tidycensus)
  library(foreign)
  library(tidyr)
  library(sampling)
  library(survey)
  library(ltm)
  library(mase)
  library(srvyr)
  library(rlang)
  library(HDInterval)
  library(readr)
  library(tigris)
})

message("All required packages are loaded.")

?everything
