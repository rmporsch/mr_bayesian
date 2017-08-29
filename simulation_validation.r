#'#'####################################
#'  Validation of Simulation with MR  #'
#'#####################################'

#devtools::install_github("MRCIEU/TwoSampleMR")

library(TwoSampleMR)
library(dplyr)
folder <- "simulation/"
files <- list.files(folder, 'GWAS')
out <- list()
for (i in 1:200) {

  gwas.files <- paste0(folder, 'simGWAS_', i, '.txt')
  solution.files <- paste0(folder, 'simSolution_', i, '.txt')

  sim <- read.table(gwas.files, head=T)
  outcome <- read.table(solution.files, head=T)

  outcome <- suppressMessages(suppressWarnings(test.mr(sim)))
  names(outcome)[1:2] <- c("from", "to")
  out[[i]] <- compare.outcomes(arcs, outcome)
}


out <- as.data.frame(do.call('rbind', lapply(out, check.mr.performance)))
names(out) <- c("notEstimated", "OverEstimated")
apply(out, 1, function(k) all(k==0))
mean(out[,2] != 0)
mean(out[,1] != 0)
