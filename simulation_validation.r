#'#'####################################
#'  Validation of Simulation with MR  #'
#'#####################################'

#devtools::install_github("MRCIEU/TwoSampleMR")

library(TwoSampleMR)
library(dplyr)
source('functions.r')
folder <- "simulation/"
files <- list.files(folder, 'GWAS')
out <- list()
for (i in 1:200) {

  gwas.files <- paste0(folder, 'simGWAS_', i, '.txt')
  solution.files <- paste0(folder, 'simSolution_', i, '.txt')

  sim <- read.table(gwas.files, head=T)
  real <- read.table(solution.files, head=T)
  real <- real[!is.na(real$value), 1:2]


  outcome <- suppressMessages(suppressWarnings(test.mr(sim)))
  names(outcome)[1:2] <- c("from", "to")
  out[[i]] <- compare.outcomes(real, outcome[,1:2])
}

undebug(compare.outcomes)


out <- as.data.frame(do.call('rbind', lapply(out, check.mr.performance)))
names(out) <- c("notEstimated", "OverEstimated")
apply(out, 1, function(k) all(k==0))
mean(out[,2] != 0)
mean(out[,1] != 0)

u  <- rnorm(100, 0, 1) + rnorm(100,0,1)
var(u)
