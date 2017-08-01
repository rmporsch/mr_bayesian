#!/usr/bin/env Rscript
source("/home/tshmak/WORK/myRpackages/General/R/TimStartup.R")
source('functions.r')
library(dplyr)
TimStartup()
Tim.load(Rplink)
Tim.load(mysimtools) # Tim.load is basically just devtools::load <- all() but applied to my directories
useWTCCC() # This runs bfile('/path/to/WTCCC bfile')
# bfile() specifies the --bfile to use with plink
set.seed(42)
n  <- 2000; p <- 2000
pdfile <- subsample(n=n, p=p)
dd <- load.to.R(pdfile)

fam <- read.table(paste0(pdfile, ".fam"))
bim <- read.table(paste0(pdfile, ".bim"))

#'#'################
#'  effect model  #'
#'#################'

# proportion of SNPs with an effect
genetic.effectmat <- matrix(c(0.01, 0.000, 0.000,
                              0.000, 0.01, 0.000,
                              0.000, 0.000, 0.01), nrow=3, ncol=3, byrow=T)

effect <- c(rep(0, 6), runif(20,0,0.8))
library(tidyverse)
for (item in 1:1000) {
  pheno.effectmat <- generate.random.causes(effect)
  report.arcs <- data.frame(pheno.effectmat, 'to'=rownames(pheno.effectmat))
  report.arcs <- gather(report.arcs, 'from', 'value', -to) %>%
    filter(value!=0) %>%
    select(from, to, value)
  arcs <- make.arcs(pheno.effectmat)
  arcs <- as.data.frame(arcs)
  names(arcs) <- c('from', 'to')
  arcs <- merge(arcs, report.arcs, all.x=T )

  out <- simulate.complex(n, p, effectmat, dd, pheno.effectmat)
  write.table(out, paste0('simulation/simData_', item, '.txt'), row.names=F, quote=F)
  write.table(arcs, paste0('simulation/simSolution_', item, '.txt'), row.names=F, quote=F)
}
zip('simulations.zip', 'simulation')
