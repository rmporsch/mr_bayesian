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

pheno.effectmat <- matrix(c(NA, 0.000, 0.000,
                            0.7, NA, 0.000,
                            0.0, 0.7, NA), nrow=3, ncol=3, byrow=T)


out <- simulate.complex(n, p, effectmat, dd, 0.2)

cor(out)
write.table(out, 'fixed_simulation_large_effect.txt', row.names=F, quote=F)

out <- list()
for (item in 1:20) {
  out[[item]] <- simulate.complex(n, p, effectmat, dd, 0.2)
}

out.cor <- lapply(out, cor)
apply(simplify2array(out.cor), 1:2, mean)
apply(simplify2array(out.cor), 1:2, sd)
