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
n  <- 5000; p <- 200
pdfile <- subsample(n=n, p=p)
dd <- load.to.R(pdfile)

fam <- read.table(paste0(pdfile, ".fam"))
bim <- read.table(paste0(pdfile, ".bim"))

#'#'################
#'  effect model  #'
#'#################'

# proportion of SNPs with an effect
genetic.effectmat <- matrix(c(0.010, 0.000, 0.000,
                              0.000, 0.01, 0.000,
                              0.000, 0.000, 0.01), nrow=3, ncol=3, byrow=T)
rownames(genetic.effectmat) <- paste0('t',1:3)
colnames(genetic.effectmat) <- paste0('t',1:3)

effect <- c(rep(0, 6), runif(20,0,0.8))
effect <- c(rep(0, 6), rep(0.2, 2))
library(tidyverse)
#for (item in 1:100) {
#  pheno.effectmat <- generate.random.causes(effect)
#  report.arcs <- data.frame(pheno.effectmat, 'to'=rownames(pheno.effectmat))
#  report.arcs <- gather(report.arcs, 'from', 'value', -to) %>%
#    filter(value!=0) %>%
#    select(from, to, value)
#  arcs <- make.arcs(pheno.effectmat)
#  arcs <- as.data.frame(arcs)
#  names(arcs) <- c('from', 'to')
#  arcs <- merge(arcs, report.arcs, all.x=T )
#
#  out <- simulate.causal.prs(n, p, effectmat, dd, pheno.effectmat)
#  write.table(out, paste0('simulation/simData_', item, '.txt'), row.names=F, quote=F)
#  write.table(arcs, paste0('simulation/simSolution_', item, '.txt'), row.names=F, quote=F)
#}
#zip('simulations_prs.zip', 'simulation')
#system('rm -r simulation/')

dir.create('simulation')
for (item in 1:200) {
  pheno.effectmat <- generate.random.causes(effect)
  report.arcs <- data.frame(pheno.effectmat, 'to'=rownames(pheno.effectmat))
  report.arcs <- gather(report.arcs, 'from', 'value', -to) %>%
    filter(value!=0) %>%
    select(from, to, value)
  arcs <- make.arcs(pheno.effectmat)
  arcs <- as.data.frame(arcs)
  names(arcs) <- c('from', 'to')
  arcs <- merge(arcs, report.arcs, all.x=T )

  out <- simulate.causal.snp(n, p, genetic.effectmat, dd, pheno.effectmat, 0.2, bim, dogwas=T)
  write.table(out[[1]], paste0('simulation/simData_', item, '.txt'), row.names=F, quote=F)
  write.table(out[[2]], paste0('simulation/simGWAS_', item, '.txt'), row.names=F, quote=F)
  write.table(arcs, paste0('simulation/simSolution_', item, '.txt'), row.names=F, quote=F)
}
#zip('simulations_snp.zip', 'simulation')
#system('rm -r simulation/')
