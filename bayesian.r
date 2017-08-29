source('functions.r')
library(bnlearn)
folder <- 'simulation/'

out <- list()
for (i in 1:200) {
  temp <- read.table(paste0(folder, 'simData_', i, '.txt'), head=T)
  truemodel <- read.delim(paste0(folder, 'simSolution_',i, '.txt'), header=T, sep = " ")[,1:2]
  out[[i]] <- bayesian_estiamtion(temp, truemodel)
}

results <- do.call('rbind', out)


mean(results$correct / results$truemodeltotal, na.rm=T)
mean(results$OverEstimation >0)
mean(results$Correct == results$TrueModelTotal & (results$OverEstimation == 0)) 

