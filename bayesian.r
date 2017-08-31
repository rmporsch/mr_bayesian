source('functions.r')
library(bnlearn)
library(dplyr)
folder <- 'simulation/'

i <- 1
out <- list()
for (i in 1:200) {
  temp <- read.table(paste0(folder, 'simData_', i, '.txt'), head=T)
  truemodel <- read.delim(paste0(folder, 'simSolution_',i, '.txt'), header=T, sep = " ")[,1:2]
  out[[i]] <- bayesian_estiamtion(temp, truemodel)
}

results <- do.call('rbind', out)


head(results)
results %>%
group_by(ModelNumber) %>%
summarise("Identified Connections"= mean(Correct == TrueModelTotal, na.rm=T),
          "OverEstimation" = mean(OverEstimation >0),
          "AllCorrect" = mean(Correct == TrueModelTotal & (OverEstimation == 0)),
          "Number of Simulations" = n()) 


combinations <- expand.grid(0:1, 0:1, 0:1)
combinations <- as.matrix(combinations)
out <- list()
mat <- matrix(NA, nrow=3, ncol=3)
mat[upper.tri(mat, diag=F)] <- 0

for (item in 1:nrow(combinations)) {
  temp <- mat
  temp[lower.tri(temp, diag=F)] <- as.vector(combinations[item,])
  out[[item]] <- temp 
}

