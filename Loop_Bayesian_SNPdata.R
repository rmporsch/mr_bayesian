library(Rgraphviz)
library(bnlearn)
library(compare)
library(plyr)
library(dplyr)
library(gdata)

Sim <- list()
Sim1 <- list()
df.true.model <- list()
df.causal.con <- list()
df.all.arcs <- list()
score_true_model <- list()
score_causal.con <- list()
score_all.con <- list()
Total_score <- list()
comparisson_1 <- list()
comparisson_2 <- list()
Score_compare_tm_sc <- list()
Score_compare_tm_ap <- list()
true_t1t2t3 <- list()
sc_t1t2t3 <- list()
sa_t1t2t3 <- list()
Score_true_t1t2t3 <- list()
Score_sc_t1t2t3 <- list()
Score_sa_t1t2t3 <- list()
comparisson_3 <- list()
comparisson_4 <- list()
Score_compare_tm_sc_t1t2t3 <- list()
Score_compare_tm_sa_t1t2t3 <- list()

Total_score <- list()


###SNPDATA

setwd(dir = "C:/Users/bbs520/Dropbox/Mendelian_Randomization/Bayesian Network/BayesianNetworkLearn/SIM/SIM4")
##### Preperate the data for analysis

for (i in 1:500){
  
  snp <- c("snp1.1" ,"snp1.2","snp1.3","snp1.4","snp1.5",
           "snp2.1","snp2.2","snp2.3","snp2.4","snp2.5",
           "snp3.1" ,"snp3.2","snp3.3","snp3.4","snp3.5")
  t1 <- c("t1")
  bl1 <- t(as.matrix(rbind(t1,snp)))
  colnames(bl1) <- NULL
  t2 <- c("t2")
  bl2 <- t(as.matrix(rbind(t2,snp)))
  colnames(bl2) <- NULL
  t3 <- c("t3")
  bl3 <- t(as.matrix(rbind(t3,snp)))
  colnames(bl3) <- NULL
  
  bl <- rbind(bl1,bl2,bl3)
  
  
  #perform the BNlearning algorithm to all simulated data
  temp <- read.delim(paste0('simData_',i, '.txt'), header=T, sep = " ")
  #temp <- temp[4:18] <- lapply(temp[4:18], as.numeric)
  temp$snp1.1 <- as.numeric(temp$snp1.1)
  temp$snp1.2 <- as.numeric(temp$snp1.2)
  temp$snp1.3 <- as.numeric(temp$snp1.3)
  temp$snp1.4 <- as.numeric(temp$snp1.4)
  temp$snp1.5 <- as.numeric(temp$snp1.5)
  
  temp$snp2.1 <- as.numeric(temp$snp2.1)
  temp$snp2.2 <- as.numeric(temp$snp2.2)
  temp$snp2.3 <- as.numeric(temp$snp2.3)
  temp$snp2.4 <- as.numeric(temp$snp2.4)
  temp$snp2.5 <- as.numeric(temp$snp2.5)
  
  temp$snp3.1 <- as.numeric(temp$snp3.1)
  temp$snp3.2 <- as.numeric(temp$snp3.2)
  temp$snp3.3 <- as.numeric(temp$snp3.3)
  temp$snp3.4 <- as.numeric(temp$snp3.4)
  temp$snp3.5 <- as.numeric(temp$snp3.5)
  
  Sim[[i]] <- hc(temp, blacklist = bl)

  
  #Select only the causal paths
  Sim1[[i]] <- lapply(Sim[[i]]$nodes, function(k) k$children)
  
  #Read in True data
  df.true.model[[i]] <- read.delim(paste0('simSolution_',i, '.txt'), header=T, sep = " ")
  df.true.model[[i]]$value <- NULL
  colnames(df.true.model[[i]])[1] <- "from"
  colnames(df.true.model[[i]])[2] <- "to"
  
  
  #Unlist causal arcs simulate data to data frame
  df.causal.con[[i]] <- data.frame(from = rep(names(Sim1[[i]]), sapply(Sim1[[i]], length)),
                                   to = unlist(Sim1[[i]]))
  
  df.all.arcs[[i]] <- as.data.frame(Sim[[i]]$arcs)
  
  score_true_model[[i]] <- nrow(df.true.model[[i]])
  score_causal.con[[i]] <- nrow(df.causal.con[[i]])
  score_all.con[[i]] <- nrow(df.all.arcs[[i]])
  comparisson_1[[i]] <- merge(df.true.model[[i]], df.causal.con[[i]])
  comparisson_2[[i]] <- merge(df.true.model[[i]], df.all.arcs[[i]])
  Score_compare_tm_sc[[i]] <- nrow(comparisson_1[[i]])
  Score_compare_tm_ap[[i]] <- nrow(comparisson_2[[i]])
  
  true_t1t2t3[[i]] <- df.true.model[[i]][grep("t", df.true.model[[i]]$from),]
  sc_t1t2t3[[i]] <- df.causal.con[[i]][grep("t", df.causal.con[[i]]$from),]
  sa_t1t2t3 [[i]] <- df.all.arcs[[i]][grep("t", df.all.arcs[[i]]$from),]
  
  Score_true_t1t2t3[[i]] <- nrow(df.true.model[[i]][grep("t", df.true.model[[i]]$from),])
  Score_sc_t1t2t3[[i]] <- nrow(df.causal.con[[i]][grep("t", df.causal.con[[i]]$from),])
  Score_sa_t1t2t3 [[i]] <- nrow(df.all.arcs[[i]][grep("t", df.all.arcs[[i]]$from),])
  
  comparisson_3[[i]] <- merge(true_t1t2t3[[i]], sc_t1t2t3[[i]])
  comparisson_4[[i]] <- merge(true_t1t2t3[[i]], sa_t1t2t3 [[i]])
  
  Score_compare_tm_sc_t1t2t3[[i]] <- nrow(comparisson_3[[i]])
  Score_compare_tm_sa_t1t2t3[[i]] <- nrow(comparisson_4[[i]])
  
  Total_score[[i]] <- cbind(Score_true_t1t2t3[[i]],Score_sc_t1t2t3[[i]],
                            Score_compare_tm_sc_t1t2t3[[i]])
  
  df.Total_score <- ldply(Total_score, data.frame)

  colnames(df.Total_score)[1] <- "True_t1t2t3"
  colnames(df.Total_score)[2] <- "Sim_causal_t1t2t3"
  colnames(df.Total_score)[3] <- "True_vs_Sim_Causal_t1t2t3"
 
  
}   

df.Total_score <-add_rownames(df.Total_score, "Model_nr")
head(df.Total_score)

## Underestimations

underestimations <- df.Total_score[which(df.Total_score$True_t1t2t3!=df.Total_score$True_vs_Sim_Causal_t1t2t3),]
underestimations1 <- as.data.frame(cbind(underestimations$Model_nr,underestimations$True_t1t2t3,underestimations$True_vs_Sim_Causal_t1t2t3))

colnames(underestimations1)[1] <- "model_nr"
colnames(underestimations1)[2] <- "True_model"
colnames(underestimations1)[3] <- "Sim_causal"


## Overestimation

overestimations <- df.Total_score[which(df.Total_score$Sim_causal_t1t2t3 > df.Total_score$True_t1t2t3),]
overestimations1 <- as.data.frame(cbind(overestimations$Model_nr,overestimations$True_t1t2t3,overestimations$Sim_causal_t1t2t3))

colnames(overestimations1)[1] <- "model_nr"
colnames(overestimations1)[2] <- "True_model"
colnames(overestimations1)[3] <- "Sim_causal"

## Correct

Correct <-  df.Total_score[which(df.Total_score$True_t1t2t3==df.Total_score$Sim_causal_t1t2t3),]
Correct1 <- as.data.frame(cbind(Correct$Model_nr,Correct$True_t1t2t3,Correct$Sim_causal_t1t2t3))


rm(list=setdiff(ls(), "Correct1","underestimations1", "overestimations1"))

keep(Correct1,underestimations1, overestimations1, sure=TRUE)

