library(Rgraphviz)
library(bnlearn)
library(compare)
library(plyr)
library(dplyr)

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

setwd(dir = "C:/Users/bbs520/Dropbox/Mendelian_Randomization/Bayesian Network/BayesianNetworkLearn/SIM/SIM3")



for (i in 1:500){
  bl1 <- matrix(c("t1","t1","t1","t1","t1","t1","t1","t1","t1","t1","t1","t1","t1","t1","t1", "snp1.1" ,"snp1.2","snp1.3","snp1.4","snp1.5","snp2.1" ,"snp2.2","snp2.3","snp2.4","snp2.5","snp3.1" ,"snp3.2","snp3.3","snp3.4","snp3.5"),ncol=2)
  bl2 <- matrix(c("t2","t2","t2","t2","t2","t2","t2","t2","t2","t2","t2","t2","t2","t2","t2", "snp1.1" ,"snp1.2","snp1.3","snp1.4","snp1.5","snp2.1" ,"snp2.2","snp2.3","snp2.4","snp2.5","snp3.1" ,"snp3.2","snp3.3","snp3.4","snp3.5"),ncol=2)
  bl3 <- matrix(c("t3","t3","t3","t3","t3","t3","t3","t3","t3","t3","t3","t3","t3","t3","t3", "snp1.1" ,"snp1.2","snp1.3","snp1.4","snp1.5","snp2.1" ,"snp2.2","snp2.3","snp2.4","snp2.5","snp3.1" ,"snp3.2","snp3.3","snp3.4","snp3.5"),ncol=2)
  bl <- rbind(bl1,bl2,bl3)
  
  #perform the BNlearning algorithm to all simulated data
  temp <- read.delim(paste0('simData_',i, '.txt'), header=T, sep = " ")
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
  
  Total_score[[i]] <- cbind(score_true_model[[i]],score_causal.con[[i]],score_all.con[[i]],Score_compare_tm_sc[[i]],
                            Score_compare_tm_ap[[i]],Score_true_t1t2t3[[i]],Score_sc_t1t2t3[[i]],Score_sa_t1t2t3 [[i]],
                            Score_compare_tm_sc_t1t2t3[[i]],Score_compare_tm_sa_t1t2t3[[i]])
  
  df.Total_score <- ldply(Total_score, data.frame)
  colnames(df.Total_score)[1] <- "True_Model"
  colnames(df.Total_score)[2] <- "Sim_causal_path"
  colnames(df.Total_score)[3] <- "Sim_all_paths"
  colnames(df.Total_score)[4] <- "True_vs_Sim_Causal"
  colnames(df.Total_score)[5] <- "True_vs_Sim_all"
  colnames(df.Total_score)[6] <- "True_t1t2t3"
  colnames(df.Total_score)[7] <- "Sim_causal_t1t2t3"
  colnames(df.Total_score)[8] <- "Sim_all_t1t2t3"
  colnames(df.Total_score)[9] <- "True_vs_Sim_Causal_t1t2t3"
  colnames(df.Total_score)[10] <- "True_vs_Sim_all_t1t2t3"
}   

print(sum(df.Total_score$True_vs_Sim_Causal)) / (sum(df.Total_score$True_Model)) *100
print(sum(df.Total_score$True_vs_Sim_all)) / (sum(df.Total_score$True_Model)) *100
print(sum(df.Total_score$True_vs_Sim_Causal_t1t2t3)) / (sum(df.Total_score$True_t1t2t3)) *100
print(sum(df.Total_score$True_vs_Sim_all_t1t2t3)) / (sum(df.Total_score$True_t1t2t3)) *100

df.Total_score <- add_rownames(df.Total_score, "Model_nr")
head(df.Total_score)

#mis identified missing arcs

missing_arcs <- df.Total_score[which(df.Total_score$True_t1t2t3!=df.Total_score$True_vs_Sim_Causal_t1t2t3),]
head(missing_arcs)
missing_arcs1 <- as.data.frame(cbind(missing_arcs$Model_nr,missing_arcs$True_t1t2t3,missing_arcs$True_vs_Sim_Causal_t1t2t3))
head(missing_arcs1)

colnames(missing_arcs1)[1] <- "model_nr"
colnames(missing_arcs1)[2] <- "True_model"
colnames(missing_arcs1)[3] <- "Sim_causal"

overestimations <- df.Total_score[which(df.Total_score$Sim_causal_t1t2t3 > df.Total_score$True_t1t2t3),]
overestimations1 <- as.data.frame(cbind(overestimations$Model_nr,overestimations$True_t1t2t3,overestimations$Sim_causal_t1t2t3))

colnames(overestimations1)[1] <- "model_nr"
colnames(overestimations1)[2] <- "True_model"
colnames(overestimations1)[3] <- "Sim_causal"
head(overestimations1)

true_t1t2t3[[3]]
sc_t1t2t3[[3]]



head(missing_arcs1)


head(df.Total_score)


head(df.Total_score)








a <- read.delim("simData_1.txt", header=T, sep =" ")

head(a)


a$snp1.1 <- as.numeric(a$snp1.1)
a$snp1.2 <- as.numeric(a$snp1.2)
a$snp1.3 <- as.numeric(a$snp1.3)
a$snp1.4 <- as.numeric(a$snp1.4)
a$snp1.5 <- as.numeric(a$snp1.5)

a$snp2.1 <- as.numeric(a$snp2.1)
a$snp2.2 <- as.numeric(a$snp2.2)
a$snp2.3 <- as.numeric(a$snp2.3)
a$snp2.4 <- as.numeric(a$snp2.4)
a$snp2.5 <- as.numeric(a$snp2.5)

a$snp3.1 <- as.numeric(a$snp3.1)
a$snp3.2 <- as.numeric(a$snp3.2)
a$snp3.3 <- as.numeric(a$snp3.3)
a$snp3.4 <- as.numeric(a$snp3.4)
a$snp3.5 <- as.numeric(a$snp3.5)

b <- hc(a)

b$arcs

















for (i in 1:100){
  bl1 <- matrix(c("t1","t1","t1", "snp1" ,"snp2","snp3"),ncol=2)
  bl2 <- matrix(c("t2","t2","t2", "snp1" ,"snp2","snp3"),ncol=2)
  bl3 <- matrix(c("t3","t3","t3", "snp1" ,"snp2","snp3"),ncol=2)
  bl <- rbind(bl1,bl2,bl3)
  
  #perform the BNlearning algorithm to all simulated data
  temp <- read.delim(paste0('simData_',i, '.txt'), header=T, sep = " ")
  temp$snp1 <- as.factor(temp$snp1)
  temp$snp2 <- as.factor(temp$snp2)
  temp$snp3 <- as.factor(temp$snp3)
  Sim[[i]] <- iamb(temp, test="cor", blacklist = bl)
  
  #Select only the causal paths
  Sim1[[i]] <- lapply(Sim[[i]]$nodes, function(k) k$children)
  
  #Read in True data
  df.true.model[[i]] <- read.delim(paste0('simSolution_',i, '.txt'), header=T, sep = " ")
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
  
  Total_score[[i]] <- cbind(score_true_model[[i]],score_causal.con[[i]],score_all.con[[i]],Score_compare_tm_sc[[i]],
                            Score_compare_tm_ap[[i]],Score_true_t1t2t3[[i]],Score_sc_t1t2t3[[i]],Score_sa_t1t2t3 [[i]],
                            Score_compare_tm_sc_t1t2t3[[i]],Score_compare_tm_sa_t1t2t3[[i]])
  
  df.Total_score <- ldply(Total_score, data.frame)
  colnames(df.Total_score)[1] <- "True_Model"
  colnames(df.Total_score)[2] <- "Sim_causal_path"
  colnames(df.Total_score)[3] <- "Sim_all_paths"
  colnames(df.Total_score)[4] <- "True_vs_Sim_Causal"
  colnames(df.Total_score)[5] <- "True_vs_Sim_all"
  colnames(df.Total_score)[6] <- "True_t1t2t3"
  colnames(df.Total_score)[7] <- "Sim_causal_t1t2t3"
  colnames(df.Total_score)[8] <- "Sim_all_t1t2t3"
  colnames(df.Total_score)[9] <- "True_vs_Sim_Causal_t1t2t3"
  colnames(df.Total_score)[10] <- "True_vs_Sim_all_t1t2t3"
}   

print(sum(df.Total_score$True_vs_Sim_Causal)) / (sum(df.Total_score$True_Model)) *100
print(sum(df.Total_score$True_vs_Sim_all)) / (sum(df.Total_score$True_Model)) *100
print(sum(df.Total_score$True_vs_Sim_Causal_t1t2t3)) / (sum(df.Total_score$True_t1t2t3)) *100
print(sum(df.Total_score$True_vs_Sim_all_t1t2t3)) / (sum(df.Total_score$True_t1t2t3)) *100


head(df.Total_score)


bl1 <- matrix(c("t1","t1","t1", "snp1" ,"snp2","snp3"),ncol=2)
bl2 <- matrix(c("t2","t2","t2", "snp1" ,"snp2","snp3"),ncol=2)
bl3 <- matrix(c("t3","t3","t3", "snp1" ,"snp2","snp3"),ncol=2)
bl <- rbind(bl1,bl2,bl3)

#perform the BNlearning algorithm to all simulated data
a <- read.delim(paste0('simData_',1, '.txt'), header=T, sep = " ")
a
a$snp1 <- as.numeric(temp$snp1)
a$snp2 <- as.numeric(temp$snp2)
a$snp3 <- as.numeric(temp$snp3)
Sim <- iamb(a, test="cor", blacklist = bl)

graphviz.plot(Sim)

b <- read.delim("simSolution_2.txt", header = T , sep = " ")
head(b)
colnames(b)


head(df.true.model)
head(df.causal.con)
















df.Total_score

TM <- sum(df.Total_score$`True Model`)
TM
SC <- sum(df.Total_score$`True vs Sim Causal`)
SC
SA <-sum(df.Total_score$`True vs Sim all`)
SA
38/53



#unlist elements lists into one dataframe


df
df.true.model[[3]]
df.all.arcs[[3]]

merge(df.true.model[[3]], df.all.arcs[[3]])





 ## order the data
  
  df.true.model <- df.true.model[order(df.true.model$from),]
  df.causal.con <- df.causal.con[order(df.causal.con$from),]
  df.all.arcs <- df.all.arcs[order(df.all.arcs$from),]
}  

  ## Score the different models
  
  score.truemodel <-nrow(model.testing$true.causal.model)
  score.causal.sim <-nrow(model.testing$sim.causal.con)
  score.all.con.sim <-nrow(model.testing$sim.all.con)
  
  score.truemodel <- as.data.frame(score.truemodel)
  score.causal.sim <- as.data.frame(score.causal.sim)
  score.all.con.sim <- as.data.frame(score.all.con.sim)
  

  # compare the simulated models with the true model
  
  compare_tm_to_sim.causal <- merge(df.true.model,df.causal.con)
  compare_tm_to_sim_allpaths <- merge(df.true.model,df.all.arcs)
  

  #score the comparison
  
  score_comp_tm_sc <- nrow(compare_tm_to_sim.causal)
  score_comp_tm_sa <- nrow(compare_tm_to_sim_allpaths)
  
  score_comp_tm_sc <- as.data.frame(score_comp_tm_sc)
  score_comp_tm_sa <- as.data.frame(score_comp_tm_sa)
  

  #cbind the different scores
  
  Total_score <- cbind(score.truemodel,score.causal.sim,score.all.con.sim, score_comp_tm_sc,score_comp_tm_sa)
  
  }



df

df <- data.frame(from = rep(names(Sim1[[1]]), sapply(Sim1[[1]], length)),
                            to = unlist(Sim1[[1]]))


df

t <- unlist(Sim[[1]])

t

list.causal.con = lapply(Sim1$nodes, function(k) k$children)




Sim1 <- iamb(Simdata1, test="cor", blacklist = bl)
graphviz.plot(Sim1, main = "SIM 1")
class(Sim1$arcs)

Sim1$nodes$prs2$mb[Sim1$nodes$prs2$mb%in%Sim1$nodes$prs2$children]
Sim1$nodes$prs2$mb[!Sim1$nodes$prs2$mb%in%Sim1$nodes$prs2$children]

df.true.model_1 <- read.delim("simSolution_1.txt", sep = " ")
df.true.model_2 <- read.delim("simSolution_2.txt", sep = " ")
df.true.model_3 <- read.delim("simSolution_3.txt", sep = " ")
df.true.model_4 <- read.delim("simSolution_4.txt", sep = " ")
df.true.model_5 <- read.delim("simSolution_5.txt", sep = " ")
df.true.model_6 <- read.delim("simSolution_6.txt", sep = " ")
df.true.model_7 <- read.delim("simSolution_7.txt", sep = " ")
df.true.model_8 <- read.delim("simSolution_8.txt", sep = " ")
df.true.model_9 <- read.delim("simSolution_9.txt", sep = " ")
df.true.model_10 <- read.delim("simSolution_10.txt", sep = " ")

########################################################################
########################################################################

Simdata1 <- read.delim("simData_1.txt", header=T, sep = " ")
Sim1 <- iamb(Simdata1, test="cor", blacklist = bl)
graphviz.plot(Sim1, main = "SIM 1")
class(Sim1$arcs)


#Select the identified causal pathways
list.causal.con = lapply(Sim1$nodes, function(k) k$children)

list.causal.con

list.causal.con
df.causal.con <- data.frame(from = rep(names(list.causal.con), sapply(list.causal.con, length)),
                 to = unlist(list.causal.con))

df.causal.con


row.names(df.causal.con) <- c() 

#Select all arcs (not neccessary all identified)
df.all.arcs <- as.data.frame(Sim1$arcs)

#select true model causal pathways
df.true.model <- read.delim("simSolution_1.txt", sep = " ")
colnames(df.true.model)[1] <- "from"
colnames(df.true.model)[2] <- "to"


#####################################################
#####################################################
#####################################################
## order the data

df.true.model <- df.true.model[order(df.true.model$from),]
df.causal.con <- df.causal.con[order(df.causal.con$from),]
df.all.arcs <- df.all.arcs[order(df.all.arcs$from),]


######################################################
######################################################
######################################################
## Score the different models

score.truemodel <-nrow(model.testing$true.causal.model)
score.causal.sim <-nrow(model.testing$sim.causal.con)
score.all.con.sim <-nrow(model.testing$sim.all.con)

score.truemodel <- as.data.frame(score.truemodel)
score.causal.sim <- as.data.frame(score.causal.sim)
score.all.con.sim <- as.data.frame(score.all.con.sim)

######################################################
######################################################
######################################################
# compare the simulated models with the true model

compare_tm_to_sim.causal <- merge(df.true.model,df.causal.con)
compare_tm_to_sim_allpaths <- merge(df.true.model,df.all.arcs)

######################################################
######################################################
######################################################
#score the comparison

score_comp_tm_sc <- nrow(compare_tm_to_sim.causal)
score_comp_tm_sa <- nrow(compare_tm_to_sim_allpaths)

score_comp_tm_sc <- as.data.frame(score_comp_tm_sc)
score_comp_tm_sa <- as.data.frame(score_comp_tm_sa)


######################################################
######################################################
######################################################
#cbind the different scores

Total_score <- cbind(score.truemodel,score.causal.sim,score.all.con.sim, score_comp_tm_sc,score_comp_tm_sa)

Total_score



score_comp_tm_sc
score_comp_tm_sa


#list all datatframes for comparisson
model.testing <- list(sim.causal.con = df.causal.con, sim.all.con= df.all.arcs, true.causal.model=df.true.model)
#order data (makes comparrison easier)
model.testing <- lapply(model.testing, function(x) x[order(x$from),])
model.testing



df.true.model
df.causal.con
a <-compare(df.true.model,df.causal.con, allowAll=TRUE)
a


new_data_frame <- merge(df.true.model,df.all.arcs)
new_data_frame


x <- c("a", "b", "c", "d", "e")
names(x) <- c("foo", "bar", "baz", "qux", "grault")

y <- c("c", "a", "d", "b")
names(y) <- c("bar", "foo", "qux", "corge")

x

y
names(x)
x[intersect(names(x), names(y))] == y[intersect(names(x), names(y))]

model.testing$true.causal.model==model.testing$sim.all.con

compare.models <- compare(model.testing$true.causal.model,model.testing$sim.causal.con)
compare.models

score.truemodel <-nrow(model.testing$true.causal.model)
score.causal.sim <-nrow(model.testing$sim.causal.con)
score.all.con.sim <-nrow(model.testing$sim.all.con)

score.truemodel <- as.data.frame(score.truemodel)
score.causal.sim <- as.data.frame(score.causal.sim)
score.all.con.sim <- as.data.frame(score.all.con.sim)





a<-as.data.frame(c(1,2,3,4,5))
colnames(a)[1] <- "a"
b<-as.data.frame(c(2,3,4,5,9))
colnames(b)[1] <- "a"
c<-as.data.frame(c(0,1,2,3,11))
colnames(c)[1] <- "c"

merge(a,b)


d <-list(a,b,c)
d

e <- lapply(d, function(k) compare(k[[1]],k[[2]]))

e[[1]]$tC

e <-unlist(compare.models$results)
e
compare.tm_scm 



t <- cbind(score.truemodel, score.causal.sim, score.all.con.sim)
c <- cbind(score.truemodel, score.causal.sim, score.all.con.sim)

x <- rbind(t,c)

totscore_tm <- sum(x$score.truemodel)

totscore_tm
head(t)

class(score.truemodel)
dim(model.testing$sim.causal.con)[,1]
dim(model.testing$sim.all.con)







a <- as.data.frame(c(1,2,3,4,5,6))
b <- as.data.frame(c(3,4,5,6,1,2))
a
colnames(a) <- "a"
colnames(b) <- "b"

c <- compare(a,b, allowAll = T )
c$tC
compare.models$true.causal.model

head(model.testing)










dl <- data.frame(ID = rep(names(causal.con), sapply(causal.con, length)),
                 Obs = unlist(causal.con))

dl

myFun <- function(data) {
  temp1 <- sapply(data, is.list)
  temp2 <- do.call(
    cbind, lapply(data[temp1], function(x) 
      data.frame(do.call(rbind, x), check.names=FALSE)))
  cbind(data[!temp1], temp2)
}

mydf <- data.frame(Date = as.Date(c("1978-01-01", "1978-01-02")),
                   V1 = c(10, 10),
                   V2 = c(11, 11))

mydf$V3 <- list(c(1:10),
                c(11:20))
mydf$V4 <- list(c(21:25),
                c(26:30))

mydf
myFun(mydf)








a = t(unlist(lapply(Sim1$nodes, function(k) k$parents)))
a
Sim1$arcs

paste0(round(100*with(df1, mean(A==B))), "%")

Sim1$arcs


causal.con
causal_sim1 <- as.data.frame(causal.con)
causal_sim1
causal_sim1$from <- rownames(causal_sim1)
head(causal_sim1)
causal_sim1 <- causal_sim1[,c(2,1)]
colnames(causal_sim1)[2] <- "to"
rownames(causal_sim1) <- c()

cols <- c("V1", "V2")
wl_sol1 <- read.delim("simSolution_1.txt", sep = " ")
wl_sol1 
wl_sol1$x <- do.call(paste,(c(wl_sol1[cols], sep="")))
wl_sol1 <- as.data.frame(wl_sol1[ , !( names( wl_sol1 ) %in% cols ) ])
colnames(wl_sol1)[1] <- "true"

cols1 <- c("from", "to")
causal_sim1
causal_sim1$x <- do.call(paste,(c(causal_sim1[cols1], sep="")))
causal_sim1
causal_sim1 <- as.data.frame(causal_sim1[ , !( names( causal_sim1 ) %in% cols1 ) ])
causal_sim1
colnames(causal_sim1)[1] <- "sim"


T1 <- list(wl_sol1, causal_sim1)
T1

t2 <- lapply(T1,function(k) paste0(round(100*length(intersect(k[[1]], k[[2]]/nrow(k[[1]]))))))

paste0(round(100*length(intersect(df1$A, df1$B))/nrow(df1)), "%")

causal.con = unlist(lapply(Sim1$nodes, function(k) k$children))

wl_sol1==causal_sim1

wl_sol1 <- wl_sol1[order(wl_sol1$V1),]
wl_sol1
c <- list(a,b)
c









out <-vector()
for (i in 1:length(Sim1$nodes)) {
  out[names(Sim1$nodes)[i]] = Sim1$nodes[[i]]$mb[Sim1$nodes[[i]]$mb%in%Sim1$nodes[[i]]$children]
}
out
names(Sim1$nodes)




comparison <- compare(wl_sol1, causal_sim1, allowAll = T)
comparison$tMpartial

(which(causal_sim1%in%wl_sol1))

a1 <- data.frame(a = 1:5, b = letters[1:5])
a2 <- data.frame(a = 1:3, b = letters[1:3])
comparison <- compare(a1,a2,allowAll=TRUE)
comparison$tM

a <- c(1,2,3,4,5)

a
b <- c(1,2,3)

b

b==a
b
c <- b[,c(2,1)]
c


head(test)
test <- causal.con$from <-rownames(causal.con)

head(test)



for(i in names(causal.con)) {
  which(Sim1$arcs[,1]%in%i)
}
Sim1$nodes$prs1

Sim1$nodes$prs2$mb[Sim1$nodes$prs2$mb%in%Sim1$nodes$prs2$children]
Sim1$nodes$prs2$mb[!Sim1$nodes$prs2$mb%in%Sim1$nodes$prs2$children]


Sim1$arcs

wl_sol1


##Solution
wl_sol1 <- read.delim("simSolution_1.txt", sep = " ")
wl_sol1
bl1 <- matrix(c("t1","t1","t1", "prs1" ,"prs2","prs3"),ncol=2)
bl2 <- matrix(c("t2","t2","t2", "prs1" ,"prs2","prs3"),ncol=2)
bl3 <- matrix(c("t3","t3","t3", "prs1" ,"prs2","prs3"),ncol=2)
bl4 <- matrix(c("prs2", "prs3"),ncol=2)
bl5 <- matrix(c("prs3", "prs2"),ncol=2)
bl <- rbind(bl1,bl2,bl3, bl4,bl5)
Sol1 <- iamb(Simdata1, test="cor", blacklist = bl, whitelist = wl_sol1)
graphviz.plot(Sol1, main = "SOL1")




## SIMULATION2
par(mfrow=c(2,1))
Simdata2 <- read.delim("simData_2.txt", header=T, sep = " ")
head(Simdata2)

bl1 <- matrix(c("t1","t1","t1", "prs1" ,"prs2","prs3"),ncol=2)
bl2 <- matrix(c("t2","t2","t2", "prs1" ,"prs2","prs3"),ncol=2)
bl3 <- matrix(c("t3","t3","t3", "prs1" ,"prs2","prs3"),ncol=2)
bl <-rbind(bl1, bl2,bl3)
Sim2 <- iamb(Simdata2, test="cor", blacklist = bl)
graphviz.plot(Sim2, main = "SIM2")


wl_sol2 <- read.delim("simSolution_2.txt", sep = " ")
wl_sol2

bl1 <- matrix(c("t1","t1","t1", "prs1" ,"prs2","prs3"),ncol=2)
bl2 <- matrix(c("t2","t2","t2", "prs1" ,"prs2","prs3"),ncol=2)
bl3 <- matrix(c("t3","t3","t3", "prs1" ,"prs2","prs3"),ncol=2)
bl4 <- matrix(c("prs1", "prs3"),ncol=2)
bl5 <- matrix(c("t2", "t3"),ncol=2)
bl6 <- matrix(c("prs3", "prs1"),ncol=2)
bl7 <- matrix(c("t3", "t2"),ncol=2)
bl8 <- matrix(c("prs3", "prs2"),ncol=2)
bl9 <- matrix(c("prs2", "prs3"),ncol=2)

bl <- rbind(bl1,bl2,bl3,bl4,bl5,bl6,bl7,bl8,bl9)

Sol2 <- iamb(Simdata2, test="cor", blacklist = bl, whitelist = wl_sol2)
graphviz.plot(Sol2, main = "SOL2")




##

par(mfrow=c(2,1))

Simdata3 <- read.delim("simData_3.txt", header=T, sep = " ")
head(Simdata3)

bl1 <- matrix(c("t1","t1","t1", "prs1" ,"prs2","prs3"),ncol=2)
bl2 <- matrix(c("t2","t2","t2", "prs1" ,"prs2","prs3"),ncol=2)
bl3 <- matrix(c("t3","t3","t3", "prs1" ,"prs2","prs3"),ncol=2)
bl <- rbind(bl1,bl2,bl3)
Sim3 <- iamb(Simdata3, test="cor", blacklist = bl)
graphviz.plot(Sim3, main = "SIM3")


wl_sol3 <- read.delim("simSolution_3.txt", sep = " ")
wl_sol3

bl1 <- matrix(c("t1","t1","t1", "prs1" ,"prs2","prs3"),ncol=2)
bl2 <- matrix(c("t2","t2","t2", "prs1" ,"prs2","prs3"),ncol=2)
bl3 <- matrix(c("t3","t3","t3", "prs1" ,"prs2","prs3"),ncol=2)
bl4 <- matrix(c("prs1", "prs3"),ncol=2)
bl5 <- matrix(c("prs1", "prs3"),ncol=2)
bl6 <- matrix(c("prs3", "prs1"),ncol=2)
bl7 <- matrix(c("prs2", "prs3"),ncol=2)
bl8 <- matrix(c("prs3", "prs2"),ncol=2)
bl9 <- matrix(c("prs2", "prs3"),ncol=2)

bl <- rbind(bl1,bl2,bl3,bl4,bl5,bl6,bl7,bl8,bl9)

Sol3 <- iamb(Simdata3, test="cor", blacklist = bl, whitelist = wl_sol3)
graphviz.plot(Sol3, main = "SOL2")


##

Simdata4 <- read.delim("simData_4.txt", header=T, sep = " ")
head(Simdata4)

bl1 <- matrix(c("t1","t1","t1", "prs1" ,"prs2","prs3"),ncol=2)
bl2 <- matrix(c("t2","t2","t2", "prs1" ,"prs2","prs3"),ncol=2)
bl3 <- matrix(c("t3","t3","t3", "prs1" ,"prs2","prs3"),ncol=2)

Sim4 <- iamb(Simdata4, test="cor", blacklist = bl)
graphviz.plot(Sim4)


Simdata5 <- read.delim("simData_5.txt", header=T, sep = " ")
head(Simdata5)

bl1 <- matrix(c("t1","t1","t1", "prs1" ,"prs2","prs3"),ncol=2)
bl2 <- matrix(c("t2","t2","t2", "prs1" ,"prs2","prs3"),ncol=2)
bl3 <- matrix(c("t3","t3","t3", "prs1" ,"prs2","prs3"),ncol=2)

Sim5 <- iamb(Simdata5, test="cor", blacklist = bl)
graphviz.plot(Sim5)

##

Simdata6 <- read.delim("simData_6.txt", header=T, sep = " ")
head(Simdata6)

bl1 <- matrix(c("t1","t1","t1", "prs1" ,"prs2","prs3"),ncol=2)
bl2 <- matrix(c("t2","t2","t2", "prs1" ,"prs2","prs3"),ncol=2)
bl3 <- matrix(c("t3","t3","t3", "prs1" ,"prs2","prs3"),ncol=2)

Sim6 <- iamb(Simdata6, test="cor", blacklist = bl)
graphviz.plot(Sim6)

par(mfrow=c(2,2)) 

Simdata7 <- read.delim("simData_7.txt", header=T, sep = " ")
head(Simdata7)

bl1 <- matrix(c("t1","t1","t1", "prs1" ,"prs2","prs3"),ncol=2)
bl2 <- matrix(c("t2","t2","t2", "prs1" ,"prs2","prs3"),ncol=2)
bl3 <- matrix(c("t3","t3","t3", "prs1" ,"prs2","prs3"),ncol=2)

Sim7 <- iamb(Simdata7, test="cor", blacklist = bl)
graphviz.plot(Sim7)


Simdata8 <- read.delim("simData_8.txt", header=T, sep = " ")
head(Simdata8)

bl1 <- matrix(c("t1","t1","t1", "prs1" ,"prs2","prs3"),ncol=2)
bl2 <- matrix(c("t2","t2","t2", "prs1" ,"prs2","prs3"),ncol=2)
bl3 <- matrix(c("t3","t3","t3", "prs1" ,"prs2","prs3"),ncol=2)

Sim8 <- iamb(Simdata8, test="cor", blacklist = bl)
graphviz.plot(Sim8)


Simdata9 <- read.delim("simData_9.txt", header=T, sep = " ")
head(Simdata9)

bl1 <- matrix(c("t1","t1","t1", "prs1" ,"prs2","prs3"),ncol=2)
bl2 <- matrix(c("t2","t2","t2", "prs1" ,"prs2","prs3"),ncol=2)
bl3 <- matrix(c("t3","t3","t3", "prs1" ,"prs2","prs3"),ncol=2)

Sim9 <- iamb(Simdata9, test="cor", blacklist = bl)
graphviz.plot(Sim9)


Simdata10 <- read.delim("simData_10.txt", header=T, sep = " ")
head(Simdata10)

bl1 <- matrix(c("t1","t1","t1", "prs1" ,"prs2","prs3"),ncol=2)
bl2 <- matrix(c("t2","t2","t2", "prs1" ,"prs2","prs3"),ncol=2)
bl3 <- matrix(c("t3","t3","t3", "prs1" ,"prs2","prs3"),ncol=2)

Sim10 <- iamb(Simdata10, test="cor", blacklist = bl)
graphviz.plot(Sim10)



