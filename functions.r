make.lists <- function (mat) {
  traits <- mat[,1]
  prs <- mat[,2]
  d <- t(combn(c(traits, prs), 2))
  d <- rbind(d, d[, c(2,1)])

  out <- vector()
  i <- 1
  for (i in 1:nrow(d)) {
    temp <- d[i,] 
    if ((temp[1] %in% traits) && (temp[2] %in% traits)) {
      out <- append(out, i)
      next
    }
    numb.trait <- which(mat==temp[1], arr.ind=T)[,'row']
    numb.prs <- which(mat==temp[2], arr.ind=T)[,'row']
    if ((temp[1] %in% prs) && (temp[2] %in% traits) && (numb.trait==numb.prs)) {
      out <- append(out, i)
    }
  }
  white <- d[out,]
  black <- d[!c(1:nrow(d))%in%out,]
  return(list(white, black))
}

prs <- function (G, sum.stat, p) {
  selected.snp <- which(sum.stat$p <= p) 
  message(paste0("Selected ", length(selected.snp), " snps"))
  out <- G[,selected.snp]%*%sum.stat$beta[selected.snp]
  return(out)
}

gwas <- function (y, X, bim, covar=NULL, label=NULL) {
  est.effect <- matrix(0, ncol=4, nrow=ncol(X))
  if (is.null(covar)) {
    for (i in 1:ncol(X)) {
      est.effect[i,] <- summary(lm(y~X[,i]))$coef[2, ]
    }
  } else {
    for (i in 1:ncol(X)) {
      est.effect[i,] <- summary(lm(y~X[,i] + covar))$coef[2, ]
    }
  }
  est.effect <- as.data.frame(est.effect)
  names(est.effect) <- c("beta", "se", "tvalue", "pval")
  est.effect$SNP  <- bim$V2
  est.effect$effect_allele  <- bim$V5
  est.effect$other_allele  <- bim$V6

  if(!is.null(label)) est.effect$Phenotype <- label
  return(est.effect)
}

define.causal <- function (p, effectDist) {
  np = 1:p
  uniq <- diag(effectDist)
  diag(effectDist) <- NA

  nn <- colnames(effectDist)
  l <- list() 
  u <- 1
  for (i in uniq) {
    temp <- np[!np%in%unlist(l)]
    l[[nn[[u]]]] <- sample(temp, floor(length(np)*i))
    u <- u + 1
  }
  u <- 1
  for (i in colnames(effectDist)) {
    temp <- np[!np%in%unlist(l)]
    temp.dist <- effectDist[,i]
    temp.dist <- temp.dist[!is.na(temp.dist)]
    for (ii in names(temp.dist)) {
      temp.s <- sample(temp, floor(length(np)*temp.dist[ii]))
      l[[ii]] <- append(l[[ii]], temp.s)
      l[[i]] <- append(l[[i]], temp.s )
      temp <- np[!np%in%unlist(l)]
    }
  }
  return(l)
}

simulate.causal.prs <- function (n, p, effectmat, dd, phenoeffect, wished.effect=0.6) {
  ca <- define.causal(p, effectmat)
  ef1 <-rep(0, p)
  ef2 <-rep(0, p) 
  ef3 <-rep(0, p) 
  efca <- list(ef1, ef2, ef3)
  for (item in 1:length(efca)) {
    efca[[item]][ca[[item]]] <- 1 
    efca[[item]][ca[[item]]] <- efca[[item]][ca[[item]]] /
    (sum(efca[[item]][ca[[item]]])/wished.effect)
  }

  t1  <-  1*(dd%*%efca[[1]]) + rnorm(n, 0, sqrt(1-wished.effect))
  t2  <-  1*(dd%*%efca[[2]]) + phenoeffect[2,1]*t1 + rnorm(n, 0, sqrt(1-wished.effect))
  t3  <-  1*(dd%*%efca[[3]]) + phenoeffect[3,1]*t1 + phenoeffect[3,1]*t2 +
  rnorm(n, 0, sqrt(1-wished.effect))


  gwas1 <- gwas(t1, dd)
  gwas2 <- gwas(t2, dd)
  gwas3 <- gwas(t3, dd)

  prs1 <- prs(dd, gwas1, 0.01)
  prs2 <- prs(dd, gwas2, 0.01)
  prs3 <- prs(dd, gwas3, 0.01)

  output <- data.frame(t1, t2, t3,
                       prs1, prs2, prs3)
  return(output)
}

qq.gwas <- function (p) {

  nn <- -log10(1:length(p)/length(p))
  p <- -log10(sort(p))
  plot(nn,p, ylab='-log10(p)', xlab='-log10(null)')
  abline(0,1, col='red')
}

make.arcs <- function (pheno.effectmat) {
  base <- matrix(c(paste0('snp',1:3), paste0('t', 1:3)), ncol=2)
  causal.path <- which(pheno.effectmat!=0, T)
  causal.path <- matrix(c(colnames(pheno.effectmat)[causal.path[,2]],
                        rownames(pheno.effectmat)[causal.path[,1]]), ncol=2)
  return(rbind(base, causal.path))
}

generate.random.causes <- function(effect.size) {
  pheno.effectmat <- matrix(0, nrow=3, ncol=3)
  pheno.effectmat[lower.tri(pheno.effectmat)] <- sample(effect.size, 3)
  diag(pheno.effectmat) <- NA
  colnames(pheno.effectmat) <- paste0('t', 1:3)
  rownames(pheno.effectmat) <- paste0('t', 1:3)
  return(pheno.effectmat)
}

simulate.causal.snp <- function (n, p, effectmat, dd, phenoeffect, wished.effect=0.6, bim, dogwas=F) {
  effectmat <- genetic.effectmat
  phenoeffect <- pheno.effectmat
  ca <- define.causal(p, effectmat)
  ef1 <-rep(0, p)
  ef2 <-rep(0, p) 
  ef3 <-rep(0, p) 
  efca <- list(ef1, ef2, ef3)
  for (item in 1:length(efca)) {
    efca[[item]][ca[[item]]] <- 1 
    efca[[item]][ca[[item]]] <- efca[[item]][ca[[item]]] /
    (sum(efca[[item]][ca[[item]]])/wished.effect)
  }

  t1  <-  1*(dd%*%efca[[1]]) + rnorm(n, 0, sqrt(1-wished.effect))
  t2  <-  1*(dd%*%efca[[2]]) + phenoeffect[2,1]*t1 + rnorm(n, 0, sqrt(1-wished.effect))
  t3  <-  1*(dd%*%efca[[3]]) + phenoeffect[3,1]*t1 + phenoeffect[3,1]*t2 +
  rnorm(n, 0, sqrt(1-wished.effect))


  if(dogwas==T) {
  gwas1 <- gwas(t1, dd, bim, label="T1")
  gwas2 <- gwas(t2, dd, bim, label="T2")
  gwas3 <- gwas(t3, dd, bim, label="T3")
  gwas1$causal <- 0
  gwas1$causal[ca[[1]]] <- 1
  gwas2$causal <- 0
  gwas2$causal[ca[[2]]] <- 1
  gwas3$causal <- 0
  gwas3$causal[ca[[3]]] <- 1
  gwas_results <- rbind(gwas1, gwas2, gwas3)
  } else {
    gwas_results <- NULL
  }

  snp1 <- as.data.frame(dd[,ca[[1]]])
  names(snp1) <- paste0('snp1.', 1:ncol(snp1))
  snp2 <- as.data.frame(dd[,ca[[2]]])
  names(snp2) <- paste0('snp2.', 1:ncol(snp2))
  snp3 <- as.data.frame(dd[,ca[[3]]])
  names(snp3) <- paste0('snp3.', 1:ncol(snp3))

  output <- data.frame(t1, t2, t3,
                       snp1, snp2, snp3)
  outout <- list(output, gwas_results)
  return(outout)
}


test.mr <- function (sim) {
  possibilities <- list(c("T1", "T2"), c("T1", "T3"), c("T2", "T3"))
  out <- list() 

  u <- 1 
  for (item in possibilities) {
    expo <- sim %>%
    filter(Phenotype==item[1], causal==1) %>%
    format_data(., type='exposure')

    outcome <-  sim %>%
    filter(Phenotype==item[2]) %>%
    format_data(., type='outcome', phenotype_col="Phenotype")

    dat <- harmonise_data(exposure_dat = expo, 
                          outcome_dat = outcome, action=1)
    out[[u]] <- mr(dat, method_list="mr_ivw")
    u <- u + 1
  }

  out <- do.call('rbind', out)
  out <- out %>% 
  select(exposure, outcome, b, se, pval) %>%
  mutate(causal=pval<=(0.05/length(possibilities)))

  return(out)
}

compare.outcomes <- function (real, estimated, causal='causal') {
  real <- arcs[!is.na(arcs$value),c("from", "to")]
  estimated <- outcome[outcome[[causal]], c("from", "to")]
  real <- mutate_each(real, funs(tolower))
  estimated <- mutate_each(estimated, funs(tolower))


  if(nrow(real)==0) names(real) <- "empty"
  if(nrow(estimated)==0) names(estimated) <- "empty"

  if (nrow(real)==0 & nrow(estimated) > 0) {
    overestimated <- estimated 
    not.estimated <- NULL
  }

  if (nrow(real)>0 & nrow(estimated) == 0) {
    not.estimated <- real 
    overestimated  <- NULL
  }

  if (nrow(real)>0 & nrow(estimated) > 0) {
    not.estimated <- anti_join(real, estimated)
    overestimated <- anti_join(estimated, real)
  }
  return(list("not.estimated"=not.estimated, "overestimated"=overestimated)) 
}

check.mr.performance <- function (k) {
  out <- vector('numeric', length=2)
  for (item in 1:length(k)) {
    if(is.null(k[item][[1]])) { 
      out[item] <- 0
    } else {
      out[item] <- nrow(k[item][[1]])
    }
  }
  return(out)
}
