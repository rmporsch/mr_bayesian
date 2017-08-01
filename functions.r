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

gwas <- function (y, X, covar=NULL) {
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
  names(est.effect) <- c("beta", "SE", "tvalue", "p")
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

simulate.complex <- function (n, p, effectmat, dd, phenoeffect, wished.effect=0.6) {
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
  base <- matrix(c(paste0('prs',1:3), paste0('t', 1:3)), ncol=2)
  causal.path <- which(pheno.effectmat!=0, T)
  causal.path <- matrix(c(colnames(pheno.effectmat)[causal.path[,2]],
                        rownames(pheno.effectmat)[causal.path[,1]]), ncol=2)
  return(rbind(base, causal.path))
}

generate.random.causes <- function(effect.size) {
  effect.size <- c(0, 0.3, 0.7, 0.3,0)
  pheno.effectmat <- matrix(0, nrow=3, ncol=3)
  pheno.effectmat[lower.tri(pheno.effectmat)] <- sample(effect.size, 3)
  diag(pheno.effectmat) <- NA
  colnames(pheno.effectmat) <- paste0('t', 1:3)
  rownames(pheno.effectmat) <- paste0('t', 1:3)
  return(pheno.effectmat)
}
