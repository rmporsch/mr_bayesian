traits <- paste0("t", 1:3)
prs <- paste0("prs", 1:3)
mat <- matrix(c(traits, prs), ncol=2)

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
make.lists(mat)
