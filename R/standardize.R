# Original script by John Dunn, retrieved from https://osf.io/n62y4/files/osfstorage

standardize <- function(standard, target, joint) {
  # returns standardized linear components from fitSR output
  # standard = fitSR output from standard data
  # target = fitSR output from target data
  # joint = fitSR output from joint standard+target data
  output <- NULL
  ns <- length(standard$z)
  ms <- length(standard$x)
  nt <- length(joint$z) - ns
  mt <- length(joint$x) - ms
  x <- cbind(rep(1, ns), joint$z[1:ns])
  b <- solve(t(x) %*% x) %*% t(x) %*% standard$z
  zs <- joint$z[1:ns] * b[2] + b[1]
  zt <- joint$z[(ns + 1):(ns + nt)] * b[2] + b[1]
  a <- standard$design
  xs <- MASS::ginv(t(a) %*% a) %*% t(a) %*% zs
  a <- target$design
  xt <- MASS::ginv(t(a) %*% a) %*% t(a) %*% zt
  output$zs <- zs
  output$zt <- zt
  output$xs <- xs
  output$xt <- xt
  return(output)
}
