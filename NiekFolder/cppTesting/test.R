install.packages("Rcpp")
install.packages("RcppArmadillo")

rm(list=ls())

sourceCpp('./test.cpp')
require(Rcpp)
require(RcppArmadillo)

sumC(1:1000)

X <- matrix(rnorm(60*60), 60, 60) 
X[1, 1] <- 0

sourceCpp('./davideCode.cpp')


W <- svd(X )$v[, 1:3]


check <- svd(X %*% t(X) %*% W)
p <- check$u %*% t(check$v)

fixW <- matrix(1, 60,3)
fixW[1:10, 1] <- 0 
fixW[60, 1]  <- 0



sourceCpp('./test.cpp')
cpp <- sparseSCAcpp(X = X, Q =3 , RIDGE = 4, LASSO = c(2, 2, 2),
                    fixW = fixW, maxItrOuterloop = 3, nStarts =2, print = TRUE)

head(cpp$P)
head(p)

