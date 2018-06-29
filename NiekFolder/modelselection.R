################################################
# test model selection
################################################
library(gtools)
source('./SPARSE_PCA_wRandomStart.R')
require(Matrix)
rm(list=ls())

#adapt the data generating to be more
#so you can give the sparcity in percentage, and you can set the amount of error
generateCommonSpecific <- function(x, nx, nfactors, p, 
                                   coefFixed = TRUE, sparsity){
    X <- matrix(rnorm(x*nx, 0, 3), nx, x)
    X <- scaleData(X)
    XTX <- t(X) %*% X
    V <- svd(X)$v
    W <- V[, 1:nfactors]
    W[1:floor(nrow(W)/2), 2] <- 0
    W[(floor(nrow(W)/2)+1):nrow(W), 3] <- 0

    if(coefFixed){
       fixW <- matrix(rnorm(dim(X)[2]*nfactors, 100), dim(X)[2], nfactors) 
       fixW[(floor(nrow(W)/2)+1):nrow(W), 3] <- 0
       fixW[1:(floor(nrow(W)/2)), 2] <- 0
    } else {
        fixW <- matrix(rnorm(dim(X)[2]*nfactors, 100), dim(X)[2], nfactors) 
    }

    #create sparsity in the columns of W based on a percentage
    for(i in 1:ncol(W)){
        spar <- quantile(abs(W[, i]), probs = sparsity[i])
        W[abs(W[, i]) < spar, i] <- 0
    }

    XandP <- betterXandP(X, W)
    Xtrue <- XandP$X
    P <- XandP$P

    E <- matrix(rnorm(x*nx, 0, 1), nx, x)
    g = sqrt((var(as.vector(Xtrue))*p) / (var(as.vector(E)) * (1 - p))) 
    X <- Xtrue + E*g

    SStrue <- var(as.vector(Xtrue))
    SSX <- var(as.vector(X))

    return(list(X = X, W = W, P = P, errorRatio = 1 - SStrue/SSX,
    percentageZeroes = apply(W, 2, countZero), fixW = fixW))
}


betterXandP <- function(X, W){
    for(i in 1:1){
        X <- scaleData(X)
        XTX <- t(X) %*% X
        SVD <- svd(XTX %*% W)   
        P <- SVD$u %*% t(SVD$v)  
        X <- X %*% W %*% t(P)
    }
    return(list(X = X, P = P))
}

generateCommonSpecific(60, 60, 3, 0.05, coefFixed = TRUE, sparsity = c(0.05, 0.5, 0.5))

W

vars <- c(5,5)
components <- 3

#get all common and distinctive structures
allCommonSpecific <- function(vars, components){

    blocks <- length(vars)

    W  <- matrix(NA, sum(vars), components)
    cd <- as.matrix(expand.grid(rep(list(0:1), blocks))[-1, ])
    commonSpecific <- combinations(n = nrow(cd), r = components,
                                  v = c(1:nrow(cd)), repeats.allowed = TRUE)
    allpossibleWmatrices <- rep(list(NA), nrow(commonSpecific))

    for(i in 1:nrow(commonSpecific)){
        for(j in 1:ncol(commonSpecific)){
            W[, j] <-  rep(cd[commonSpecific[i,j], ], times = vars)
            allpossibleWmatrices[[i]] <- W
    }
    }
    return(allpossibleWmatrices)
}

allpossibleWmatrices

allCommonSpecific(c(10, 10), 3)


# Eigenvector crossvalidation method
EigenVectorCV2 <- function(inputDat, ridge, lasso, nrFolds, nfactors, fixW, nScale = 0){

    folds <- rep_len(1:nrFolds, nrow(inputDat))
    #folds <- sample(folds)
    cvError  <- matrix(NA, nrow(inputDat), ncol(inputDat))
    MSEkthFold <- rep(NA, nrFolds)

    for(a in 1:nrFolds){

        # actual split of the data
        fold <- which(folds == a)
        trainDat <- inputDat[-fold, ]
        testDat <- inputDat[fold, , drop = FALSE]

		#Scale training and the testing data
		means <- apply( trainDat, 2, mean )		
		trainDat <- t( t(trainDat) - means)
		sdTrain <- apply(trainDat, 2, function(x) sqrt(sum( x^2, na.rm = TRUE ) / (length(na.omit(x)) - nScale )))
		trainDat <- t( t(trainDat)/sdTrain )
		
		testDat <-  t( t(testDat) - means)
		testDat <- t( t(testDat)/sdTrain ) 
		
		rm(means, sdTrain)		
        
        #model estimation
        res <- sparsePCAMultistart(trainDat, nfactors = nfactors, RIDGE = ridge, LASSO = lasso,
                                    percentage, fixW = fixW,
                                    maxItrOuterloop = 100000,
                                    percentageMode = FALSE)

        #Eigenvector crossvalidation Bro Kiers
        pred <- matrix(NA, nrow(testDat), ncol(testDat))
        print(a)
        for(i in 1:ncol(inputDat)){
            TMinusVariableJ <- testDat[, -i] %*% res$W[-i, ]
            pred[, i] <- TMinusVariableJ %*% res$P[i, ] 
            }
        cvError[fold, ] <- (testDat - pred)^2
        MSEkthFold[a]  <-  mean(cvError[fold, ]) 
        }

    return(list(cvError = cvError, MSE = mean(MSEkthFold), MSEkthFold = MSEkthFold, 
                stdError = sd(MSEkthFold) / sqrt(nrFolds)))
}




source('./SPARSE_PCA_wRandomStart.R')


dat$X
dat$W

dat$errorRatio

fit <- sparsePCAMultistart(dat$X, nfactors = 3, RIDGE = 0.1, LASSO = 0.1,
                                    percentage = dat$percentageZeroes, fixW = dat$fixW,
                                    maxItrOuterloop = 100000,
                                    percentageMode = FALSE)




dat <- generateCommonSpecific(1000, 100, 3, p = 0.05,
                              coefFixed = TRUE, sparsity = c(0.6, 0.6, 0.6))

dat$X <- scale(dat$X, scale = FALSE)
aap <- sparseSCAcpp(dat$X, Q = 3, RIDGE = 0.1, LASSO = rep(0.1 , 3), fixW = dat$fixW, 
             maxItrOuterloop = 1000, nStarts = 1, print = TRUE, tol = 10^-10)



fit$loss
mod

require(Rcpp)
require(RcppArmadillo)

sourceCpp('./cppTesting/test.cpp')


dat$X <- scale(dat$X, scale = FALSE)
aap <- sparseSCAcpp(dat$X, Q = 3, RIDGE = 0.1, LASSO = rep(0.1 , 3), fixW = dat$fixW, 
             maxItrOuterloop = 1000, nStarts = 1, print = TRUE, tol = 10^-10)


fit <- sparsePCAMultistart(dat$X, nfactors = 3, RIDGE = 0.1, LASSO = 0.1,
                                    percentage = dat$percentageZeroes, fixW = dat$fixW,
                                    maxItrOuterloop = 100000,
                                    percentageMode = FALSE)

system.time( 
aap <- sparseSCAcpp(dat$X, Q = 3, RIDGE = 0.1, LASSO = rep(0.1 , 3), fixW = dat$fixW, 
             maxItrOuterloop = 1000, nStarts = 1, print = FALSE, tol = 10^-10)
)

system.time(
fit <- sparsePCAMultistart(dat$X, nfactors = 3, RIDGE = 0.1, LASSO = 0.1,
                                    percentage = dat$percentageZeroes, fixW = dat$fixW,
                                    maxItrOuterloop = 100000,
                                    percentageMode = FALSE)
)




head(aap$W, 40)
head(fit$W, 40)

dat$W



dat <- generateCommonSpecific(500, 100, 3, p = 0.05,
                              coefFixed = TRUE, sparsity = c(0.6, 0.6, 0.6))
dat$X

dat$fixW


lasso <- seq(0, 0.1, by = 0.001)
folds <- 10
res <- rep(NA, length(lasso)) 
stdErrorRes <- rep(NA, length(lasso)) 

for(i in 1:length(lasso)){
    error <- EigenVectorCV2(dat$X, ridge = 0.1, lasso = lasso[i] ,
                            nrFolds = folds, nfactors = 3, fixW = dat$fixW)
    res[i] <- error$MSE
    stdErrorRes[i] <- error$stdError
}
plot(res)

x <- 1:length(lasso)
plot(x, res,
    ylim=range(c(res - stdErrorRes, res + stdErrorRes)),
    pch=19, xlab="Measurements", ylab="Mean +/- SD",
    main="Scatter plot with std.dev error bars"
)
# hack: we draw arrows but with very special "arrowheads"
arrows(x, res - stdErrorRes, x, res + stdErrorRes, length=0.05, angle=90, code=3)


#get the models within the One StdError of the model with the smallest MSE
withinOneStdError <- res < res[which.min(res)] + stdErrorRes
lasso[withinOneStdError]
best <- max(lasso[withinOneStdError])
best

dat <- generateCommonSpecific(500, 100, 3, p = 0.05, coefFixed = TRUE, sparsity = c(0.6, 0.6, 0.6))
fit <- sparsePCAMultistart(dat$X, nfactors = 3, RIDGE = 0.1, LASSO = 0.006,
                                    percentage = dat$percentageZeroes, fixW = dat$fixW,
                                    maxItrOuterloop = 100000,
                                    percentageMode = FALSE)


correctlyClassified(dat$W, fit$W)
tuckerCongruence(dat$W, fit$W)

head(dat$W)
head(fit$W)


CVerror <- EigenVectorCV2(dat$X, ridge = 0, lasso = 0.05, nrFolds = 10, nfactors = 3, fixW = dat$fixW)

for(

res <- seq(0, 1, by = 0.01)



allCombn <- allCommonSpecific(c(5, 5), 3)
length(allCombn)

dat <- generateCommonSpecific(x = 40 , nx = 80, nfactors = 3, p = 0.05, coefFixed = TRUE, 
                       sparsity = c(.6, .6,.6))

dat

allCombn

folds <- 10
res <- rep(NA, folds) 
stdErrorRes <- rep(NA, folds)

for(i in 1:length(allCombn)){
    error <- EigenVectorCV2(dat$X, ridge = 0, lasso = 0 , nrFolds = folds, nfactors = 3, fixW = allCombn[[i]])
    res[i] <- error$MSE
    stdErrorRes[i] <- error$stdError
}

x <- 1:length(allCombn)
plot(x, res,
    ylim=range(c(res - stdErrorRes, res + stdErrorRes)),
    pch=19, xlab="Measurements", ylab="Mean +/- SD",
    main="Scatter plot with std.dev error bars"
)
# hack: we draw arrows but with very special "arrowheads"
arrows(x, res - stdErrorRes, x, res + stdErrorRes, length=0.05, angle=90, code=3)

#get the models within the One StdError of the model with the smallest MSE
withinOneStdError <- res < res[which.min(res)] + stdErrorRes
allCombn[withinOneStdError]

withinOneStdError

minimum <- which.min(res)
order(res)

allCombn
?order

minimum

allCombn[minimum]
minimum


dat$fixW
dim(dat$X)
test <- EigenVectorCV2(dat$X, ridge = 0, lasso = 0 , nrFolds = folds, nfactors = 3, fixW = allCombn[[1]])

test



