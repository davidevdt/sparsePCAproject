####################
# Functions to get:
# - All common and distinctive structures
# - EigenVector crossvalidation
# - analysis of the family data


library(gtools)
source('./SPARSE_PCA_wRandomStart.R')

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
        }
        allpossibleWmatrices[[i]] <- W

    }
    return(allpossibleWmatrices)
}



# Eigenvector crossvalidation method
EigenVectorCV2 <- function(inputDat, ridge, lasso, nrFolds, nfactors, fixW){

    folds <- rep_len(1:nrFolds, nrow(inputDat))
    #folds <- sample(folds)
    cvError  <- matrix(NA, nrow(inputDat), ncol(inputDat))

    for(a in 1:nrFolds){

        # actual split of the data
        fold <- which(folds == a)
        trainDat <- inputDat[-fold, ]
        testDat <- inputDat[fold, , drop = FALSE]
        
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
        }
    return(list(cvError = cvError))
}



#########################
# Family data analysis


#to obtain the data see README
#load('./family_data.RData')

str(family_data)

child <- family_data$child
mom <- family_data$mom
dad <- family_data$dad

namesDadMomChild <- c(names(dad), names(mom), names(child))

ncol(child)
ncol(dad)
ncol(mom)

rownames(child) <- 1:nrow(child)
rownames(mom) <- 1:nrow(mom)
rownames(dad) <- 1:nrow(dad)

dadMomChild <- cbind(dad, mom, child)
dadMomChild <- scale(as.matrix(dadMomChild))

#scree plot
eigenV <- eigen(cov(scale(as.matrix(dadMomChild))))$values 
plot(eigenV / sum(eigenV))


########################################################################

set.seed(1)

#get all common and distinctive structures
allWs3comp <- allCommonSpecific(vars = c(8, 8, 7), components = 3)


#LOOCV on all common distinctive structures
res3comp <- matrix(NA, length(allWs3comp), 2)
for(i in 1:length(allWs3comp)){
    out  <- EigenVectorCV2(inputDat = dadMomChild, ridge = 0, fixW = allWs3comp[[i]], 
                           lasso = 0, nrFolds = 195, nfactors = 3)
    res3comp[i] <- mean(out$cvError) 
}


bestModel3comp <- which.min(res3comp[, 1])
allWs3comp[[bestModel3comp]]

#plot MPRESS all common distinctive structures
plot(1:nrow(res3comp), res3comp[, 1],
    #ylim = range(c(min(res3comp[, 1] - res3comp[, 2]), max(res3comp[, 1] + res3comp[, 2]))),
    xlab = lasso,
    main = "MPSE of all possible models, the best model  with a 1 sd-error bars",
    col = color2,
    pch = rep(1, nrow(res3comp))
)

#tune lasso best model
lasso3 <- seq(0, 1, by = 0.01)
lasso3

lasso3comp <- matrix(NA, length(lasso3), 2)

for(i in 1:length(lasso3)){
    out  <- EigenVectorCV2(inputDat = dadMomChild, ridge = 0, lasso = lasso3[i], nrFolds = 195, nfactors = 5, fixW = allWs3comp[[bestModel3comp]])
    lasso3comp[i, 1] <- mean(out$cvError)
}


plot(lasso3comp[, 1])
bestLasso <- which.min(lasso3comp[, 1])
lasso3[bestLasso]


res <- sparsePCAMultistart(dadMomChild, nfactors = 3, RIDGE = 0, LASSO = lasso3[bestLasso],
                            percentage, fixW = allWs3comp[[bestModel3comp]],
                            maxItrOuterloop = 100000,
                            percentageMode = FALSE)

rownames(res$W) <- namesDadMomChild
res$W


#run the models for in Table 1
res1 <- sparsePCAMultistart(dadMomChild, nfactors = 3, RIDGE = 0, LASSO = 0,
                            percentage, fixW = matrix(1, ncol(dadMomChild), 3),
                            maxItrOuterloop = 100000,
                            percentageMode = FALSE)

res2 <- sparsePCAMultistart(dadMomChild, nfactors = 3, RIDGE = 0, LASSO = 0,
                            percentage, fixW = allWs3comp[[bestModel3comp]],
                            maxItrOuterloop = 100000,
                            percentageMode = FALSE)

res3 <- sparsePCAMultistart(dadMomChild, nfactors = 3, RIDGE = 0, LASSO = lasso3[bestLasso],
                            percentage, fixW = allWs3comp[[bestModel3comp]],
                            maxItrOuterloop = 100000,
                            percentageMode = FALSE)

#function to get the variance accounted for
VAFw <- function(inputDat, analysisResults){
    #get sum of squares of the data
    SS <- sum(scale(inputDat, scale = FALSE)^2)
    SSmodel <- sum((inputDat - inputDat %*% analysisResults$W %*% t(analysisResults$P))^2)
    return(1 - SSmodel/SS)
}


#VAFS for the three models
round(VAFw(dadMomChild, res1)*100, 1)
round(VAFw(dadMomChild, res2)*100, 1)
round(VAFw(dadMomChild, res3)*100, 1)

table1 <- cbind(res1$W, res2$W, res3$W)
table1 <- round(table1, 3)
rownames(table1) <- namesDadMomChild
colnames(table1) <- rep(c("T1", "T2", "T3"), 3)
table1


#save.image(file = './resPaper.RData')
#load('./resPaper.RData')

#Best model after LOOCV-error
allWs3comp[[bestModel3comp]]

#
plot(x = seq(0, 1, by = 0.01), lasso3comp[, 1],
     main = 'MPRESS given the lasso tuning parameter lambda',
     ylab = 'MPRESS',
     xlab = 'Lambda',
     type = 'l',
     xaxt = 'n' )

axis(side = 1, at = c(0, lasso3[bestLasso], 0.4, 0.6, 0.8, 1))
abline(v = lasso3[bestLasso], col = "gray", lty = "dashed")

table1 <- cbind(res1$W, res2$W, res3$W)
table1 <- round(table1, 3)
rownames(table1) <- namesDadMomChild
colnames(table1) <- rep(c("T1", "T2", "T3"), 3)
table1





