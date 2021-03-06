#################################################################
# Function to create data according to sparse weights loadings stucture
# X = XPt(P)
# with common and distinctive sparse P
##################################################################
library(MASS)
library(gtools)

#divide the columns of matrix by 2norm
divideByNorm <- function(mat){
    A <- 1 / matrix(apply(mat, 2, function(x){sqrt(sum(x^2))}),
                  nrow(mat), ncol(mat), byrow = T)
    return(mat * A)
}

#function to get all common and specific block structures
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

#determine the row indices of groups
#probably I do not need this function anymore
determineGroups <- function(groups){
	out <- list()
	add <- c(0, cumsum(groups)[-length(groups)]) 
	for(i in 1:length(groups)){
		out[[i]] <- 1:groups[i] + add[i]
	}
	return(out)
}

#Start of the data generating functions
#Functions to create data
################ This is where the magic happens ####################

#orthogonalize two columns of a matrix by using gramm-schmid
#only use intersection of the non-zero coefficients of the two vectors
#to enforce that the zero's do not change
orthogonalize <- function(A, index1, index2){
    int <- intersect(which(A[, index1] != 0), which(A[, index2] != 0))

    u <- A[int, index1]
    v <- A[int, index2]

    newv <- v - as.numeric((t(u) %*% v) / (t(u) %*% u)) * u
    A[int, index2] <- newv

    return(A)
}

#create an orthonormal matrix 
#where the zero's are in fixed position
makeP <- function(A){
    counter <- 0
    #while matrix A is not orthonormal and max iterations is not reached
    #do orthogonalize all columns with respect to each other
    #Do this until the columns are orthogonal, 
    #sometimes multiple passes are needed
    while(TRUE != all.equal(t(A) %*% A, diag(ncol(A))) &&  counter < 1000 ){
        for(i in 2:ncol(A)){
            for(j in 1:(i-1)){
                A <- orthogonalize(A, j, i)
            }
        }
        A <- divideByNorm(A)
        counter  <- counter + 1
        print(counter)
    }
    if(counter < 1000){
        return(list(A=A, status=1))
    } else {
        return(list(A=A, status=0))
    }
}

#make a fixed percentage of the non zero coefficients zero (at random) 
sparsify <- function(comdis, sparsity){
    amounts <- round(apply(comdis, 2, function(x){sum(x != 0)}) * sparsity)
    TF <- apply(comdis, 2, function(x){x != 0})
    where <- apply(TF, 2, which)
    for(i in 1:length(where)){
        comdis[sample(where[[i]], amounts[i]), i]  <- 0
    }
    return(comdis)
}

#create a dataset
makeDat <- function(n, p, ncomp, comdis, variances){
    #generate random P and fix the zero structure
    #and make P square
    P <- matrix(rnorm(p*ncomp), p, ncomp)
    P[comdis == 0] <- 0
    P <- cbind(P, matrix(rnorm(p*(p-ncomp)), p, p-ncomp))

    result <- makeP(P)
    if(result$status == 1){
        P <- result$A
        Sigma <- P %*% diag(variances) %*% t(P)
        X <- mvrnorm(n, mu = rep(0, p), Sigma, empirical=F)
        return(list(X=X, P=P, Sigma=Sigma))
    } else {
        print("failed")
        return(NA)
    }
}

#function to create variances for the components
#you have to supply the variances of the components you are interested in
#the variances of the other non interesting components are on a log scale 
#starting with mean(variances of interesting components) /2
#these variances then get scaled such that error ratio is gotten. 
makeVariance <- function(varianceOfComps, p, error){
    ncomp <- length(varianceOfComps)
    varsOfUnimportantComps <- exp(seq(log(0.0001), log(mean(varianceOfComps)/2),
                                      length.out = p-ncomp))[(p-ncomp):1]

    x <- (-error*sum(varianceOfComps) / (error-1)) / sum(varsOfUnimportantComps)
    
    return(c(varianceOfComps, x * varsOfUnimportantComps))
}
