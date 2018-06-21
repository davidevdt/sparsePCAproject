#############################################
###	eigenCV : EigenVector Cross-Validation###
#############################################
source('./Functions.R')
source('./lrSca.R')

eigenCV <- function(X, Q, RIDGE, LASSO, fixW, nFolds, nScale = 0, tol = 10^-8, maxIter = 250, nStart = 100, showIt = FALSE){
	
	# Dataset permutation and Folds Creation
	X <- X[sample( 1:nrow(X), nrow(X), rep = FALSE ) ,]
	folds <- rep_len( 1:nFolds, nrow(X) )
	
	cvError <- matrix(NA, nrow(X), ncol(X))
	
	meanErr <- rep(0, nFolds)

    for(n in 1:nFolds){
	
		if(showIt){
			cat('Processing fold # ', n, '\n')
		}

        # actual split of the data
        fold <- which(folds == n)
        trainDat <- X[ -fold, ]
		testDat <- X[ fold, , drop = FALSE ]
		
		#Scale training and the testing data
		means <- apply( trainDat, 2, mean )		
		trainDat <- t( t(trainDat) - means)
		sdTrain <- apply(trainDat, 2, function(x) sqrt(sum( x^2, na.rm = TRUE ) / (length(na.omit(x)) - nScale )))
		trainDat <- t( t(trainDat)/sdTrain )
		
		testDat <-  t( t(testDat) - means)
		testDat <- t( t(testDat)/sdTrain ) 
		
		rm(means, sdTrain)		
        
		# Model estimation
		res <- lrSca(trainDat, Q, LASSO = LASSO, RIDGE = RIDGE, tol = tol, maxIter = maxIter, fixW = fixW, nStart = nStart, scaleData = FALSE, showIt = FALSE)
		
        #Eigenvector crossvalidation Bro Kiers
        pred <- matrix( NA, nrow(testDat), ncol(testDat) )
		
        
        for( j in 1:ncol(X) ){
		
            TMinusVariableJ <- testDat[, -j] %*% res$W[-j, ]
            pred[, j] <- TMinusVariableJ %*% res$P[ j, ]
			
		}
		
        cvError[fold, ] <- (testDat - pred)^2
		meanErr[n] <- mean(cvError[fold, ])
	}
	
	return( list( cvError = cvError, press = sum(cvError), mpress = mean(cvError), sdMpress = sd(meanErr) ))

}
