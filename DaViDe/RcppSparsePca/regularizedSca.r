#########################################################
#####		regularizedSca : Lasso-Ridge SCA		#####
#########################################################

# First parse and compile the source code through : 
# Rcpp::sourceCpp( 'estimateSca.cpp' )

regularizedSca <- function(X, Q, LASSO = 0, RIDGE = 0, tol = 10^-8, maxIter = 250, fixW = NULL, nStart = 100, scaleData = TRUE, nScale = 0, showIt = FALSE, svdStart = TRUE){
	
	if( scaleData ){
		X <- scaleData(X, value = nScale)
	}
	
	I <- dim(X)[1]
	J <- dim(X)[2]
	XTX <- t(X) %*% X
	XTXdiag <- diag(XTX)
	
	if(length(LASSO) == 1){
		LASSO = rep(LASSO, Q)
	}
	
	if( svdStart ){
		W0  <- svd(X, Q, Q)$v	
		W0 <- W0[, 1:Q]		
	}else{
		W0 <-  matrix(rnorm(J * Q, sd = 0.5), J, Q)
	}

	W0 <- as.matrix(W0)

	
	if( !is.null(fixW)){
		W0[fixW == 0] <- 0
	}
	

	results <- estimateSca(	X, XTX, Q, LASSO, RIDGE, tol, maxIter, fixW, nStart, showIt, W0, I, J, XTXdiag )
		
	try(if( !results$converged ) 
		warning(' regularizedSca did not converge. Please re-run by increasing maxIter or the convergence criterion <tol> . ', call. = FALSE, immediate. = FALSE, noBreaks. = FALSE, domain = NULL) 				)
	
	
	return(results)
	
}
