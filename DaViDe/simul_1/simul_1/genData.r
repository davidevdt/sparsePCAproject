#####################################################################
#####		Dataset Creation according to a specific W matrix	#####
#####################################################################

# The function takes as input : 

# I, number of subjetxt
#J, number of variables
#startingX , an optional initial data matrix;
# sparsity, the % of desired sparsity 
# percNoise, the desired level of noise; 
# maxIter, the number of iterations for the re-calculation of X and P
# tol, the criterion for stopping the iterations (based on sum(W_target - W[i])^2 at iteration i) 
# seed, set the seed for random number generators
# sdX and sdNoise, standard deviations for generating X and noise
# fixW , a pattern matrix (it matters only when it is equal to 0)
# scaling, nScale : should the matrix X be scaled? If so, with what correcting factor for the sd? 


# Output : 
# X, the final X matrix; 
# W, the final W matrix 
# P, the final P matrix
# noiseLevel, the generated % of noise on X  

# Here I also leave W_found (the W matrix estimated by svd before adding noise) as an output, just for debugging  


genData <- function( I, J, Q, startingX = NULL, sparsity = 0, percNoise = 0.05, maxIter = 20, tol = 10^-8, seed = NULL, sdX = 1, sdNoise = 1, fixW = NULL, scaling = FALSE, nScale = 0, showIt = TRUE ){
	
	
	# Set seed 
	if( !is.null(seed) ){
		set.seed( seed )
	}
	
	### Generate X 
	if( is.null(startingX) ){
	
		# If a starting matrix is not provided generate data randomly
		X0 <- matrix(rnorm( (I * J) , 0, sd = sdX ), I, J)
		
	}else{
	
		X0 <- startingX
	
	}

	if( scaling == TRUE ){
		X0 <- scaleData(X0, value = nScale)
	}

	 

	# Initialize W_target 
	SVD_1 <- svd(X0, Q, Q)
	W_target <- SVD_1$v
		
	rm( SVD_1 )

	if( !is.null(fixW) ){
		W_target[ fixW == 0 ] <- 0
	}
	
	# Insert sparsity (Niek's method)
	if( sparsity > 0 ){
	
		if( is.null(fixW) ){
		
			W_target[ abs(W_target) < quantile(abs(W_target), sparsity) ] <- 0
		
		}else{
		
			W_target[ abs(W_target[ W_target!= 0]) < quantile(abs(W_target[ W_target!= 0]), sparsity) ] <- 0
		
		}
	
	}


	
	sumDiff <- 0
	
	### Iterate and generate X0 and P0 until the calculated W and W_target are close to each other

	for( i in 1:maxIter){
	
		if(showIt & (i %% 10) == 0){
			cat( 'Iteration ', i , ' , diff =  ', sumDiff,'\n' )
		} 
		
		
		XTX <- t(X0) %*% X0
		SVD_2 <- svd(XTX %*% W_target)
		P0 <- SVD_2$u %*% t(SVD_2$v)
		
		
		#######		QuEsTioN	#####
		# Why are the values of X0, and thus the ones of XTX,  increasing at each iteration? 
		cat('Diag(XTX) = ', diag(XTX), '\n')
		
		rm( SVD_2 )	
		
		# Update X0
		X0 <- X0 %*% W_target %*% t(P0)

		SVD_1 <- svd(X0, Q, Q)
		W_new <- SVD_1$v
		rm( SVD_1 )
		
		
		# Should I leave this one? I think so 
		if( !is.null(fixW) ){
			W_new[ fixW == 0 ] <- 0
		}
		
		# Calculate SS of differences between W_new and W_target (modified tuckerCongruence function)  
		sumDiff <- tuckerCongruence(W_new, W_target, diffMode = TRUE)$sumSquares		
		# When the sum of differences is below the set tol, we can break the for loop
		if( sumDiff < tol ){
			break
		}
		
		
	}

	try( if( i >= maxIter ) 
		message( 'Warning : convergence not reached after ', maxIter, ' iterations.' )	)
		
	
	Xtrue <- X0 
	rm(X0)
	
	### Add Noise ###
	
	E <- matrix(rnorm( (I * J), 0, sd = sdNoise), I, J )
	
	SSX <- varMatrix( scale(Xtrue, scale = FALSE ) )
	SSN <- varMatrix( scale(E, scale = FALSE ) )
	
	# Adjusting factor
	if( percNoise == 1 ){
		percNoise <- .999999999999999	# To allow calculation of k
	}
	k <- sqrt( (percNoise * SSX) / ( SSN - ( percNoise * SSN ) ) )
	
	# Add the adjusted noise to the signal
	Xtrue <- Xtrue + ( k * E )
	
	SSX <- varMatrix( scale(Xtrue, scale = FALSE ) )	
	SSN <- varMatrix( scale( k * E, scale = FALSE ) )
	
	# Final proportion of noise generated 
	propNoiseGen <- SSN / SSX 
	
	return( list( X = Xtrue, W = W_target, P = P0, noiseLevel = propNoiseGen, W_found = W_new, convergence = sumDiff ) )
	
}















