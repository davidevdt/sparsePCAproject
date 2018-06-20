####################################
### Scaling (borrowed from Niek) ###
####################################

scaleData <- function(X, value = 0){

	X <- scale(X, scale = FALSE)   
	
	#attr(X, "scaled:center") <- NULL 
	
	sdX <-  apply(X, 2, function(x) sqrt(sum( x^2, na.rm = TRUE ) / (length(na.omit(x)) - value )))  	
	sdX <- matrix(sdX, nrow(X), ncol(X), byrow = T)          
	
	X <- X * (1 / sdX)      #divide each entry in X by its sd
	
	return(X)
}


###############################################
### Tucker Congruence (borrowed from Niek)	###
###############################################
# New part : diffMode = TRUE, so that the function can calculate SS of differences between matrices. 


tuckerCongruence <- function(A, B, diffMode = FALSE){                                                                 
    combinationList <- combinat::permn(1:(ncol(A)))
	fixA <- A                                                                   
    changeInSign <- c(1, -1)
    placeHolder <- rep(NA, 2)
	
	if( diffMode == FALSE ){
		res <- -Inf 
        for(combination in combinationList){		
			A <- fixA[, combination]
			tucCon <- rep(NA, ncol(A))
			for( i in 1:ncol(A) ){
                for(j in c(1, 2)){
					placeHolder[j] <- (changeInSign[j]*A[, i]) %*% B[, i]  / sqrt(sum(A[, i]^2) * sum(B[, i]^2)) 
                }                                 
                tucCon[i] <- max(placeHolder)
            }
			candidate <- mean(tucCon)
			
            if(candidate > res){
				res <- candidate
				perm <- combination 
            }        
			
        }      
		return( list(Tucker = res, perm = perm) )
	}else{  	# I allow the function to do the same, but with sum of squared differences between closest columns 
		res <- Inf 
		for(combination in combinationList){
			A <- fixA[, combination]
			diffs <- rep(NA, ncol(A))
			#signs <- rep(NA, ncol(A))
			for( i in 1:ncol(A) ){
                for(j in c(1, 2)){
					placeHolder[j] <- sum(((changeInSign[j]*A[, i]) - B[, i])^2)  
                }                                 
                diffs[i] <- min(placeHolder)
				#signs[i] <- changeInSign[which(placeHolder == min( placeHolder ))]
            }
			candidate <- sum(diffs)
			
			if( candidate < res ){
				res <- candidate 
				perm <- combination
			}			
		}
		return( list(sumSquares = res, perm = perm) )
	}
} 






######################
### SS of a matrix ###
######################

varMatrix <- function(A){
	V = sum( A^2 )
	return( V )
}