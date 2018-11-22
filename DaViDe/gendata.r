# The function generates a 3-components structures with 3 blocks of variables, 2 of which are correlated with each other on one component. 
genStructuredData <- function( n, p1, p2, p3, propNoise, vars, prop3, coefs1 = NULL, coefs2 = NULL, coefs3 = NULL ){
	p <- p1+p2+p3
	
	if( is.null(coefs1) ){
		coefs1 <- c(runif( p1,-1,1 ),rep(0,p2+p3))
	}
	if( is.null(coefs2) ){
		coefs2 <- c(rep(0,p1),runif( p2,-1,1 ),rep(0,p3))
	}
	if( is.null(coefs3) ){
		coefs3 <- c(rep(0,p1),runif( p2+p3,-1,1 ))
	}	
	
	coefs <- cbind(coefs1,coefs2,coefs3)

	# 1st part : create Xsignal 
	V1 <- rnorm(n,0,sqrt(vars[1]))
	V2 <- rnorm(n,0,sqrt(vars[2]))
	V3 <- rnorm(n,0,sqrt(vars[3]))
	Xsignal <- matrix(0,n,p)


	Xsignal[,1:p1] <- V1
	Xsignal[,(p1+1):(p1+p2)] <- V2
	Xsignal[,(p1+p2+1):p] <- V3

	TmpMat1 <- t(t(Xsignal[,1:p1]) * coefs1[1:p1])
	TmpMat2 <- (1-prop3)*(t(t(Xsignal[,(p1+1):(p1+p2)]) * coefs3[(p1+1):(p1+p2)])) + (prop3)*(t(t(Xsignal[,(p1+p2+1):p]) * coefs3[(p1+p2+1):p]))
	TmpMat3 <- t(t(Xsignal[,(p1+1):(p1+p2)]) * coefs2[(p1+1):(p1+p2)])


	Xsignal <- cbind(TmpMat1,TmpMat2,TmpMat3)
	s = svd( Xsignal )
	#round(s[[1]],2)
	#round(s[[3]],2)[,1:3]				# Check out W-matrix Structure

	## 2nd part : artificially round to 0 smallest values from each component
	W <- s[[3]][,1:3]
	W[(p1+1):p,1] <- 0
	W[1:p1,2] <- 0
	W[(p1+p2+1):(p),2] <- 0
	W[1:p1,3] <- 0


	# 3rd part: Re-construct Xsignal with new W matrix 
	P <- svd((t(Xsignal)%*%Xsignal)%*%W)[[2]]			#Procrustes rotation
	Xnew <- Xsignal %*% W %*% t(P)
	#round(svd(Xnew)[[1]],2)
	#round(svd(Xnew)[[3]],2)[,1:3]			
	Wnew <- svd(Xnew)[[3]][,1:3]

	# 4th part : Add Noise 
	E <- matrix( rnorm(n*p,0,1), n, p )
	g = sqrt((var(as.vector(Xnew))*propNoise) / (var(as.vector(E)) * (1 - propNoise)))
	Xfinal <- Xnew + E*g

	SStrue <- var(as.vector(Xnew))
	SSX <- var(as.vector(Xfinal))


	res <- list( X=Xfinal, W=Wnew, P=P, errorRatio = 1 - (SStrue/SSX)  )
	return( res )
}


# Try out 
#set.seed(42)
#set.seed(runif(1))
#n = 1e+2					# The larger the sample size, the closer the wanted values will go to 0
#p1 <- 20
#p2 <- 20
#p3 <- 20 


#coefs1 <- c(runif( p1,-1,1 ),rep(0,p2+p3))
#coefs2 <- c(rep(0,p1),runif( p2,-1,1 ),rep(0,p3))
#coefs3 <- c(rep(0,p1),runif( p2+p3,-1,1 ))

#propNoise <- 0.1
#vars <- c(5,0.1,3)
#prop3 <- 0.5				# Play around with this parameter and see how W changes

#genX <- genStructuredData( n, p1, p2, p3, propNoise, vars, prop3  )
#round( svd(res$X)[[1]] ,2)
#round( svd(res$X)[[3]][,1:3] ,2)
#e <- eigen( t(res$X)%*%res$X )
#plot( e[[1]], type = "l" )
