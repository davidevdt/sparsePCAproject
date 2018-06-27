/************************************
*************************************
*		estimateSca - Inner Cycle  	*
*************************************
************************************/ 

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]



// [[Rcpp::export]]
Rcpp::List estimateSca(	arma::mat X, arma::mat XTX, int Q, arma::vec LASSO, double RIDGE, double tol, int maxIter, 
					Rcpp::Nullable<Rcpp::NumericMatrix> fixW, int nStart, bool showIt, arma::mat W0, int I, int J, 
					arma::vec XTXdiag ){
	
	// RNG initialization  - I need to understand better how to set this part
	//arma::arma_rng::set_seed_random();
	//Rcpp::Environment base_env("package:base");
    	//Rcpp::Function set_seed_r = base_env["set.seed"];
    	//set_seed_r(std::floor(std::fabs(seed)));

	
	// Initialize objects
	Rcpp::List ret; 
	
	arma::mat P; 	
	arma::mat W; 

	arma::vec c;  	
	arma::vec lossValue; 
	arma::vec tmp(2);
	tmp.fill( 0.0 );
	arma::vec loss; 
		
	double w ; 
	double w_new; 
	double cp; 
	double minLoss; 
	double minLossGlobal = arma::datum::inf; 
	int sgn; 
	int n; 
	int itCounter;
	int j; 
	int q; 
	bool converged; 
	
	
	// Multi-Start Cycle 
	for( n = 0; n<nStart; n++){
		
		lossValue = arma::vec( maxIter + 1);
		lossValue.fill( arma::datum::inf ); 
		converged = FALSE; 
		W = W0 + (arma::randn(J, Q) * 0.001); 
		
		if( fixW.isNotNull( ) ){			
			W.elem( find(Rcpp::as<arma::mat>(fixW) == 0.0) ).zeros();			
		}
		
		for( itCounter = 1; itCounter < (maxIter + 1); itCounter++ ){
			
			// Update P :
			arma::mat U; 
			arma::vec D; 
			arma::mat V; 
			arma::svd(U, D, V, (XTX * W));
				
			U = U.cols(0, (Q-1)); 	
				
			P = U * V.t();
			
			
			// Update W via coordinate descent :
			for( q = 0; q < Q; q++ ){
				
				c =  X * ( P.col(q) - W.col(q) );

				for( j = 0; j < J; j++ ){
					
					if( W0(j, q) != 0.0 ){
						
						w = W( j, q );
						cp = (1 / double(I)) * (dot((X.col(j)), c.t()) + (w * XTXdiag(j)));	 
						
						if( cp > 0 ){
							sgn = 1; 
						}else if( cp < 0 ){
							sgn = -1; 
						}else{
							sgn = 0; 
						}
						
						tmp(0) = std::abs(cp) - LASSO(q);						
						w_new = double(sgn) * max( tmp );
						w_new /= ( ((1 / double(I)) * XTXdiag(j)) + RIDGE );

						if( abs(cp) < LASSO(q)){
							c += (w * X.col(j) );
						}else{
							c += ((w - w_new) * X.col(j));
						}
						
						W(j, q) = w_new; 
					}					
				}				
			}		
			// End coordinate descent 
			
			// Update loss function 
			lossValue( itCounter ) = ( 1 / double(I * 2) ) * accu( pow( X - ( X * W * P.t() ) , 2 ) ) + ( (1 / double(2)) * RIDGE * accu( pow( W, 2) ) ) + dot(LASSO, sum( abs(W), 0 )); 
			
			
			// Evaluate condition to stop inner loop 
			if( lossValue( itCounter - 1 ) - lossValue( itCounter ) <= tol ){
				converged = TRUE;
				break; 
			}
			
		}
		
		loss = lossValue( find(lossValue != arma::datum::inf) );
		minLoss = loss( loss.size() - 1 ); 
		
		
		if( showIt ){ 
			if( converged ){
				Rcpp::Rcout << "Start # " << n + 1 << " has converged in " << itCounter << " iterations; loss = " << minLoss << "\n" ; 
			}else{
				Rcpp::Rcout << "Start # " << n + 1 << " has not converged after " << itCounter << " iterations. \n" ; 
			}
		}
		
		
		if( minLoss < minLossGlobal ){
			
			minLossGlobal = minLoss; 
			
			ret["W"] = W; 
			ret["P"] = P; 
			ret["loss"] = minLoss; 
			ret["converged"] = converged; 		
		}	
	
	}
	
	return ret; 	
	
}
