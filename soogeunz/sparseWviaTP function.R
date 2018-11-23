# sparseWviaTP function #

# this is a function that takes the Davide's idea to generate the TP' first,
# then work on sparse weights

# As Davide's idea is equal to generating data through sparse P, 
# I have generalized the function such that it allows any numbers of variables or components

# The function involves iteration. Here is how it works:
  # input:
  # n = sample size
  # R = number of components
  # Wstructure = matrix that specifies the zero-non-zero pattern of the W
  # randomly = TRUE / FALSE. if TRUE, generates random values from uniform dist.
  # to replace the non-zero elements provided in Wstructure matrix
  # eigs = vector that specifies eigenvalues
  # seed = seed
  # maxiter = number of iterations

# steps:
  # 1. Generate T matrix by sampling from normal distribution and orthogonalizing
  # 2. Define a sparse P structure by using Wstructure
  # 3. with T, P and eigenvalues, we now have X = TP'
  # (eigenvalues specify the amount of variance in each component)

# what do we want at this point? We want either:
  # (A) X matrix that fulfills abs(X - XWP' ) == 0, while having W sparse and P column-orthogonal
  # (B) X matrix that actually has eigenvectors with zero-elements. Note that eigenvectors = W = P in standard PCA
  # we try to achieve either of these by the following steps.

# back to the steps:
  # 3. with T, P and eigenvalues, we now have X1 = TP'

  # for iteration, we now define Xiter <- X1. 
    # ITERATE {
  # 4. calculate the right singular vectors "V" from Xiter = UDV'. save the squared error between the elements of the V, at the indices where your pre-defined Wstructue is zero.  
  # 5. apply the zero-pattern provided by Wstructure on V. At indices where Wstructure is zero, force V elements to be zero. This sparse V matrix is now called "Wsparse"
  # 6. now we have X1 and Wsparse. Using these, calculate P by Procrustes rotation. This P is now called "Pnew".
  # 7. calculate the squared error of (X - X1 * Wsparse * Pnew'). save this squared error.
  # 8. define "Xnew" <- X1 * Wsparse * Pnew'. This is your new dataset
  # 9. define Xiter <- Xnew, and do the steps 4 to 8 again, each time saving the squared error at step 5 and step 7.
    # }
  # 10. this process does not converge. However, after many iterates like 10000, you can find which iterate led to a very small error, and use that dataset.
    
# you can try the function by using example codes below the function definition

normalize <- function(MATRIX) {apply(MATRIX,2,function(x){x/norm(as.matrix(x),"F")})}

sparseWviaTP <- function(n, R, Wstructure, randomly, eigs, seed, maxiter){
  
  # you can just provide the zero-nonzero pattern
  # the function replaces the nonzeros with values sampled from uniform dist
  if (randomly){
    set.seed(seed)
    nr <- sum(Wstructure!=0)
    nonzeros <- runif(n = nr, min = -1, max = 1)
    Wstructure[Wstructure!=0] <- nonzeros
  } 
  
  zeropattern <- (Wstructure == 0)
  
  P <- Wstructure
  
  # specifying the T matrix (norm 1: U matrix)
  Tmat <- MASS::mvrnorm(n = n, mu = rep(0,R), Sigma = diag(R))
  
  # orthogonalizing the Tmat
  Tmat <- qr.Q(qr(Tmat))
  
  # eigenvalues (vector)
  S <- eigs
  
  i <- 1
  V_hist <- list()
  zero_hist <- c()
  disc_hist <- c()
  W_hist <- list()
  P_hist <- list()
  X_hist <- list()
  
  # X1 generated
  X1 <- Tmat %*% diag(S) %*% t(P)
  
  X_hist[[i]] <- X1
  
  # the initial run
  svd1 <- svd(X1)
  V1 <- svd1$v[,1:R]
  
  V_hist[[i]] <- V1
  
  # the initial difference with the zero-elements
  zero1 <- sum(abs(V1[zeropattern]))
  
  zero_hist[i] <- zero1
  
  # forcing the zero-pattern, creating sparse W
  W1 <- V1
  W1[zeropattern] <- 0
  
  W_hist[[i]] <- W1
  
  # normalize W1: 
  # due to forcing near-zero values to zero, columns have lengths different from 1
  W1 <- normalize(W1)
  
  # re-calculate the P through Procrustes rotation (using sparse W)
  svdp1 <- svd(t(X1) %*% X1 %*% W1)
  P1 <- svdp1$u %*% t(svdp1$v)
  
  # the initial discrepancy (X - XWP')
  disc1 <- sum((X1 - X1 %*% W1 %*% t(P1))^2)
  
  disc_hist[i] <- disc1
  
  # then we have our new X to start working with..
  X2 <- X1 %*% W1 %*% t(P1)
  
  while (i < maxiter){
    i <- i + 1
    
    X_hist[[i]] <- X2
    
    # check if the zero-elements of the new eigenvectors are closer than X1?
    svd2 <- svd(X2)
    V2 <- svd2$v[,1:R]
    
    V_hist[[i]] <- V2
    
    # difference with the zero-elements?
    zero2 <- sum(abs(V2[zeropattern]))
    
    zero_hist[i] <- zero2
    
    # forcing the zero-pattern on eigenvector, making sparse W
    W2 <- V2
    
    W2[zeropattern] <- 0
    
    # normalize W: 
    W2 <- normalize(W2)
    
    W_hist[[i]] <- W2
    
    # re-calculate the P through Procrustes rotation
    svdp2 <- svd(t(X2) %*% X2 %*% W2)
    P2 <- svdp2$u %*% t(svdp2$v)
    
    P_hist[[i]] <- P2
    
    # (X-XWP')^2 discrepancy
    disc2 <- sum((X2 - X2 %*% W2 %*% t(P2))^2)
    
    disc_hist[i] <- disc2
    
    # new X for the next iteration
    X3 <- X2 %*% W2 %*% t(P2)
    
    X2 <- X3
    
  }
  result <- list(X_hist = X_hist, 
                 W_hist = W_hist,
                 V_hist = V_hist,
                 P_hist = P_hist,
                 zero_hist = zero_hist,
                 disc_hist = disc_hist)
  return (result)
}

# test 1: structure that we already talked about
Wstructure1 <- matrix(0,30,3)
Wstructure1[1:10,1] <- 1
Wstructure1[11:20,2] <- 1
Wstructure1[11:30,3] <- 1
# you can observe that the third component is common, while the first two are distinctive

hii <- sparseWviaTP(n = 100, R = 3, Wstructure = Wstructure1, randomly = TRUE, eigs = c(50,30,20), seed = 11, maxiter = 10000)

hii$disc_hist
# squared distance (X - XWP') at each iterate

hii$zero_hist
# squared distance (V[zero pattern] - Wstructure[zero pattern]) at each iterate

min(hii$disc_hist)
min(hii$zero_hist)

disc_min <- which.min(hii$disc_hist)
zero_min <- which.min(hii$zero_hist)

dat_disc <- hii$X_hist[[disc_min]]
W_disc <- hii$W_hist[[disc_min]]
P_disc <- hii$P_hist[[disc_min]]

# this part is really important. 
# this informs us about the scaling of our X matrix
# often, you might be able to find the X matrix that has a really really small squared error of (X - XWP'),
# but that might be due to the fact that your X is simply scaled very small.
# Here, you can observe that this X is definitely scaled similar to the original scaling
sum(apply(dat_disc^2,2,sum))
sum(apply(hii$X_hist[[1]]^2,2,sum))

svd_disc <- svd(dat_disc)
round(svd_disc$v[,1:3],10)
svd_disc$d[1:4]

sum((dat_disc - dat_disc %*% W_disc %*% t(P_disc))^2)
# squared error

W_disc
# sparse W

P_disc
# P very similar to W

# using the zero-distance criterion..
dat_zero <- hii$X_hist[[zero_min]]
W_zero <- hii$W_hist[[zero_min]]
P_zero <- hii$P_hist[[zero_min]]

sum(apply(dat_zero^2,2,sum))
sum(apply(hii$X_hist[[1]]^2,2,sum))

svd_zero <- svd(dat_zero)
round(svd_zero$v[,1:3],10)
svd_zero$d

sum((dat_zero - dat_zero %*% W_zero %*% t(P_zero))^2)


# test 2: a different structure
Wstructure2 <- matrix(0,30,3)
Wstructure2[1:10,1] <- 1
Wstructure2[11:30,2] <- 1
Wstructure2[11:20,3] <- 1
# this time the second component is common


hii2 <- sparseWviaTP(n = 100, R = 3, Wstructure = Wstructure2, randomly = TRUE, eigs = c(50,30,20), seed = 11, maxiter = 10000)

hii2$disc_hist
# squared distance (X - XWP') at each iterate

hii2$zero_hist
# squared distance (V[zero pattern] - Wstructure[zero pattern]) at each iterate

min(hii2$disc_hist)
min(hii2$zero_hist)

disc_min <- which.min(hii2$disc_hist)
zero_min <- which.min(hii2$zero_hist)

dat_disc <- hii2$X_hist[[disc_min]]
W_disc <- hii2$W_hist[[disc_min]]
P_disc <- hii2$P_hist[[disc_min]]

# Observe here.. 
sum(apply(dat_disc^2,2,sum))

apply(dat_disc^2,2,sum)
# each column has become really really small

head(dat_disc[,1:5])

# this is the sum of squares of the initial X specification
sum(apply(hii2$X_hist[[1]]^2,2,sum))

svd_disc <- svd(dat_disc)
round(svd_disc$v[,1:3],10)

# eigenvalues also very small
svd_disc$d[1:4]

sum((dat_disc - dat_disc %*% W_disc %*% t(P_disc))^2)
# squared error

W_disc
# sparse W

P_disc
# P very similar to W

# therefore, although the X - XWP' discrepancy has reduced, this X is not feasible

# basically, i think that this function (or a sparse W via X = TP' maneuver) only works for certain structures only
