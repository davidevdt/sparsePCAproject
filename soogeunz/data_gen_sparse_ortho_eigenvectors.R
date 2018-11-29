# data generation through sparse and orthogonal eigenvectors = right singular vectors #

# this is a trick that allows specification of sparse and column-orthogonal matrix #

# (original) PCA always has column-orthogonal W and P which are equal to each other.
# with this trick we can make sparse W and P which fulfill that property.

# this procedure relies on SVD  #
# recall that PCA is the same as SVD. 
# SVD: X = UDV'
# PCA: X = XWP' = TP'

# from the above equations, we know that:
# V = W = P
# XW = T = UD

# we make use of the X = UDV' decomposition to specify our X
# we need to define U, D, V that meet the properties of SVD:
# U'U = I (left singular vectors)
# V'V = I (right singular vectors)
# D diagonal matrix with singular values

# importantly, we need to find V which is also sparse on top of column-orthogonality!

n <- 100
p <- 30
R <- 3 # number of components

# ** step 1. specify the sparse and orthogonal matrix V **
# (i named it V. it can either be W or P depending on your interpretation.)
# V = W = P
V <- matrix(0,p,R)
V[1:10,1] <- 1
V[11:20,2] <- 1
V[11:30,3] <- 1

V
# some common-dinstinctive processes in here

# not column-orthogonal
t(V) %*% V

# extract a block of the V matrix
# look at the rows that have elements in multiple columns
# in our V matrix, this is the rows 11:20
block <- V[11:20,2:3]

# orthogonalize this block
block_ortho <- qr.Q(qr(block))

# put this block back into the V matrix
V[11:20,2:3] <- block_ortho

t(V) %*% V
# you can observe that the off-diagonals are all zeros now

# this is a function that normalizes each column vector to unit vector
normalize <- function(MATRIX) {apply(MATRIX,2,function(x){x/norm(as.matrix(x),"F")})}

V <- normalize(V)

# V'V = I
t(V) %*% V


# ** step 2. specify the U and D **
# UD = T = XW
# randomly generate U from multivariate normal
# (columns of U are not correlated - remember that PCA components are uncorrelated)
set.seed(11)
U <- MASS::mvrnorm(n = 100, mu = rep(0,R), Sigma = diag(R), empirical=TRUE)

# orthogonalizing the U
# by doing this, we achieve U'U = I which is a property of left singular vectors
# (also means that our components are uncorrelated)
U <- qr.Q(qr(U))

t(U) %*% U

# now we specify the D matrix
# within SVD: this is the diagonal matrix with singular values
# this defines the amount of variance each corresponding principal component has
D <- diag(c(50, 30, 20))

# so now we have them all:
# V: V'V = I and sparse
# U: U'U = I 
# D 

X <- U %*% D %*% t(V)

svd1 <- svd(X)

round(svd1$v[,1:3] ,6)
# yes

round(prcomp(X)$rotation[,1:3], 6)
