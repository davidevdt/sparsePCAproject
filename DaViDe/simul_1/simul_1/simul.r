# Set your R environment to the 'simul_1' folder 

source('functions.r') 	# auxiliary functions
source('genData.r')		# The data generating function

# Parameters for the simulation
I = 10
J = 20
Q = 3			# 3 components example
startingX = NULL
sparsity = 0.1 
percNoise = 0.05
maxIter = 1				# Let's start with 1 iteration...
tol = 10^-5
seed = NULL
sdX = 1
sdNoise = 1
fixW = NULL
scaling = FALSE
seed = 42 		# " Let 42 be the Answer to the Ultimate Question of Life, the Universe, and Everything "
nScale = 0

gd <- genData( I = I, J = J, Q = Q, startingX = startingX, sparsity = sparsity, percNoise = percNoise, maxIter = maxIter, tol = tol, seed = seed, sdX = sdX, sdNoise = sdNoise, fixW = fixW, scaling = scaling, nScale = nscale )

gd$W
gd$W_found
gd$convergence			

# As you can see, the 'convergence' value is not too bad after only 1 iteration, but above the 'tol' paramter we set for tolerance. Also, the sought W and W_found don't seem to differ too much in term of squared differences (the 0's in the goal matrix are values more or less close to 0 in the W matrix estimated before inserting noise). So let's see how things evolve after a bunch more iterations...


maxIter = 100

gd <- genData( I = I, J = J, Q = Q, startingX = startingX, sparsity = sparsity, percNoise = percNoise, maxIter = maxIter, tol = tol, seed = seed, sdX = sdX, sdNoise = sdNoise, fixW = fixW, scaling = scaling, nScale = nscale )

gd$W
gd$W_found
gd$convergence	

# Yes, convergence criterion is now worse. Furthermore, some terms that should be close to 0 in the estimated W (W_found) are now farther from 0 than after 1 iteration. What I found more interesting, though, is that the values in the XTX matrix (as a result of an increase in the values of the X matrix) [I put the print-screen of the diag of XTX to ease visualization.	]. Look when we increase the iterations to 1000 : 


maxIter = 1000

gd <- genData( I = I, J = J, Q = Q, startingX = startingX, sparsity = sparsity, percNoise = percNoise, maxIter = maxIter, tol = tol, seed = seed, sdX = sdX, sdNoise = sdNoise, fixW = fixW, scaling = scaling, nScale = nscale )

gd$W
gd$W_found
gd$convergence	


# The convergence is exploding, as well as the values of XTX and the ones of X :

gd$X 

# ...the hell is this? 
# Eventually, as the number of iterations goes to infinity, the probability for the svd() function to calculate the SVD of XTX*W decreases : 

maxIter = 1000000


gd <- genData( I = I, J = J, Q = Q, startingX = startingX, sparsity = sparsity, percNoise = percNoise, maxIter = maxIter, tol = tol, seed = seed, sdX = sdX, sdNoise = sdNoise, fixW = fixW, scaling = scaling, nScale = nscale )

# Furthermore, as you can see from the print-screens, the squared diefferences don't seem to increase monotonically. 


# So my question is guys : what am I doing wrong? Do you think the squared differences between the W_target and the W_found matrices are the right thing to do? Even if I change the criterion (e.g., Tucker congruence), this would not change the way how the X is updated. Am I doing something idiot or missing something? Can you give a look at my code to see if I am doing it in the right way? 
