
#################################################
##                                             ##
## Code for manuscript titled:                 ##
##                                             ##
## Network Approaches for Physical Activity    ##
## and Mental Health Research                  ##
##                                             ## 
##                                             ##  
## Christopher Huong: chrishuong95@gmail.com   ##
## Emily E. Bernstein                          ##
## Joshua E. Curtiss                           ##
## Denver M. Brown                             ##
##                                             ##
## 3/26/2024                                   ##
##                                             ##
#################################################


library(graphicalVAR)
library(qgraph)
library(corpcor)
library(mlVAR)



# Graphical VAR model from single-subject multivariate time series -------------------------------------------------

# Specify a hypothetical true data generating model. 
# Data are:
# Continuous
# Multivariate Gaussian
# Centered at 0

varnames <- c("Depr", "Anhed", "Fatigue", "PA")

# graphical VAR model: y[t] = Beta*y[t-1] + e[t] 
# where e[t] ~ N(0, Sigma)

# Sigma is the variance-covariance matrix after accounting for lag-1 temporal effects
# and is used to estimate the contemporaneous network (conditional dependencies within time points)
Sigma <- matrix(c(1.91,  0.69, 1.01,  0.87,
                  0.69,  1.39, 0.14, -0.14,
                  1.01,  0.14, 1.91,  1.22,
                  0.87, -0.14, 1.22,  1.91),nrow=4,byrow=T)

# Get PCC (partial contemporaneous correlations) by standardizing 
# the precision matrix (Kappa) and multiplying the off-diagonal elements by -1
Kappa <- solve(Sigma)
PCC <- -Kappa / (sqrt(diag(Kappa)) %*% t(sqrt(diag(Kappa))))
diag(PCC) <- 1
PCC
qgraph(PCC, layout="circle", labels=varnames)

# PA and fatigue are strongly associated within time points,
# patient interview may indicate that strenuous activity results
# in immediate feelings of fatigue.
# 
# A negative edge between PA and anhedonia may suggest
# anhedonia leads to sudden declines in PA in this individual
# 
# A positive edge between PA and depressed mood may arise
# from the individual using PA as a coping strategy when mood is low




# Beta encodes lagged regression effects
# To aid in interpretability we obtain the partial directed correlations
# by standardizing and transposing the beta matrix
# (thus rows indicate the node of origin, and columns the destination )
Beta <- matrix(c())





# PDC = partial directed correlations
# conditional dependencies between time points
PDC <- matrix(c(0.3,  0.02, -0.5,  0.3,
                0.0,  0.3,   0.02, 0.3,
                0.0, -0.5,   0.2,  0.5,
                0.0,  0.0,  -0.3,  0.3), nrow=4,byrow=T)

qgraph(PDC, layout="circle", labels=varnames)


# Simulate observations ---------------------------------------------------
# decomposing the grpahicalVARsim() function






# graphical VAR model: y[t] = Beta*y[t-1] + e[t] 
# where e[t] ~ N(0, Sigma)

# Sigma is the variance-covariance matrix after accounting for lag-1 temporal effects
# and is used to estimate the contemporaneous network
# from the partial correlation coefficients (conditional dependencies within time points)
Sigma <- matrix(c(1.7, 0.4, -0.3, -1.3,
                  0.4, 1.5, 0.5, -0.7,
                  -0.3, 0.5, 1.3, 0.5,
                  -1.3, -0.7, 0.5, 2.1), nrow=4,ncol=4, byrow=T)  

# Check if positive definite
is.positive.definite(Sigma)
# Kappa is thus the precision matrix (inverse of Sigma)
Kappa <- solve(Sigma)
# Get PCC (partial contemporaneous correlations) by standardizing Kappa and multiplying by -1
# PCC = - kappa[ij] / sqrt(kappa[ii] * kappa[jj])
PCC <- -Kappa / (sqrt(diag(Kappa)) %*% t(sqrt(diag(Kappa))))
diag(PCC) <- 1
PCC
qgraph(PCC, layout="circle", labels=varnames)

# Beta encodes lagged regression effects
# To aid in interpretability we obtain the partial directed correlations
# by standardizing and transposing the beta matrix
# (thus rows indicate the node of origin, and columns the destination )
PDC <- t(beta / sqrt( diag(solve(kappa)) %o% diag(kappa) + beta^2))



qgraph(PDC, layout="circle", labels=varnames)
qgraph(pcc, layout="circle", labels=varnames)

Beta <- matrix(c())



PCC <- matrix(c(1.0, 0.8, -0.5, 0.8,
                0.8, 1.0, 0.3, 0.05,
                -0.5, 0.3, 1.0, -0.5,
                0.8, 0.05, -0.5, 1.0), ncol=4,nrow=4,byrow=T)

PCC <- -Kappa / (sqrt(diag(Kappa)) %*% t(sqrt(diag(Kappa))))
diag(PCC) <- 1

qgraph(PCC, layout="circle", labels=varnames)


find_kappa <- function(PCC) {
  sqrt_diag <- sqrt(diag(PCC))
  return(-PCC * (sqrt_diag %*% t(sqrt_diag)))
}

# Find Kappa
Kappa <- find_kappa(PCC)

PCC <- -Kappa / (sqrt(diag(Kappa)) %*% t(sqrt(diag(Kappa))))
diag(PCC) <- 1

qgraph(PCC, layout="circle", labels=varnames)








PDC <- matrix(c(0.3,  0.02, -0.5,  0.3,
                0.0,  0.3,   0.02, 0.3,
                0.0, -0.5,   0.2,  0.5,
                0.0,  0.0,  -0.3,  0.3), nrow=4,byrow=T)



PDC <- t(beta / sqrt( diag(solve(kappa)) %o% diag(kappa) + beta^2))


# To compute the partial correlation coefficients

is.positive.definite()

eigen(randomsim$kappa)
is.positive.definite(randomsim$kappa)
eigen(randomsim$beta)
is.positive.definite(randomsim$beta)


ntime <- 50



















