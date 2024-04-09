library(mlVAR)
library(graphicalVAR)
library(qgraph)
library(dplyr)
library(MASS)
library(corpcor)

# Simulate N > 1 multivariate timeseries data
# where there are between-person differences in stationary means,
# but not in temporal lagged-effects


nVar <- 4 
nTime <- 100
nPerson <- 30
varnames <- c("Depr", "Stress", "Fatigue", "PA")

# Specify lagged relationships
Phi <- matrix(c(0.3, 0.0, 0.0, 0.0,
                0.6, 0.0, 0.0, 0.0,
                0.6, 0.3, 0.3,-0.6,
                -0.6, 0.0,-0.6, 0.0), nrow=4,byrow=T)

truePDCnetwork <- qgraph(Phi, layout="circle",
                         theme="colorblind", labels=varnames)


eigen(Phi)$values

# Specify a contemporaneous network
truePCC <- matrix(c(1.0, 0.0, 0.6, 0.0,
                    0.0, 1.0, 0.3, 0.4,
                    0.6, 0.3, 1.0,-0.4,
                    0.0, 0.4,-0.4, 1.0), nrow=4, byrow=T)

Sigma <- pcor2cor(truePCC)
# True partial contemporaneous correlations (PCC)
# Sigma is inverted (to obtain the precision matrix) and standardized 
# to produce the partial correlation matrix
truePCCnetwork <- qgraph(Sigma, graph="pcor",
                         layout="circle", theme="colorblind", labels=varnames)

# Check if positive-definite (eigenvalues > 0)
eigen(Sigma)$values > 0

# Specify a variance-covariance matrix for between-person differences
trueBetween <- matrix(c(1.0,  0.6,  0.0, -0.4,
                        0.6,  1.0,  0.5,  0.0,
                        0.0,  0.5,  1.0, -0.4,
                       -0.4,  0.0, -0.4,  1.0), nrow=4, byrow=T)

Omega <- pcor2cor(trueBetween) 
trueBetweenNetwork <- qgraph(Omega, graph="pcor",
                             layout="circle", theme="colorblind", labels=varnames)


# Initialize empty data frames for innovations
U <- list()
for(i in 1:nPerson){
  U[[i]] <- matrix(0, nrow=nTime, ncol=nVar)
}
# Simulate innovations (dynamic errors) from multivariate Gaussian distribution 
# covariance matrix Sigma

for(i in 1:nPerson){
  U[[i]] <- mvrnorm(n=nTime, mu=rep(0, nVar), Sigma=Sigma)
}

# Initialize person specific data frames
Data <- list()
for(i in 1:nPerson){
  Data[[i]] <- matrix(0, nTime, nVar)
  Data[[i]][1, ] <- U[[i]][1, ] 
}
# Lagged effects corresponding to transition matrix Phi
for(i in 1:nPerson){
  for (t in 2:nTime){ 
   Data[[i]][t, ] <- Data[[i]][t-1,] %*% Phi + U[[i]][t, ]
   } 
}




# Add differences in within-person stationary means for between-subject networks
# from multivariate distribution with variance-covariance matrix Omega
# implying the following conditional dependencies:
# PA ~ Depr (negative)
# PA ~ Fatigue (negative)
# Stress ~ Depr (positive)
# Stress ~ Fatigue (positive)

mu <- mvrnorm(nPerson, mu=rep(0, nVar), Sigma=Omega)
mu <- 3*mu
for(i in 1:nPerson){
  Data[[i]] <- Data[[i]] + matrix(mu[i, ], nrow=nTime, ncol=nVar, byrow=T)
}


# Unlist and add id variable
Data <- do.call(rbind, Data) |> as.data.frame()
names(Data) <- varnames
Data$id <- rep(1:nPerson, each=nTime)


# Fit N = 1 unregularized VAR ---------------------------------

# Standardize to make coefficients comparable (this is done in graphicalVAR by default)
p1_std <- apply(p1, 2, scale)

# Fit unregularized VAR
# Initialize null matrices
unreg_PDC <- matrix(0, nVar, nVar)
unreg_resid <- matrix(0, nrow(p1_std)-1, nVar)
# GLM loop to predict each var at time=t-1, from itself and every other var at t
for(i in 1:nVar){
  res <- glm(p1_std[-1,varnames[i]] ~ p1_std[-nrow(p1_std),varnames],
             family="gaussian")
  unreg_PDC[, i] <- res$coefficients[-1]
  unreg_resid[, i] <- res$residuals
}

# Since data were standardized, coefficients represent partial correlations
unreg_PDC |> qgraph(layout="circle", theme="colorblind", labels=varnames)
unreg_PCC <- cov(unreg_resid) |> cor2pcor() |> qgraph(layout="circle", theme="colorblind", labels=varnames)





# Fit individual and pooled graphical VAR ---------------------------------


# Individual GVAR networks
individualGVAR <- list()
for(i in 1:nPerson){
  individualGVAR[[i]] <- graphicalVAR(Data[Data$id==i,], 
                                      nLambda=10, vars=varnames)
}

# Temporal and contemporaneous effects of participant 1
qgraph(individualGVAR[[1]]$PDC, layout="circle", theme="colorblind", labels=varnames)
qgraph(individualGVAR[[1]]$PCC, layout="circle", theme="colorblind", labels=varnames)


# Pooled graphical LASSO
wp_cen <- list()
for(i in 1:nPerson){
  wp_cen[[i]] <- apply(Data[Data$id==i, -(nVar+1)], 2, scale)
}
wp_cen <- do.call(rbind, wp_cen) |> as.data.frame()
names(wp_cen) <- varnames
wp_cen$id <- rep(1:nPerson, each=nTime)

pooledGVAR <- graphicalVAR(wp_cen[, -(nVar+1)], nLambda=10)
pooledGVAR$PDC |> qgraph(layout="circle", theme="colorblind", labels=varnames)
pooledGVAR$PCC |> qgraph(layout="circle", theme="colorblind", labels=varnames)




# Fit N > 1 multilevel VAR ---------------------------------------------------------------


res <- mlVAR(Data, vars=varnames,
             idvar = "id", 
             scale=T, #grand-mean centered and scaled
             contemporaneous = "correlated",
             temporal="fixed")

# Average population temporal effects
plot(res, "temporal", layout="circle")
# Average population contemporaneous effects
plot(res, "contemporaneous", layout="circle")
# Between-subject effects
plot(res, "between", layout="circle")








