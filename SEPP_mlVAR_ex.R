##=============================================##
##=============================================##
##                                             ##
## Supplementary code for manuscript titled:   ##
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
## 4/9/2024                                    ##
##============================================ ##
##=============================================##


library(mlVAR)
library(graphicalVAR)
library(qgraph)
library(dplyr)
library(MASS)
library(corpcor)
library(ggplot2)


# Specify true models (synthetic networks) to simulate data from-----------------------------------
# Simulate n > 1 multivariate timeseries data where there are between-person differences in stationary means,
# but not in temporal and contemporaneous effects

nVar <- 4 
varnames <- c("Depr", "Stress", "Fatigue", "PA")
nTime <- 200
nPerson <- 30

# Specify lagged relationships (partial directed correlations)
# Rows = node of origin, columns = destination
# Diagonals = autocorrelations
# Equidistant time intervals
Phi <- matrix(c(0.3, 0.0, 0.0, 0.0,
                0.6, 0.0, 0.0, 0.0,
                0.6, 0.3, 0.3,-0.6,
               -0.6, 0.0,-0.6, 0.0), nrow=4,byrow=T)

truePDCnetwork <- qgraph(Phi, layout="circle",
                         theme="colorblind", labels=varnames,
                         title="True temporal network")

# Specify a contemporaneous network (partial contemporaneous correlations)
# i,j = 0 implies conditional independence between i & j
truePCC <- matrix(c(0.0, 0.0, 0.6, 0.0,
                    0.0, 0.0, 0.3, 0.4,
                    0.6, 0.3, 0.0,-0.4,
                    0.0, 0.4,-0.4, 0.0), nrow=4, byrow=T)

# Recall the VAR model:
# y[t] = Phi * y[t-1] + e
# e ~ N(0, Sigma)
# The variance-covariance matrix Sigma is encodes the PCC
diag(truePCC) <- 1
# Obtain Sigma from the true partial correlation matrix by flipping signs of off-diagonal elements, taking the inverse, and standardizing
Sigma <- pcor2cor(truePCC) 

Sigma |> qgraph(graph="pcor", layout="circle", theme="colorblind",
                labels=varnames, title="True contemporaneous network")

# Check if positive-definite (eigenvalues > 0)
eigen(Sigma)$values > 0

# Specify a variance-covariance matrix for between-person differences
trueBetween <- matrix(c(1.0,  0.6,  0.0, -0.4,
                        0.6,  1.0,  0.5,  0.0,
                        0.0,  0.5,  1.0, -0.4,
                       -0.4,  0.0, -0.4,  1.0), nrow=4, byrow=T)

Omega <- pcor2cor(trueBetween) # Omega encodes the between-subjects network

Omega |> qgraph(graph="pcor", layout="circle", theme="colorblind",
                labels=varnames, title="True between-subjects network")

# Simulate data -----------------------------------------------------------
# Initialize empty data frames for innovations (dynamic errors)
U <- list()
for(i in 1:nPerson){
  U[[i]] <- matrix(0, nrow=nTime, ncol=nVar)
}

# Simulate innovations from multivariate Gaussian distribution with covariance matrix Sigma (contemporaneous relations)
for(i in 1:nPerson){
  U[[i]] <- mvrnorm(n=nTime, mu=rep(0, nVar), Sigma=Sigma)
}

# Create person specific data frames in a list
Data <- list()
for(i in 1:nPerson){
  Data[[i]] <- matrix(0, nTime, nVar)
  Data[[i]][1, ] <- U[[i]][1, ] # Initialize first row
}

# Lagged effects corresponding to transition matrix Phi
for(i in 1:nPerson){
  for (t in 2:nTime){ 
   Data[[i]][t, ] <- Data[[i]][t-1,] %*% Phi + U[[i]][t, ]
   } 
}

# Add differences in within-person stationary means to produce a between-subject networks
# from multivariate Gaussian distribution with variance-covariance matrix Omega
mu <- mvrnorm(nPerson, mu=rep(0, nVar), Sigma=Omega)
mu <- 3*mu # Since small sample in this working example (n=30), increase signal
for(i in 1:nPerson){
  Data[[i]] <- Data[[i]] + matrix(mu[i, ], nrow=nTime, ncol=nVar, byrow=T)
}

# Unlist and add id variable
Data <- do.call(rbind, Data) |> as.data.frame()
names(Data) <- varnames
Data$id <- rep(1:nPerson, each=nTime)




# Fit n = 1 unregularized VAR on random individual ---------------------------------
n1 <- Data[Data$id==sample(1:nPerson, 1), ]
id <- n1$id[1]
# Standardize to make coefficients comparable (this is done in graphicalVAR by default)
n1_std <- apply(n1[, -(nVar+1)], 2, scale)

# Fit unregularized VAR
unreg_PDC <- matrix(0, nVar, nVar)
unreg_resid <- matrix(0, nrow(n1_std)-1, nVar)
# GLM loop to predict each variable at time=t-1, from itself and every other var at t
for(i in 1:nVar){
  res <- glm(n1_std[-1,varnames[i]] ~ n1_std[-nrow(n1_std),varnames],
             family="gaussian")
  unreg_PDC[, i] <- res$coefficients[-1]
  unreg_resid[, i] <- res$residuals
}

# Since data were standardized, coefficients represent partial correlations
unreg_PDC |> qgraph(layout="circle", theme="colorblind", labels=varnames,
                    title=paste("Unregularized VAR for subject ", id))

cov(unreg_resid) |> cor2pcor() |> qgraph(layout="circle", theme="colorblind", labels=varnames,
                                         title=paste("Unregularized contemporaneous network for subject", id))


# True edges are obtained, but there are also many spurious edges. 
# We want to obtain a sparse network where edges that are likely to be spurious are shrunk to 0
# Thus we apply lasso regularization, with tuning parameters (lambda) controlling sparsity in Beta (matrix encoding lagged effects) and Kappa (matrix encoding contemporaneous effects)
# Tuning parameters are selected using a grid search by minimizing the Extended Bayesian Information Criterion
# The EBIC hyperparameter (gamma) is specified to prefer discovery (gamma=0) or parsimony (gamma=0.5) 

# Graphical VAR with LASSO regularization is readily implemented in the graphicalVAR package

# Time intervals are assumed equally spaced
# Yet in typical EMA structure, there are no assessments at night time
# Thus, specify beepvar and dayvar so the first beep of the day is not 
# regressed on the last beep of the previous day
n1 <- as.data.frame(n1)
n1$beep <- rep(1:5, length.out=nrow(n1))
n1$day <- rep(1:(nTime/5), each=5)


n1_GVAR <- graphicalVAR(n1, nLambda=30, vars=varnames,
                        beepvar="beep", dayvar="day")

n1_PDC <- n1_GVAR$PDC |> qgraph(layout="circle", theme="colorblind", labels=varnames,
                                title=paste("Regularized VAR for subject ", id))

n1_PCC <- n1_GVAR$PCC |> qgraph(layout="circle", theme="colorblind", labels=varnames,
                                title=paste("Regularized contemporaneous network for subject ", id))

# Centrality metrics for individual networks
centralityTable(n1_GVAR$PDC) |>
  filter(measure %in% c("InStrength", "OutStrength")) |>
  ggplot(aes(x = node, y = value, group = measure)) +
  geom_line(aes(linetype = measure), linewidth = 1) +
  labs(x = '', y = 'z-score') +
  scale_linetype_discrete(name = "Centrality") +
  coord_flip() +
  theme_bw()

# Visually compare to a randomly generated VAR
rand_VAR <-randomGVARmodel(Nvar = nVar, probKappaEdge=0.5, probKappaPositive=0.5,
                           probBetaEdge=0.5, probBetaPositive=0.5)


par(mfrow=c(2,2))
plot(n1_PDC)
plot(n1_PCC)
rand_VAR$PDC |> qgraph(layout="circle", theme="colorblind", labels=varnames,
                   title="Randomly generated PDC")
rand_VAR$PCC |> qgraph(layout="circle", theme="colorblind", labels=varnames,
                   title="Randomly generated PCC")

par(mfrow=c(1,1))


# Fit N > 1 multilevel VAR ---------------------------------------------------------------
# Two-step multilevel VAR as described in Epskamp, S., Waldorp, L. J., Mõttus, R., & Borsboom, D. (2018). The Gaussian graphical model in cross-sectional and time-series data. Multivariate behavioral research, 53(4), 453-480.
# And implemented in the mlVAR package
res <- mlVAR(Data, vars=varnames,
             idvar = "id", 
             estimator = "lmer", # two-step multilevel VAR
             scale=T, # grand-mean centered and scaled
             contemporaneous = "correlated",
             temporal="correlated") # correlated random-effects, works well with < 8 variables

summary(res)
res$output$temporal

# Average (fixed) population temporal effects
temporal <- plot(res, "temporal", theme="colorblind", layout="circle",
                 title=paste("Temporal network of n =",nPerson,"subjects"))
# Average population contemporaneous effects
contemporaneous<- plot(res, "contemporaneous", theme="colorblind", layout="circle",
                       title=paste("Contemporaneous network of n =",nPerson,"subjects"))
# Between-subject effects (differences between stationary means), power to estimate depends on how many nPersons
between <- plot(res, "between", theme="colorblind", layout="circle",
                title=paste("Between-subjects network of n =",nPerson,"subjects"))

# Centrality
centralityPlot(temporal, include=c("InStrength", "OutStrength"), scale="z-scores")
centralityPlot(contemporaneous, include=c("ExpectedInfluence"), scale="z-scores")

# Individual subjects
# In multilevel estimation, individual-level estimates are shrunk towards fixed-effects
rand_subject <- sample(1:nPerson, 1)

plot(res, "temporal", theme="colorblind", layout="circle",
     subject=rand_subject,
     title=paste("Individual estimation of subject",rand_subject))










