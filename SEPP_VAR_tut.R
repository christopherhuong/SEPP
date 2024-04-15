#################################################
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
## 3/26/2024                                   ##
##                                             ##
#################################################


library(MASS) 
library(graphicalVAR)
library(qgraph)
library(corpcor)
library(ggplot2)
library(dplyr)


# Specify a data generating model
# Simulate an individuals data
# then fit the an unregularized and regularized graphical VAR model
nVar <- 4 
nTime <- 200 
varnames <- c("Depr", "Stress", "Fatigue", "PA")

# Specify a matrix of hypothetical lagged relations
# Rows = node of origin, columns = destination
# Diagonals = autocorrelations
Phi <- matrix(c(0.3, 0.0, 0.0, 0.0,
                0.6, 0.0, 0.0, 0.0,
                0.6, 0.3, 0.3,-0.6,
               -0.6, 0.0,-0.6, 0.0), nrow=4,byrow=T)

truePDC <- qgraph(Phi, layout="circle",
                  theme="colorblind", labels=varnames)


eigen(Phi)$values

# Specify a contemporaneous variance-covariance matrix 
Sigma <- matrix(c(1.0,  0.1,  0.5, -0.1,
                  0.1,  1.0,  0.2,  0.3,
                  0.5,  0.2,  1.0, -0.2,
                 -0.1,  0.3, -0.2,  1.0), nrow=4,byrow=T)

# True partial contemporaneous correlations (PCC)
# Sigma is inverted (to obtain the precision matrix) and standardized 
# to produce the partial correlation matrix
truePCC <- qgraph(Sigma, graph="pcor",
                  layout="circle", theme="colorblind", labels=varnames)

# Check if positive-definite (eigenvalues > 0)
eigen(Sigma)$values > 0

# Simulate innovations (dynamic errors) from multivariate Gaussian distribution with covariance matrix Sigma
set.seed(999)
U <- mvrnorm(n = nTime, mu = rep(0, nVar), Sigma = Sigma) 

p1 <- matrix(0, nrow=nTime, ncol=nVar) 

# Produce data set with 1-lagged relationships according to Phi
p1[1,] <- U[1,] 
for (t in 2: nTime){ 
  p1[t,] <- p1[t-1,] %*% Phi + U[t,] 
} 

colnames(p1) <- varnames
# Could add more noise to all observations from rnorm()

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
p1 <- as.data.frame(p1)
p1$id <- 1
p1$beep <- rep(1:5, length.out=nrow(p1))
p1$day <- rep(1:40, each=5)


res <- graphicalVAR(p1, nLambda=20,
                    vars=varnames,
                    beepvar="beep",
                    dayvar="day")

res$PDC |> qgraph(layout="circle", theme="colorblind", labels=varnames, fade=F)
res$PCC |> qgraph(layout="circle", theme="colorblind", labels=varnames, fade=F)

# Directed network included a spurious autoregressive effect of PA,
# but otherwise much improved recovery of true network structure

# Centrality
df <- centralityTable(res$PDC) |>
  filter(measure %in% c("InStrength", "OutStrength"))

df |>
  ggplot(aes(x = node, y = value, group = measure)) +
  geom_line(aes(linetype = measure), linewidth = 1) +
  labs(x = '', y = 'z-score') +
  scale_linetype_discrete(name = "Centrality") +
  coord_flip() +
  theme_bw()








# Fit individual and pooled graphical VAR on entire data set ---------------------------------

# Recycle bin from mlVAR_tut
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

pooledGVAR <- graphicalVAR(wp_cen[, -(nVar+1)], nLambda=20)
pooledGVAR$PDC |> qgraph(layout="circle", theme="colorblind", labels=varnames,
                         title="Pooled GVAR PDC")
pooledGVAR$PCC |> qgraph(layout="circle", theme="colorblind", labels=varnames,
                         title="Pooled GVAR PCC")


