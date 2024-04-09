
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


library(graphicalVAR)
library(mlVAR)
library(qgraph)
library(dplyr)
library(ggplot2)


# Code are heavily adapted from the graphicalVAR package

# Epskamp S (2021). _graphicalVAR: Graphical VAR for
# Experience Sampling Data_. R package version 0.3,
# <https://CRAN.R-project.org/package=graphicalVAR>.

# https://github.com/SachaEpskamp/graphicalVAR


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
trueSigma <- matrix(c(1.91,  0.69, 1.01,  0.87,
                      0.69,  1.39, 0.14, -0.14,
                      1.01,  0.14, 1.91,  1.22,
                      0.87, -0.14, 1.22,  1.91),nrow=4,byrow=T)

# Get PCC (partial contemporaneous correlations) by standardizing 
# the precision matrix (Kappa) and multiplying the off-diagonal elements by -1
trueKappa <- solve(trueSigma)
truePCC <- -trueKappa / (sqrt(diag(trueKappa)) %*% t(sqrt(diag(trueKappa))))
diag(truePCC) <- 1
truePCC
qgraph(truePCC, layout="circle", labels=varnames)

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
trueBeta <- matrix(c(0.10, 0.60, -0.90, -0.81,
                     0.30, 0.10,  0.00,  0.00,
                     0.00, 0.00,  0.20, -0.75,
                     0.00,-0.40, -0.70, -0.20), nrow=4,byrow=T)

# PDC = partial directed correlations
# conditional dependencies between time points
truePDC <- t(trueBeta / sqrt( diag(trueSigma) %o% diag(trueKappa) + trueBeta^2))

qgraph(truePDC, layout="circle", labels=varnames)




# Simulate N=1 observations from hypothetical true model ---------------------------------------------------

set.seed(2024)
p1 <- graphicalVARsim(nTime = 100,
                      beta = trueBeta,
                      kappa = trueKappa)


p1 <- as.data.frame(p1)
colnames(p1) <- varnames


# Fit the unregularized VAR model ------------------------------------------------------
p1_unreg <- p1

# Center and scale
p1_unreg[, varnames] <- scale(p1_unreg[, varnames], center=T,scale=T)

# Compute current and lagged data
p1_current <- p1_unreg[-1,, drop=F] |> as.matrix()
p1_lagged <- cbind(1, p1_unreg[-nrow(p1_unreg),,drop=F]) |> as.matrix()


SigmaHat <- cov(p1_unreg[, varnames])
L1Hat <- cov(p1_lagged[, varnames], p1_current[, varnames])
beta <- t(solve(SigmaHat) %*% L1Hat)
kappa <- solve(SigmaHat - beta %*% SigmaHat %*% t(beta))

unreg_PCC <- -kappa / (sqrt(diag(kappa)) %*% t(sqrt(diag(kappa))))
diag(unreg_PCC) <- 0
unreg_PDC <- t(beta / sqrt( diag(solve(kappa)) %o% diag(kappa) + beta^2))

qgraph(unreg_PCC, layout="circle", labels=varnames, fade=F)
qgraph(unreg_PDC, layout="circle", labels=varnames, fade=F)




# Many spurious edges. We want to obtain a sparse network where edges that are likely to be spurious are shrunk to 0
# Thus we apply lasso regularization, with tuning parameters (lambda) controlling sparsity in Beta and Kappa
# Tuning parameters are selected by minimizing the Extended Bayesian Information Criterion
# The EBIC hyperparameter (gamma) is specified to prefer discovery (gamma=0)
# or parsimony (gamma=0.5) 


# Graphical VAR with LASSO regularization is readily implemented in the graphicalVAR package

# Optional to have structure mimic EMA, collected 5x per day for 20 days
# Time intervals are assumed equally spaced
# Specify beepvar and dayvar so the first beep of the day is not 
# regressed on the last beep of the previous day
# p1$beep <- rep(1:5, 20)
# p1$day <- rep(1:20, each=5)
# str(p1)

set.seed(2024)
gvar1 <- graphicalVAR(p1, vars=varnames,
                      nLambda=20, gamma=0.5, 
                      scale=T,lags=1)

gvar_pcc <- qgraph(gvar1$PCC, layout="circle", labels=varnames)
gvar_pdc <- qgraph(gvar1$PDC, layout="circle", labels=varnames)



# The true edge between PA-Anhed in the contemporaneous model was not obtained
# The small autoregressive effects were also not obtained
# Otherwise, we can see better recovery of the true PCC and PDC models



# Centrality --------------------------------------------------------------


# In-strength and out-strength: 
# sum of magnitudes of edges entering and exiting a node 

df <- centralityTable(gvar_pdc) %>%
  filter(measure %in% c("InStrength", "OutStrength"))

df %>%
  ggplot(aes(x = node, y = value, group = measure)) +
  geom_line(aes(linetype = measure), linewidth = 1) +
  labs(x = '', y = 'z-score') +
  scale_linetype_discrete(name = "Centrality") +
  coord_flip() +
  theme_bw()



# Multilevel VAR -----------------------------------------------------------

# Simulate a random model and observations
set.seed(1)
ml_sim <- mlVARsim(nPerson=50,
                   nNode=4,
                   nTime=100,
                   lag=1)

names(ml_sim$Data) <- c(varnames, "ID")
# Mimic EMA data 5x/day for 20 days
ml_sim$Data$beep <- rep(1:5)
ml_sim$Data$day <- rep(1:20, each=5)


str(ml_sim$Data)


# Fit the model
mlvar_mod <- mlVAR(ml_sim$Data,
                   vars=varnames,
                   lags=1,
                   idvar="ID",
                   dayvar = "day",
                   beepvar = "beep",
                   temporal="correlated" #correlated random effects
                   )

print(mlvar_mod)

# Parameter estimates
summary(mlvar_mod) 

# Compare simulated model with fitted model
plot(ml_sim, "temporal", title = "True temporal relationships", layout = "circle")
plot(mlvar_mod, "temporal", title = "Estimated temporal relationships", layout = "circle", labels=varnames)


plot(ml_sim, "contemporaneous", title = "True contemporaneous relationships",
     layout = "circle")
plot(mlvar_mod, "contemporaneous", title = "Estimated contemporaneous relationships",
     layout = "circle")











