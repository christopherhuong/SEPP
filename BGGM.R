

library(bootnet)
library(qgraph)
library(MASS)
library(BGGM)

# Implementing two Bayesian approaches to compare Gaussian graphical models (GGMs)
# From:
# Williams, D. R., Rast, P., Pericchi, L. R., & Mulder, J. (2020). Comparing Gaussian graphical models with the posterior predictive distribution and Bayesian model selection. Psychological methods, 25(5), 653.


# Simulate 2 random GGMs
nVar=5
vars <- c("depr_mood", "fatigue", "conc_issues", "anhedonia", "stress")
nSample = 1e3
set.seed(123)
PA <- rbinom(n=nSample, size=1, p=0.5)
GGM_0 <- genGGM(nVar, p=0.5, propPositive=0.8)
GGM_1 <- genGGM(nVar, p=0.5, propPositive=0.8)


# To simulate multivariate Gaussian distributed data, get variance-covariance matrix Sigma
# by setting diagonal elements=1, reversing the signs of the off-diagonals, taking the inverse, and standardizing 
Sigma0 <- solve(diag(ncol(GGM_0)) - GGM_0) |> cov2cor()
Sigma1 <- solve(diag(ncol(GGM_1)) - GGM_1) |> cov2cor()


# Check if covariance matrix is positive semi definite
eigen(Sigma0)$values >= 0
eigen(Sigma0)$values >= 0


# Simulate data for control and treatment groups
data_0 <- mvrnorm(n=nSample, mu=rep(0, nVar), Sigma=Sigma0) |>
  cbind(rep(0, nSample))
# Direct effect of PA on depr_mood and fatigue
data_1 <- mvrnorm(n=nSample, mu=c(0.3, 0.3, rep(0, nVar-2)), Sigma=Sigma1) |>
  cbind(rep(1, nSample))

data <- rbind(data_0, data_1) |> as.data.frame()


# Estimate GGM on all data, and GGMs for each condition using lasso regularization
# with lambda selected by minimizing the EBIC, with hyperparameter gamma set to 0.5
ggm_fit_all <- estimateNetwork(data=data, default="EBICglasso", tuning=0.5)
ggm_fit_0 <- estimateNetwork(data=data[data$V6==0, 1:nVar], default="EBICglasso", tuning=0.5)
ggm_fit_1 <- estimateNetwork(data=data[data$V6==1, 1:nVar], default="EBICglasso", tuning=0.5)
ggm_list <- list(ggm_all=ggm_fit_all, ggm_0=ggm_fit_0, ggm_1=ggm_fit_1)

# Plot
par(mfrow=c(3,1))
for(v in 1:3){
  qgraph(ggm_list[[v]]$graph, title=names(ggm_list[v]))
}

# To compare 















