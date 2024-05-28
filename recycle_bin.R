
# use Beta and Theta from estimated network as "true" data generating parameters to simulate data from
nPerson <- 50
nTime <- 100
nVar <- length(vars2)
# Initialize empty data frames for innovations (dynamic errors)
U <- list()
for(i in 1:nPerson){
  U[[i]] <- matrix(0, nrow=nTime, ncol=nVar)
}

# Simulate innovations from multivariate Gaussian distribution with covariance matrix Theta (encodes contemporaneous relations)
for(i in 1:nPerson){
  U[[i]] <- MASS::mvrnorm(n=nTime, mu=rep(0, nVar), Sigma=Theta)
}

# Create person specific data frames in a list
Data <- list()
for(i in 1:nPerson){
  Data[[i]] <- matrix(0, nTime, nVar)
  Data[[i]][1, ] <- U[[i]][1, ] # Initialize first row
}

# Lagged effects corresponding to Beta matrix
for(i in 1:nPerson){
  for (t in 2:nTime){ 
    Data[[i]][t, ] <- Data[[i]][t-1,] %*% Beta + U[[i]][t, ]
  } 
}

# Unlist and add id variable
Data <- do.call(rbind, Data) |> as.data.frame()
names(Data) <- vars2
Data$id <- rep(1:nPerson, each=nTime)


mod3 <- mlVAR(Data, vars=vars2,
              idvar = "id", 
              estimator = "lmer", # two-step multilevel VAR
              scale=T, # grand-mean centered and scaled
              contemporaneous = "correlated",
              temporal="correlated") # correlated random-effects, works well with < 8 variables


plot(mod3, "temporal", layout="circle")
plot(mod3, "contemporaneous", layout="circle")



