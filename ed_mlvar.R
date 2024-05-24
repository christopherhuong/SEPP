

library(dplyr)
library(mlVAR)
library(imputeTS)
library(qgraph)

data_full <- readRDS("G:/My Drive/R Code Examples/kolar-2024-exercise moderates eating pathology networks/01-data/data_levinson.RData")

data_full <- data_full[, c("unique_id", "dh5", "dh60", "dh58", "dh56", "dh29", "dh52",
                               "dh54", "dh15", "dh51", "dh61", "dh53", 
                               "dh55", "dh59", "day", "time_point")] %>%
  mutate(across(c("dh60", "dh58", "dh56", "dh29", "dh52", "dh54", "dh15", "dh51", "dh61", "dh53", "dh55", 
                  "dh59", "time_point"), as.numeric))# change factor to numeric



data_full$dh5 <- data_full$dh5 %>%
  recode(
    Yes = '1',
    No = '0'
  ) %>%
  as.character() %>%
  as.numeric()


data_full_imp <- na_kalman(data_full, model = "StructTS", smooth = TRUE)


colnames(data_full_imp) <- c("ID", "exerc", "binge", "bodycheck", "compens", "deserve",
                         "diet", "feargain", "guilty", "overeat", "restrict", 
                         "thinner", "vomit", "weigh", "day", "time_point")


# All ED symptom vars
# vars1 <-  c("exerc", "binge", "bodycheck", "compens", "deserve",
#             "diet", "feargain", "guilty", "overeat", "restrict", 
#             "thinner", "vomit", "weigh")
# 
# mod1 <- mlVAR(data_full_imp, vars=vars1, idvar="ID", lags=1, dayvar="day", beepvar="time_point",
#              estimator="lmer", temporal="orthogonal", contemporaneous="orthogonal")
# 
# 
# plot(mod1, "temporal", title="mlvar")

# For visual demonstration purposes, retain variables that have direct connections to exerc. 
# Note that choices of node in/exclusion affect edge estimates,
# and thus careful consideration should be taken when deciding which nodes to include in estimating a symptom network

vars2 <- c("exerc", "guilty", "feargain", "bodycheck", "weigh", "restrict")

mod2 <- mlVAR(data_full_imp, vars=vars2, idvar="ID", lags=1, dayvar="day", beepvar="time_point",
              estimator="lmer", temporal="correlated", contemporaneous="correlated")

plot(mod2, "temporal", title="temporal network", layout="circle")
plot(mod2, "contemporaneous", title="contemporaneous network", layout="circle")
# Retrieve weighted adjacency matrix from estimated beta matrix to plot directed networks.
# Note that the beta matrix must be transposed when plotting directed networks, as per convention
# (rows = node of origin, columns = destination )
Beta <- mod2$results$Beta$mean |> as.data.frame() |> as.matrix() |> t()
# Set non-significant edges to 0
p_mat_Beta <- mod2$results$Beta$P |> as.data.frame() |> as.matrix() |> t()

for(i in 1:length(vars2)){
  for(j in 1:length(vars2)){
    Beta[i,j] <- ifelse(p_mat_Beta[i,j] < 0.05, Beta[i,j], 0)
  }
}

# Remove self-loops for clearer plotting
diag(Beta) <- 0


# Retreive contemporaneous network Kappa
Kappa <- plot(mod2, "contemporaneous", nonsig="hide", DoNotPlot=T)
Kappa <- Kappa$Arguments$input

# Inverting Kappa and flipping the signs of the off-diagonal elements obtains the variance-covariance matrix of the residuals, Theta 
Theta <- solve(Kappa)
Theta <- -Theta
diag(Theta) <- 1


# use Beta and Theta from estimated network as "true" data generating parameters to simulate data from
nPerson <- 100
nTime <- 200
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






















