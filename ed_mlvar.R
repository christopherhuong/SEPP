

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
vars1 <-  c("exerc", "binge", "bodycheck", "compens", "deserve",
            "diet", "feargain", "guilty", "overeat", "restrict", 
            "thinner", "vomit", "weigh")

mod1 <- mlVAR(data_full_imp, vars=vars1, idvar="ID", lags=1, dayvar="day", beepvar="time_point",
             estimator="lmer", temporal="orthogonal", contemporaneous="orthogonal")


plot(mod1, "temporal", title="mlvar")

# For visual demonstration purposes, retain variables that have direct connections to exerc. 
# Note that choices of node in/exclusion affect edge estimates,
# and thus careful consideration should be taken when deciding which nodes to include in estimating a symptom network

vars2 <- c("exerc", "guilty", "feargain", "bodycheck", "weigh", "restrict")

mod2 <- mlVAR(data_full_imp, vars=vars2, idvar="ID", lags=1, dayvar="day", beepvar="time_point",
              estimator="lmer", temporal="orthogonal", contemporaneous="orthogonal")

# Retrieve weighted adjacency matrix from estimated beta matrix.
# Note that the beta matrix must be transposed when plotting directed networks, as per convention
# (rows = node of origin, columns = destination )
wadjmat <- mod2$results$Beta$mean |> as.data.frame() |> as.matrix() |> t()
# Set non-significant edges to 0
p_mat <- mod2$results$Beta$P |> as.data.frame() |> as.matrix() |> t()

for(i in 1:length(vars2)){
  for(j in 1:length(vars2)){
    wadjmat[i,j] <- ifelse(p_mat[i,j] < 0.05, wadjmat[i,j], 0)
  }
  
}

# remove self-loops for clear plotting
diag(wadjmat) <- 0
qgraph(wadjmat, labels=vars2, title="qgraph", theme="colorblind", layout="circle")


# Individual subjects
# In multilevel estimation, individual-level estimates are shrunk towards fixed-effects
# Pull out 4 random subjects
set.seed(123)
random_subjects <- sample(1:length(unique(data_full_imp$ID)), 4)

par(mfrow=c(2,2))
for(i in 1:length(random_subjects)){
  plot(mod2, "temporal", theme="colorblind", layout="circle",
       subject=random_subjects[i],
       title=paste("Individual estimation of subject",random_subjects[i]))
}






















