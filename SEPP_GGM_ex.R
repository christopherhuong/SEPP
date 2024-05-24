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



library(huge)
library(qgraph)
library(bootnet)

setwd("G:/My Drive/INCH Projects/[2023] SEPP Special Issue/datacode")


# stretch <- read.csv('StretchingGroupData.csv')
# cycle <- read.csv('CyclingGroupData.csv')
#
#
# stretch <- huge.npn(stretch, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE)
# stretch <- as.data.frame(apply(stretch, use='pairwise.complete.obs', 2, as.numeric))
#
# cycle <- huge.npn(cycle, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE)
# cycle <- as.data.frame(apply(cycle, use='pairwise.complete.obs', 2, as.numeric))
#
#
names = c("Perserv", "Neg", "Critic", "Brood", "Replay")
#
#
# stretch_glasso <- estimateNetwork(stretch, default="EBICglasso", tuning=.5)
# cycle_glasso <- estimateNetwork(cycle, default="EBICglasso", tuning=.5)
#
#
# stretch_GGM <- stretch_glasso$graph
# save(stretch_GGM, file="stretch_GGM.RData")
# cycle_GGM <- cycle_glasso$graph
# save(cycle_GGM, file="cycle_GGM.RData")

# Load  GGMs estimated on data from Bernstein, E. E., Heeren, A., & McNally, R. J. (2020). A network approach to understanding the emotion regulation benefits of aerobic exercise. Cognitive Therapy and Research, 44, 52-60.
load("stretch_GGM.RData") # stretch control
load("cycle_GGM.RData") # exercise intervention

# Plot 
stretch_GGM |> qgraph(title = "Stretching", layout="circle", edge.labels=FALSE,
                       labels=names, theme = "colorblind")


cycle_GGM |> qgraph(title = "Cycling", layout="circle", edge.labels=FALSE,
                     labels=names, theme = "colorblind")


# Simulate data from estimated networks using GGMsim() from https://github.com/RiaHoekstra/simulation_functions
source("https://raw.githubusercontent.com/RiaHoekstra/simulation_functions/main/GGMsim.R")
set.seed(123)
stretch_sim <- GGMsim(n=113, omega=stretch_GGM)
cycle_sim <- GGMsim(n=113, omega=cycle_GGM)

# # GGMsim() decomposed:
# omega <- stretch_GGM
# # Obtain variance-covariance matrix Sigma
# sigma <- (diag(ncol(omega)) - omega) |> # flip signs of off-diagonals, and add 1's to diagonal
#   solve() |> # take inverse
#   cov2cor() # standardize
# stretch_sim <- MASS::mvrnorm(n=113, mu=rep(0, ncol(omega)), Sigma=sigma)


# Fit networks ------------------------------------------------------------

stretch_glasso <- estimateNetwork(stretch_sim, default="EBICglasso", tuning=.5)
cycle_glasso <- estimateNetwork(cycle_sim, default="EBICglasso", tuning=.5)


# Plot --------------------------------------------------------------------
L <- averageLayout(stretch_glasso, cycle_glasso)

stretch_plot <- qgraph(stretch_glasso$graph, title = "Stretching",
                       layout=L, labels=names, theme = "colorblind")


cycle_plot <- qgraph(cycle_glasso$graph, title = "Cycling",
                     layout=L, labels=names, theme = "colorblind")



all_cen <- list("Stretching"=stretch_glasso, "Cycling"=cycle_glasso)
cen_plot <- centralityPlot(all_cen, include = "Strength", labels = names, scale = "z-scores")


pdf(file="Fig1.pdf", width = 6, height = 12)
par(mfrow = c(2,1))
plot(stretch_plot)
plot(cycle_plot)
plot(cen_plot)
dev.off()


# Reliability check -------------------------------------------------------

# Procedure outlined in Epskamp, S., & Fried, E. I. (2018). A tutorial on regularized partial correlation networks. Psychological methods, 23(4), 617.
bootnet()





# Network Comparison Test -------------------------------------------------
library(BGGM)

t <- ggm_compare_ppc(stretch_sim, cycle_sim,
                     test="global",
                     iter=500)

t$ppp_jsd
p <- plot(t)
p$plot_jsd


# Moderated network -------------------------------------------------------










