

library(huge)
library(qgraph)
library(bootnet)
library(mgm)
library(NetworkComparisonTest)


setwd("G:/My Drive/INCH Projects/[2023] SEPP Special Issue/datacode")





stretch <- read.csv('StretchingGroupData.csv')
cycle <- read.csv('CyclingGroupData.csv')

stretch <- huge.npn(stretch, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE)
stretch <- as.data.frame(apply(stretch, use='pairwise.complete.obs', 2, as.numeric))

cycle <- huge.npn(cycle, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE)
cycle <- as.data.frame(apply(cycle, use='pairwise.complete.obs', 2, as.numeric))



load("stretch_GGM.RData") # stretch control
load("cycle_GGM.RData") # exercise intervention


source("https://raw.githubusercontent.com/RiaHoekstra/simulation_functions/main/GGMsim.R")
set.seed(123)
stretch_sim <- GGMsim(n=113, omega=stretch_GGM)
cycle_sim <- GGMsim(n=113, omega=cycle_GGM)


stretch_sim <- as.data.frame(stretch_sim)
cycle_sim <- as.data.frame(cycle_sim)




stretch_sim$exer <- 2
cycle_sim$exer <- 1

d <- rbind(stretch_sim, cycle_sim)


names(d) = c("Perserv", "Neg", "Critic", "Brood", "Replay", "Exer")



d_mgm <- mgm(data=d,
             type=c(rep("g", 5), "c"),
             level = c(rep(1, 5), 2),
             moderators=6,
             lambdaSel = "EBIC",
             lambdaGam = 0.25,
             ruleReg = "OR")

l_mgm_cond <- list()
for(g in 1:2) l_mgm_cond[[g]] <- condition(d_mgm, values=list("6" = g))


v_max <- rep(NA, 2)
for(g in 1:2) v_max[g] <- max(l_mgm_cond[[g]]$pairwise$wadj)





pdf("mgm.pdf", width = 9, height = 3)
par(mfrow=c(1,2))
for(g in 1:2) {
  qgraph(input = l_mgm_cond[[g]]$pairwise$wadj, 
         edge.color = l_mgm_cond[[g]]$pairwise$edgecolor, 
         layout = "circle", mar=c(2,3,5,3),
         maximum = max(v_max), vsize=16, esize=23, 
         edge.labels  = TRUE, edge.label.cex = 3)
  mtext(text = paste0("Group ", g), line=2.5)
}
dev.off()






# d_mgm <- mgm(data=d, 
#                 type=c(rep("g", 5), "c"),
#                 level = c(rep(1, 5), 2)
#                 )
# 
# qgraph(d_mgm$pairwise$wadj,
#        edge.color=d_mgm$pairwise$edgecolor)
# 
# 
# 
# d$tot <- rowSums(d[, 1:5])
# 
# lm(tot ~ exer, d) |> summary()
# 
# 



# stretch_glasso <- estimateNetwork(stretch, default="EBICglasso", tuning=.5)
# cycle_glasso <- estimateNetwork(cycle, default="EBICglasso", tuning=.5)
# 
# nct1 <- NCT(stretch_glasso, cycle_glasso,
#             test.edges=F,
#             test.centrality=T)
# nct1
# summary(nct1)
# 









