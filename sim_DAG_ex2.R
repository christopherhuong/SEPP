
library(qgraph)
library(mgm)
library(bootnet)
library(dplyr)


vars <- c("sed", "depr", "fatigue", "conc")
# Sample size of each group
N <- 250
# 1 = PA, 0 = Control
PA.1 <- rep(1, N)
PA.0 <- rep(0, N)
PA <- c(PA.1, PA.0)
# Depressed mood = -0.8*PA + e~N(0,1)
depr <- rnorm(N, -0.8*PA, 1)
# When PA=1: Fatigue = -0.8*PA  + e~N(0,1)
# When PA=0: Fatigue = -0.8*PA + 0.4*depr + e~N(0,1)
fatigue <- -0.8*PA 



# Simulate sedentary=1 & active=0
sed <- rbinom(N, 1, p=0.5) 
# depr = sed*0.8 + e~N(0,1)
depr <- rnorm(N, sed*0.8, 1)
# fatigue = sed*0.8 + depr*sed*0.4 + e~N(0,1)
# depr -> fatigue only when sed=1
fatigue <- rnorm(N, sed*0.8 + depr*sed*0.4, 1)
# when sed=0: conc = fatigue*0.4 + e~N(0,1)
# when sed=1: conc = fatigue*0.8 + e~N(0,1)
conc <- ifelse(sed==0, rnorm(1, fatigue*0.4, 1), rnorm(1, fatigue*0.8, 1))


d <- data.frame(sed=sed, depr=depr, fatigue=fatigue, conc=conc)
names(d) <- vars
d_0 <- d[d$sed==0, 2:4]
d_1 <- d[d$sed==1, 2:4]

d <- as.matrix(d)


mgm_all <- mgm(data=d, type=c("c", "g", "g", "g"), level = c(2, 1, 1, 1),
               lambdaSel = "EBIC", lambdaGam = 0.5)

qgraph(mgm_all$pairwise$wadj, edge.color=mgm_all$pairwise$edgecolor,
       layout="spring", labels=vars, edge.labels=T, title="All")



glasso_0 <- estimateNetwork(d_0, default="EBICglasso", tuning=0.5)

glasso_1 <- estimateNetwork(d_1, default="EBICglasso", tuning=0.5)

# Plot all graphs
layout(t(1:3))
qgraph(mgm_all$pairwise$wadj, edge.color=mgm_all$pairwise$edgecolor,
       layout="spring", labels=vars, edge.labels=T, title="All")
qgraph(glasso_0$graph, layout="circle", labels=vars[2:4], title="active", edge.labels=T)
qgraph(glasso_1$graph, layout="circle", labels=vars[2:4], title="sedentary", edge.labels=T)


# Specify PA as moderator
mnm <- mgm(data=d,
           type=c(rep("g", 4),"c"),
           level = c(rep(1, 4),2),
           moderators = 5,
           lambdaSel = "CV",
           ruleReg = "AND")


l_mgm_cond <- list()
for(g in 1:2) l_mgm_cond[[g]] <- condition(mnm, values=list("5" = g))


v_max <- rep(NA, 2)
for(g in 1:2) v_max[g] <- max(l_mgm_cond[[g]]$pairwise$wadj)




par(mfrow=c(1,2))
for(g in 1:2) {
  qgraph(input = l_mgm_cond[[g]]$pairwise$wadj, 
         edge.color = l_mgm_cond[[g]]$pairwise$edgecolor, 
         layout = "circle", labels=vars,
         maximum = max(v_max), vsize=16, esize=23, 
         edge.labels  = TRUE, edge.label.cex = 3)
  mtext(text = paste0("Group ", g), line=2.5)
}




par(mfrow=c(1,2))
for(g in 1:2) {
  qgraph(input = l_mgm_cond[[g]]$pairwise$wadj, 
         edge.color = l_mgm_cond[[g]]$pairwise$edgecolor, 
         layout = "circle", labels=vars,
         maximum = max(v_max), vsize=16, esize=23, 
         edge.labels  = TRUE, edge.label.cex = 3)
  mtext(text = paste0("Group ", g), line=2.5)
}




















