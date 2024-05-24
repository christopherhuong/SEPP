

library(simDAG)
library(bootnet)
library(mgm)
library(qgraph)



# Initialize an empty DAG for PA = 1
dag_2 <- empty_dag()
dag_2 <- dag_2 + node("PA", type="rbernoulli", p=1)
# Define relations between parent and child nodes
dag_2 <- dag_2 +
  # depr ~ 0 + PA*-0.3 + e~N(0,0.2)
  node("depr", type="gaussian", parents=c("PA"),
       betas=c(-0.3), intercept=0, error=0.2) +
  # fatigue ~ 0 + PA*-0.3 + e~N(0,0.2)
  node("fatigue", type="gaussian", parents=c("PA"),
       betas=c(-0.3), intercept=0, error=0.2) +
  # irritable ~ 0 + fatigue*0.3 + e~N(0, 0.2)
  node("irritable", type="gaussian", parents=c("fatigue"),
       betas=c(0.3), intercept=0, error=0.2) +
  # conc_issue ~ 0 + fatigue*0.3 + e~N(0, 0.2)
  node("conc_issue", type="gaussian", parents=c("fatigue"),
       betas=c(0.3), intercept=0, error=0.2)


# Initialize an empty DAG for PA = 0
# Differences: 
# add depr -> fatigue
# doubled effect of fatigue -> irritable
# doubled effect of fatigue -> conc_issue
dag_1 <- empty_dag()
dag_1 <- dag_1 + node("PA", type="rbernoulli", p=0)
# Define relations between parent and child nodes
dag_1 <- dag_1 +
  # depr ~ 0 + PA*-0.3 + e~N(0,0.2)
  node("depr", type="gaussian", parents=c("PA"),
       betas=c(-0.3), intercept=0, error=0.2) +
  # fatigue ~ 0 + PA*-0.3 + depr*0.3 + e~N(0,0.2)
  node("fatigue", type="gaussian", parents=c("PA", "depr"),
       betas=c(-0.3, 0.3), intercept=0, error=0.2) +
  # irritable ~ 0 + fatigue*0.3 + e~N(0, 0.2)
  node("irritable", type="gaussian", parents=c("fatigue"),
       betas=c(0.6), intercept=0, error=0.2) +
  # conc_issue ~ 0 + fatigue*0.3 + e~N(0, 0.2)
  node("conc_issue", type="gaussian", parents=c("fatigue"),
       betas=c(0.6), intercept=0, error=0.2)


# Plot
vars <- colnames(dag2matrix(dag_2))
adjmat_2 <- dag2matrix(dag_2)
adjmat_1 <- dag2matrix(dag_1)[vars, vars]
  
  
par(mfrow=c(2, 1))
qgraph(adjmat_2, layout="circle", title="PA = 1", labels=vars)
qgraph(adjmat_1, layout="circle", title="PA = 0", labels=vars)

# Simulate data from DAGs
set.seed(123)
nSample <- 300
sim_2 <- sim_from_dag(dag=dag_2, n_sim = nSample)
sim_1 <- sim_from_dag(dag=dag_1, n_sim = nSample)
# rbind
sim_d <- rbind(sim_2, sim_1)
# Compute sum score of symptoms
sim_d$depr_tot <- rowSums(sim_d[, 2:5])
# PA condition = 2, control = 1
sim_d$PA <- ifelse(sim_d$PA == T, 2, 1)

cov(sim_d) |> round(2)
cov(sim_d[sim_d$PA==2]) |> round(2)
cov(sim_d[sim_d$PA==1]) |> round(2)

# Linear effect of PA on depression sum-score
lm(depr_tot ~ PA, data=sim_d) |> summary()



mnm_d <- sim_d[, 1:5] |> as.matrix()

# Plot
mgm_all <- mgm(data=mnm_d, type=c("c", rep("g", 4)), level = c(2, rep(1, 4)),
             lambdaSel = "EBIC", lambdaGam = 0.5)

par(mfrow=c(1,1))
qgraph(mgm_all$pairwise$wadj, edge.color=mgm_all$pairwise$edgecolor, layout="circle", labels=vars,
       title="all")



# Specify PA as moderator
mnm <- mgm(data=mnm_d,
             type=c("c", rep("g", 4)),
             level = c(2, rep(1, 4)),
             moderators = 1,
             lambdaSel = "CV",
             ruleReg = "AND")


l_mgm_cond <- list()
for(g in 1:2) l_mgm_cond[[g]] <- condition(mnm, values=list("1" = g))


v_max <- rep(NA, 2)
for(g in 1:2) v_max[g] <- max(l_mgm_cond[[g]]$pairwise$wadj)



par(mfrow=c(1,2))
for(g in 1:2) {
  qgraph(input = l_mgm_cond[[g]]$pairwise$wadj, 
         edge.color = l_mgm_cond[[g]]$pairwise$edgecolor, 
         layout = "circle", mar=c(2,3,5,3),
         maximum = max(v_max), vsize=16, esize=23, 
         edge.labels  = TRUE, edge.label.cex = 3)
  mtext(text = paste0("Group ", g), line=2.5)
}




















