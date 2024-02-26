
library(huge)
library(qgraph)
library(bootnet)


setwd("G:/My Drive/INCH Projects/[2023] SEPP Special Issue/datacode")



# Figure 1 ----------------------------------------------------------------


stretch <- read.csv('StretchingGroupData.csv')
cycle <- read.csv('CyclingGroupData.csv')


stretch <- huge.npn(stretch, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE)
stretch <- as.data.frame(apply(stretch, use='pairwise.complete.obs', 2, as.numeric))

cycle <- huge.npn(cycle, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE)
cycle <- as.data.frame(apply(cycle, use='pairwise.complete.obs', 2, as.numeric))


names = c("Perserv", "Neg", "Critic", "Brood", "Replay")


stretch_glasso <- estimateNetwork(stretch, default="EBICglasso", tuning=.5)
cycle_glasso <- estimateNetwork(cycle, default="EBICglasso", tuning=.5)
  
L <- averageLayout(stretch_glasso, cycle_glasso)

stretch_plot <- qgraph(stretch_glasso$graph,
                       title = "Stretching",
                       layout=L, edge.labels=FALSE,
                       labels=names, theme = "colorblind")


cycle_plot <- qgraph(cycle_glasso$graph, 
                     title = "Cycling",
                     layout=L, edge.labels=FALSE,
                     labels=names, theme = "colorblind")



all_cen <- list("Stretching"=stretch_glasso, "Cycling"=cycle_glasso)
cen_plot <- centralityPlot(all_cen, include = "Strength", labels = names, scale = "z-scores")



pdf(file="Fig1.pdf", width = 6, height = 12)
par(mfrow = c(2,1))
plot(stretch_plot)
plot(cycle_plot)
plot(cen_plot)
dev.off()



# Figure 2 ----------------------------------------------------------------
library(tidyverse)
library(qgraph)
library(mlVAR)
library(psych)


load("G:/My Drive/INCH Projects/[2023] EMA + PA/Data/combined_data_11_25_23.RData")


length(unique(d$participant)) # 15 participants from initial EMA collection
# for demonstration purposes, duplicate and combine data for N = 30 participants to return an artificially dense network

d1 <- d
d1$participant <- paste0(d1$participant, "A")

d3 <- rbind(d, d1)
length(unique(d3$participant)) # 30 participants



colnames(d3) <- c("participant", "day", "beep", "email", "time", "lesson",
                 "Tired", "Stress", "Pos", "Conc", "Sad", "Anx", "Anhed", "Inferior", "Alone",
                 "t", "PA")

# only retain variables related to PA for plotting clarity
varlabs <- c("Tired",
             # "Stress",
             "Pos",
             "Conc",
             "Sad",
             "Anx",
             "Anhed",
             # "Inferior",
             "Alone",
             "PA")

names <- c("Feeling tired or having little energy", "I found it difficult to relax",
           "Feeling positive or enthusiastic", "Trouble concentrating on things",
           "Feeling down, depressed, or hopeless", "Feeling nervous, anxious, or on edge",
           "Little interest or pleasure in doing things", "Feeling bad about yourself",
           "Feeling alone", "Physical Activity")



m1 <- mlVAR(d3,
            vars=varlabs,
            idvar="participant",
            dayvar="day",
            beepvar="beep",
            lags = 1,
            temporal = "orthogonal",
            contemporaneous = "orthogonal",
            nCores = 8)


temp <- getNet(m1, "temporal", nonsig = "hide")
cont <- getNet(m1, "contemporaneous", layout = "spring", nonsig = "hide", rule = "and")

L <- averageLayout(cont, temp)

#remove auto regressive effects for plotting clarity
diag(temp) <- 0



pdf("G:/My Drive/INCH Projects/[2023] SEPP Special Issue/datacode/Fig2.pdf", width=10, height=6)
par(mfrow=c(1,2))
temp_plot <- qgraph(temp, layout = "circle",
                    title="Temporal, Lag-1", theme='colorblind', negDashed=FALSE,
                    legend.cex=0.4, details=TRUE, legend=F, nodeNames = names)


cont_plot <- qgraph(cont, layout = "circle",
                    title="Contemporaneous", theme='colorblind', negDashed=FALSE,
                    legend.cex=0.4, details=TRUE, legend=F, nodeNames = names)
dev.off()




# Figure 3 ----------------------------------------------------------------
library(tidyverse)
library(graphicalVAR)

# import daily phq2 data
# https://www.synapse.org/#!Synapse:syn10848316/wiki/548732
phq2 <- read.csv("G:/My Drive/INCH Projects/[2023] SEPP Special Issue/datacode/phq2.csv")
# import smartphone-based GPS mobility data
# https://www.synapse.org/#!Synapse:syn10848316/wiki/587778
passive_v2 <- read.csv("G:/My Drive/INCH Projects/[2023] SEPP Special Issue/datacode/passive_v2_mobility.csv")

phq2 <- phq2 %>% select(participant_id, dt_yesterday, phq2_1, phq2_2)
colnames(phq2) <- c("id", "date", "phq1", "phq2")

passive_v2 <- passive_v2 %>% select(participant_id, dt_passive, hours_active, hours_of_sleep)
colnames(passive_v2) <- c("id", "date", "active", "sleep")

sum(is.na(phq2))# 14 NA's
phq <- na.omit(phq2)
sum(is.na(passive_v2)) #zero NA's

phq2$date <- as.Date(phq2$date, format = "%Y-%m-%d")
passive_v2$date <- as.Date(passive_v2$date, format = "%Y-%m-%d")

d <- passive_v2 %>% left_join(phq2, by = c("id", "date"))

sum(is.na(d))
d <- na.omit(d)

d <- d %>% 
  group_by(id) %>% filter(n() > 50)

d <- d %>% arrange(id, date) %>% group_by(id) %>%  
  complete(date = seq.Date(min(date), max(date), by = "day")) %>%  
  ungroup() 


vars = c("Actv", "Slp", "Depr", "Anhe")

# retain participant EN05297 for demonstration purposes
p <- d %>% filter(id == "EN05297")


p_net <- graphicalVAR(p, nLambda = 50, verbose = T, gamma = 0,
                       scale = T, vars = c("active", "sleep", "phq1", "phq2"))


# nodenames = c("The cumulative time spent in the active velocity bin",
#               "The estimated hours of sleep accrued during the previous night",
#               "Yesterday, were you bothered by any of the following problems: feeling down, depressed, or hopeless",
#               "Yesterday, did you have little interest or pleasure in doing things?")

pdf(file="G:/My Drive/INCH Projects/[2023] SEPP Special Issue/datacode/Fig3.pdf", width = 6, height = 6)
plot(p_net, "PDC", layout="circle", labels = vars, theme = "colorblind")
dev.off()








