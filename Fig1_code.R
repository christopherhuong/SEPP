
library(huge)
library(qgraph)
library(bootnet)
library(ggplot2)
library(gridExtra)


setwd("G:/My Drive/INCH Projects/[2023] SEPP Special Issue/datacode")


# Load in data
stretch <- read.csv('StretchingGroupData.csv')
cycle <- read.csv('CyclingGroupData.csv')

# Nonparanormal transformation
stretch <- huge.npn(stretch, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE)
stretch <- as.data.frame(apply(stretch, use='pairwise.complete.obs', 2, as.numeric))

cycle <- huge.npn(cycle, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE)
cycle <- as.data.frame(apply(cycle, use='pairwise.complete.obs', 2, as.numeric))


names = c("Perserv", "Neg", "Critic", "Brood", "Replay")
nodenames =c("Perseveration; “How much have you been thinking about your performance?”",
             "Negativity; “How negative have your thoughts been?",
             "Self-criticism; “How much have you criticized yourself about your performance?",
             "Brooding; “How much have you thought about how upset you felt?”",
             "Replaying; “To what extent did you replay parts of what happened in your mind?")

# Estimate glasso networks with EBIC hyperparameter set to 0.5
stretch_glasso <- estimateNetwork(stretch, default = "EBICglasso", tuning = .5)
cycle_glasso <- estimateNetwork(cycle, default = "EBICglasso", tuning = .5)

L <- averageLayout(stretch_glasso, cycle_glasso)

# Plot
stretch_plot <- qgraph(stretch_glasso$graph,
                       title = "Stretching",
                       layout=L, labels=names, theme = "colorblind")


cycle_plot <- qgraph(cycle_glasso$graph, 
                     title = "Cycling",
                     layout=L, labels=names, theme = "colorblind")

# Strength centrality plots
all_cen <- list("Stretching"=stretch_glasso, "Cycling"=cycle_glasso)
cen_plot <- centralityPlot(all_cen, include = "Strength", labels = names, scale = "z-scores") +
  xlab("z-scores")

# Assess accuracy of edge weights using bootstrapped confidence intervals
# boot_stretch <- bootnet(stretch_glasso, type = "nonparametric", nBoots = 1000)
# boot_cycle <- bootnet(stretch_glasso, type = "nonparametric", nBoots = 1000)
# 
# bootplot_stretch <- plot(boot_stretch, labels = T, order = "id") +
#   ggtitle("Stretching")
# 
# bootplot_cycle <- plot(boot_cycle, labels = F, order = "id") +
#   ggtitle("Cycling") 
# 



# Figures  ------------------------------------------------------------------

# Networks
pdf(file="Fig1_A.pdf", width = 8, height = 4)
par(mfrow = c(1,2))
plot(stretch_plot)
plot(cycle_plot)
dev.off()

# Centrality plot
pdf(file="Fig1_B.pdf", width=4, height=4)
plot(cen_plot)
dev.off()

# Bootstrapped edge-weights
# pdf(file="Fig2_A.pdf", width=8, height=4)
# grid.arrange(bootplot_stretch, bootplot_cycle,
#              ncol=2, widths=c(2, 1.7))
# dev.off()
# 











