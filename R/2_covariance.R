# 0 - Load libraries and data ####

library(readr)
library(Morpho)
install.packages("evolqg")
library(evolqg) # random skewers covariance test
library(dplyr)
library(ggplot2)
install.packages("ggforce")
library(ggforce)
library(ggrepel)
library(glue) # makes string manipulation work like python

# load feesh
load("Zfish_Shh_GMM/Zfish_Shh_GMM.RData")
Zfish_Classifiers <- read_csv("./data/Zfish_Classifiers.csv")

# check dims, and assume everything is in order
dim(PCA_df)
dim(Zfish_Classifiers)

# 1. Histograms - Aponte ####

# do some regression
treatmentLabels <- unique(Zfish_Classifiers$Treatment)
Zfish_Classifiers$Treatment <- factor(Zfish_Classifiers$Treatment, levels = treatmentLabels[c(6, 2:5, 1)])

Zfish_Classifiers$Treatment
# PCA_df <- matrix(as.numeric(PCA_df), nrow = 117, ncol = 117, byrow = T)

# factor level regression ####
bigLM <- lm(as.matrix(PCA_df[,-117]) ~ Zfish_Classifiers$Treatment)

factorScores <- RegScore(bigLM)

boxplot(factorScores[,5] ~ Zfish_Classifiers$Treatment)

# continuous####
continuousTreatment <- rep(0, 117)
continuousTreatment[Zfish_Classifiers$Treatment  == "DMSO"] <- 0
continuousTreatment[Zfish_Classifiers$Treatment  == "Cyc20"] <- 20
continuousTreatment[Zfish_Classifiers$Treatment  == "Cyc40"] <- 40
continuousTreatment[Zfish_Classifiers$Treatment  == "Cyc60"] <- 60
continuousTreatment[Zfish_Classifiers$Treatment  == "Cyc80"] <- 80
continuousTreatment[Zfish_Classifiers$Treatment  == "Cyc100"] <- 100
as.numeric(continuousTreatment)

bigContinuousLM <- lm(as.matrix(PCA_df[,-117]) ~ continuousTreatment)

continuousScores <- RegScore(bigContinuousLM)

#boxplot(continuousScores ~ Zfish_Classifiers$Treatment)
ggDf <- data.frame(scores = continuousScores, treatment = Zfish_Classifiers$Treatment)
ggDf$treatment <- factor(ggDf$treatment, levels = c("DMSO", "Cyc20", "Cyc40", "Cyc60", "Cyc80", "Cyc100"))
treatmentColors <- c("DMSO" = "black", "Cyc20" = "blue", "Cyc40" = "purple", "Cyc60" = "pink", "Cyc80" = "orange", "Cyc100" = "red")

p <- ggplot(ggDf, aes(x = treatment, y = scores, fill = treatment)) + 
  geom_boxplot() +
  scale_fill_manual(values = treatmentColors) + 
  labs(x = "Cyclopamine level (units)", y = "Cycopamine effect\nshape scores") +
  theme_classic() +
  theme(legend.position = "NA")

print(p)
ggsave(filename = "./figs/cycRegressionScores.pdf", plot = p, height = 6.5, width = 10, units = "cm", dpi = "print")


# how do we test the covariance structure ####
covDMSO <- cov(PCA_df[Zfish_Classifiers$Treatment  == "DMSO",-117])
cov20 <- cov(PCA_df[Zfish_Classifiers$Treatment  == "Cyc20",-117])
cov40 <- cov(PCA_df[Zfish_Classifiers$Treatment  == "Cyc40",-117])
cov60 <- cov(PCA_df[Zfish_Classifiers$Treatment  == "Cyc60",-117])
cov80 <- cov(PCA_df[Zfish_Classifiers$Treatment  == "Cyc80",-117])
cov100 <- cov(PCA_df[Zfish_Classifiers$Treatment  == "Cyc100",-117])

#not that interesting. All treatment levels are different from control. 
RandomSkewers(covDMSO, cov100)

# Measure the dispersion around each group vector 
mahalAngle <- function(PCA, refGruppo = "DMSO", tarGruppo, nameList = Zfish_Classifiers$Treatment, refVec = NULL){
  
  # standardize PC scores 
  standardizedScores <- PCA$x
  for(i in 1:ncol(PCA$x)) standardizedScores[,i] <-  standardizedScores[,i]/PCA$sdev[i]
  
  standardizedScores <- PCA$x
  
  # calculate the mean difference vector
  refMean <- colMeans(standardizedScores[nameList  == refGruppo,])
  tarMean <- colMeans(standardizedScores[nameList  == tarGruppo,])
  diffVec <- refMean - tarMean
  
  if(is.null(refVec) == F) diffVec <- refVec
  
  # calculate the angle between each individual and the vector
  angleValues <- rep(NA, sum(nameList == tarGruppo))
  for(i in 1:length(angleValues)){
    tarInd <- standardizedScores[which(nameList  == tarGruppo)[i],]
    diffInd <- refMean - tarInd
    angleValues[i] <- Morpho::angle.calc(diffVec, diffInd)
  }
  
  return(angleValues * (180/pi))
}

# 1.1 Dispersion around their own trajectory ####

# I didn't make a ggplot version of this b/c it's a bit confusing
# It's every group only compared to its own trajectory
# the only thing you really get from it is 100 is more consistent than everyone else
mah20 <- mahalAngle(PCA_head, tarGruppo = "Cyc20")
mah40 <- mahalAngle(PCA_head, tarGruppo = "Cyc40")
mah60 <- mahalAngle(PCA_head, tarGruppo = "Cyc60")
mah80 <- mahalAngle(PCA_head, tarGruppo = "Cyc80")
mah100 <- mahalAngle(PCA_head, tarGruppo = "Cyc100")

plot(density(mah100), col = adjustcolor(1, alpha.f = .3), xlim = c(0,120))
polygon(density(mah20), col = adjustcolor(2, alpha.f = .3))
polygon(density(mah40), col = adjustcolor(3, alpha.f = .3))
polygon(density(mah60), col = adjustcolor(4, alpha.f = .3))
polygon(density(mah80), col = adjustcolor(5, alpha.f = .3))

# 1.2 Dispersion around each treatment level trajectory ####
treatmentLevel <- c('20', '40', '60', '80', '100')
for(i in treatmentLevel){
  refMean <- colMeans(PCA_head$x[Zfish_Classifiers$Treatment  == "DMSO",])
  tarMean <- colMeans(PCA_head$x[Zfish_Classifiers$Treatment  == glue('Cyc{i}'),])
  diffVec <- refMean - tarMean
  
  mah20 <- mahalAngle(PCA_head, tarGruppo = "Cyc20", refVec = diffVec)
  mah40 <- mahalAngle(PCA_head, tarGruppo = "Cyc40", refVec = diffVec)
  mah60 <- mahalAngle(PCA_head, tarGruppo = "Cyc60", refVec = diffVec)
  mah80 <- mahalAngle(PCA_head, tarGruppo = "Cyc80", refVec = diffVec)
  mah100 <- mahalAngle(PCA_head, tarGruppo = "Cyc100", refVec = diffVec)
  
  # make it pretty with ggplot
  ggDf <- data.frame(angle = c(mah20, mah40, mah60, mah80, mah100), treatment = c(rep("20", length(mah20)), rep("40", length(mah40)), rep("60", length(mah60)), rep("80", length(mah80)), rep("100", length(mah100))))
  ggDf$treatment <- factor(ggDf$treatment, levels = treatmentLevel)
  treatmentColors <- c("20" = "blue", "40" = "purple", "60" = "pink", "80" = "orange", "100" = "red")
  
  p <- ggplot(ggDf, aes(x = angle, group = treatment, fill = treatment)) +
    geom_density(alpha = 0.75) + 
    xlim(c(0,140)) +
    scale_fill_manual(values = treatmentColors) + 
    labs(x = glue("Angle to cyclopamine\n {i} trajectory"), y = "Density", fill = "Cylopamine\ntreatment") +
    theme_classic() +
    theme(legend.position = "top")
  
  print(p)
  ggsave(filename = glue("figs/angle_plots/cyc{i}_angle_histogram.pdf"), plot = p, height = 6.5, width = 10, units = "cm", dpi = "print")
  
}

# 2. Angle plots - Cass ####

# Reference shape
refMean = colMeans(PCA_head$x[Zfish_Classifiers$Treatment  == "DMSO",])
# Target groups
test_groups = paste0("Cyc", seq(20,100,20))

# Loop, calculate angles, and plot
for (i in 1:length(test_groups)) {
  # Select target
  tarMean = colMeans(PCA_head$x[Zfish_Classifiers$Treatment  == test_groups[i],])
  # Create axis
  diffVec = refMean - tarMean
  # Calculate angles for each group
  DMSO = mahalAngle(PCA_head, tarGruppo = "DMSO", refVec = diffVec)
  mah20 = mahalAngle(PCA_head, tarGruppo = "Cyc20", refVec = diffVec)
  mah40 = mahalAngle(PCA_head, tarGruppo = "Cyc40", refVec = diffVec)
  mah60 = mahalAngle(PCA_head, tarGruppo = "Cyc60", refVec = diffVec)
  mah80 = mahalAngle(PCA_head, tarGruppo = "Cyc80", refVec = diffVec)
  mah100 = mahalAngle(PCA_head, tarGruppo = "Cyc100", refVec = diffVec)
  # Join in data frame
  angles = data.frame(Treatment = c(rep("DMSO", length(DMSO)), rep("Cyc20", length(mah20)), rep("Cyc40", length(mah40)), rep("Cyc60", length(mah60)), rep("Cyc80", length(mah80)), rep("Cyc100", length(mah100))), 
                      degrees = c(DMSO, mah20, mah40, mah60, mah80, mah100))
  # Convert to radians
  angles = angles %>% mutate(radians = degrees*(pi/180))
  
  # Calculate means
  mean_degrees = tapply(angles$degrees, angles$Treatment, mean, na.rm = TRUE)
  mean_radians = tapply(angles$radians, angles$Treatment, mean, na.rm = TRUE)
  # Combine results into a data frame
  mean_angles = data.frame(Treatment = names(mean_degrees), mean_degrees = mean_degrees, mean_radians = mean_radians)
  # Plot radius
  rad = 0.5
  # Calculate x and y coordinates
  angles = angles %>% mutate(x = cos(radians)*(rad), 
                             y = sin(radians)*(rad))
  mean_angles = mean_angles %>% 
    mutate(x = cos(mean_radians)*(rad+0.05), 
           y = sin(mean_radians)*(rad+0.05))
  # Add labels
  labels = mean_angles %>% select(Treatment, mean_degrees, mean_radians) %>% 
    mutate(x = cos(mean_radians)*(rad+0.1), 
           y = sin(mean_radians)*(rad+0.1),
           label = paste0(format(round(mean_degrees, digits = 1), nsmall = 1), "°"))
  # Score angle plot
  angle_plot = ggplot(mean_angles) +
    geom_arc(aes(x0 = 0, y0 = 0, r = rad, start = -0.5*pi, end = 0.5*pi), color = "#9E9E9E") +
    xlim(-rad-0.25, rad+0.25) +
    ylim(-0.05, rad+0.15) +
    annotate("segment", x = 0, y = 0, xend = 0, yend = rad+0.1, lineend = "round", color = "#9E9E9E") +
    annotate("segment", x = -rad-0.1, y = 0, xend = rad+0.1, yend = 0, lineend = "round", color = "#9E9E9E") +
    geom_segment(aes(x = 0, y = 0, xend = x, yend = y, col = Treatment), linewidth = 0.75, lineend = "round") +
    #geom_point(data = angles, aes(x = x, y = y, col = Treatment), alpha = 0.6, size = 1.5) +
    #geom_text_repel(data = labels, aes(x = x, y = y, label = label), colour = c("red", "blue", "purple", "pink", "orange", "black"), 
     #               point.padding = 0, hjust = 0, min.segment.length = 5, size = 3, force = 0.0005, fontface = "bold", direction = "y", box.padding = 0.2, ) +
    annotate("text", label = "0°", x = 0.65, y = 0, size = 3, colour = "#9E9E9E") +
    #annotate("text", label = "90°", x = 0, y = 0.65, size = 3, colour = "#9E9E9E") +
    annotate("text", label = "180°", x = -0.7, y = 0, size = 3, colour = "#9E9E9E") +
    scale_colour_manual(values = c("red", "blue", "purple", "pink", "orange", "black")) +
    theme_light(base_size = 9) +
    theme(
      legend.position = "none",
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 9),
      legend.background = element_rect(fill = NA),
      legend.key = element_rect(fill = NA),
      panel.border = element_rect(color = "#E0E0E0"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank())
  print(angle_plot)
  ggsave(filename = paste0("figs/angle_plots/", tolower(test_groups[i]), "_anchored_angle_plot.pdf"), plot = angle_plot, height = 4, width = 7.5, units = "cm", dpi = "print")
}



