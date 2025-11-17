# sunfishDiet.R
# Script designed to read in the diet data (from gut content analyses) and use 
# it to test for differences in diets among the sunfish species in Montreal 
# and in Lake Opinicon. The script does a few things at first to see how diets
# can be summarised. It then performs tests to determine differences in diet using 
# permanova/adonis2. It then runs a PCA on the diet points matrix allowing us to 
# visualise the variation in diet items per species and estimate diet overlaps 
# between them as we did for the shape data.

################# 1. HOUSEKEEPING AND DATA READ IN ############
# Clearing instances and resetting environment
rm(list=ls())

# Loading libraries
{
  library(parallel)
  library(ggplot2)
  library(devtools)
  library(pairwiseAdonis)
  library(tidyverse)
  library(vegan)
  library(lessR)
  library(SIBER)
  library(ggforce)
  library(forcats)
  }

# if not installed - this should be done before hand.
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

# Setting the core allocation to run parallel computations
options(mc.cores = parallel::detectCores())
detectCores(logical = T) -> cores

# Setting working directory
setwd('/Users/denis/Library/CloudStorage/GoogleDrive-denisroy1@gmail.com/My Drive/sunfish/Scripts/Diet items/')

# Setting up the file name so we can filter later
# metfile <- "dietsdr_mtl.csv"
metfile <- "dietsdr_opi.csv"

# Read in metfile
diet<-read.csv(metfile, header = T, stringsAsFactors = T, na.strings = NA)

# Viewing data to see if read in and attach
head(diet)
attach(diet)

# Extract names of the fish used for diet based on their TPSID. 
# This is important because here the TPSID species was confirmed 
# using genetics. We can store ids in a vector
ids <- diet$NEW_TPSID
ids

# Extract and assign site factor to each sample using the IDs
site <- as.factor(substr(ids,1,3))
site

# Extract and assign species for each individual using the IDs (from genetics).
# Here, because the "opi" and "mtl" file have different length names, can use 
# a ifelse statement to filter which file is being used.
ifelse(substr(metfile,9,11) == "opi", species <- as.factor(substr(ids,7,9)),
       species <- as.factor(substr(ids,9,11)))
species

# Making sure data have no NAs
any_missing <- apply(diet[, -c(1:5)], 2, function(x) any(is.na(x)))
any_missing

################# 2. PRELIMINARIES - Unique diet items #######################
# Looking at the number of unique diet items in the different individuals
# over entire data
udi <- apply(diet[, -c(1:5)], 1, function(row) length(unique(row[row > 0])))
udi

# Combine the udi, ids, and species information from the diet data
ibi<-cbind.data.frame(udi=udi, fishid=ids, site = site, species=species)

# Calculate the mean and 95%CI of diet items in overall data 
mudi <- mean(ibi$udi)
sdudi <- sd(ibi$udi)
n <- nrow(ibi)
SEM <- sdudi/sqrt(n)

# Determine the critical value for 121 "mtl"/61 "opi" degrees of freedom.
# Use qt function to get the critical value for a 95% confidence interval
# The 95% confidence interval corresponds to 97.5th percentile for the 
# two-tailed test
critval <- qt(0.975, df = n - 1)

# Calculate the margin of error for 95% CIs
ME <- critval * SEM

# Calculate the confidence interval
CI_lower <- mudi - ME
CI_upper <- mudi + ME

# Print the overall results
CI <- c(CI_lower, CI_upper)
paste("mean # diet items =", round(mudi,2), " ", 
      "95% CI", round(CI[1],2), " ", round(CI[2],2))

# Tabulate the data to get a feel for the number of diet items per species
speciesdit <- table(ibi$species, ibi$udi)
speciesdit

# Make proportions table
speciesditp <- prop.table(speciesdit,1)
speciesditp
addmargins(speciesditp)

# Plot the data to see the frequency of the number of unique items 
# in all the samples
ggplot(ibi, aes(udi)) + 
  geom_bar(stat="count", fill = "skyblue3", col = "black") +
  labs(x = "# of unique diet items", y = "Frequency") +
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 1)) +
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, by = 5)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(linewidth = 1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3, "cm"),
        legend.position = "top",
        legend.text = element_text(size = 20),
        legend.title = element_blank())

# Plots the data a bit differently to see what the mean and spread of the # of 
# unique items per species

## WARNING
# This makes it look like we have a spread of different diet items per species 
# but these are just counts (a bit misleading)
ggplot(ibi, aes(x=species, y=udi)) +
  geom_jitter(color = "darkorchid4", size = 4, alpha = 0.7, width = 0.2, height = 0.1) +
  stat_summary(fun.data = mean_sdl, colour = "gray49", na.rm = T, size = 1, linewidth = 1.2) +
  labs(x = "", y = "# of unique diet items") + 
  scale_y_continuous(limits = c(-1, 10), breaks = seq(0, 10, by = 1)) +
  theme_classic()+
  theme(axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(linewidth = 1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3, "cm"),
        legend.position = "top",
        legend.text = element_text(size = 20),
        legend.title = element_blank())

# This shows that on average, there are more items in BLG than in 
# HYB and PKS. But this is not likely significant.

################# 3. POINTS METHOD FOR DIET #######################
# Use prop.table to calculate the proportions of diet items/ind 
# essentially devising volumetric based data and making a diet proportion (dp) 
# data frame
dp <- prop.table(as.matrix(diet[,-c(1:5)]), 1)
dp <- replace(dp, is.nan(dp), 0.000)

# Set the max number of point categories needed
ptsn<-max(udi)+1
ptsn

# Set an empty vector
pts <- vector()

# Setting the first two values to 0 and 1, respectively
pts[1]<-0
pts[2]<-1

# Short loop to make the points for up to 9 categories (# unique items)
for (a in 3:ptsn) {
  pts[a]<-pts[a-1]*2
}

# Look at the points available per individual
pts <- rev(pts)
pts

# Function replacing non-zero elements of dp with ranks
replace_with_rank <- function(row) {
  nozeroel <- row[row != 0]
  ranks <- rank(-nozeroel)
  vals <- ranks
  row[row != 0] <- vals
  return(row)
}

# Apply the function "replace_with_rank" to the transposed
# matrix of dp (diet proportions) by row
dp2 <- t(apply(dp, 1, replace_with_rank))

# rename dp2 to dp3 for use in points allocation
dp3 <- dp2

# Assign the points to the ranks of dp2 starting with 1 = 128, 2 = 64,... 
for (b in 1:length(pts)) {
  dp3[dp2 == b] <- pts[b] 
}

# Look at and compare dp2 and dp3 to make sure its worked
dp3

# rReplace the pts in dp3 to the calculated point proportion/ind
dp4 <- proportions(dp3, 1)

# replace any individuals that are nan to 0.000
dp4 <- replace(dp4, is.nan(dp4), 0.000)

# Combine dp4 matrix with fish_ID, species, and site
dp5 <- cbind.data.frame(fish_id=ids, species=species, 
                      site=site, dp4)

table(dp5$species)

# Short function to get the mean of the columns of data.frame
bycolmean <- function(x) {
  tdf <- apply(x, 2, mean)
  return(tdf)
}

# Use the by cmd to split the data by species and use
# the bycolmean to get means for each diet item (could do proportions too)
dietmeans<-by(dp5[,-c(1:3)], dp5$species, bycolmean)
dietmeans

# Return the data back into a single matrix/array
speciesdiet <- do.call(cbind, dietmeans)
speciesdiet

# Get the proportions of species (spp) diet in individuals
psppdiet<-prop.table(speciesdiet,2)
addmargins(psppdiet)

### Convert to long format ###

# Set matrix as data frame to manipulate 
psppdiet <- as.data.frame(psppdiet)

## Set the new dataframe columns to specific species names
colnames(psppdiet) <- c("BLG", "HYB", "PKS")

## Set a new column to also have the names of all the diet items
psppdiet$Row <- rownames(psppdiet)

## Use tidyverse to restructure the species diet df to a long format 
plsppdiet <- psppdiet %>%
  gather(key = "species", value = "diet", -Row)

# display the data.frame to make sure it appears as expected
plsppdiet

# Define your custom colors
cuscols <- c(
  "waterplant"= "palegreen",
  "hydrachnidia"= "yellow",
  "cladoceran"= "thistle",
  "ostracoda" = "darkorchid",
  "nematode" = "darkgreen",
  "crustacean" = "cyan",
  "rotifera" = "skyblue",
  "annelida" = "firebrick4",
  "non_pred_larv" = "#E1A76B",
  "pred_larv" = "#654321",
  "malacostraca" = "darkorange",
  "adult_insect" = "dodgerblue2",
  "eggs" = "deeppink2",
  "gastropods" = "bisque4",
  "bivalve" = "navyblue")

diet_labels <- c(
  "waterplant" = "Aquatic Plants",
  "hydrachnidia" = "Water Mites",
  "cladoceran" = "Water Fleas",
  "ostracoda" = "Class Ostracoda",
  "nematode" = "Roundworms",
  "crustacean" = "Unidentified Crustaceans",
  "rotifera" = "Planktonic wheel animals",
  "annelida" = "Segmented worms",
  "non_pred_larv" = "Non-predatorial insect larvae",
  "pred_larv" = "Predatorial insect larvae",
  "malacostraca" = "Amphipods and Isopods",
  "adult_insect" = "Adult Insects",
  "eggs" = "Eggs",
  "gastropods" = "Snails",
  "bivalve" = "Mussels and Clams"
)

### WARNING ### Anotation of the superscripts on top depend on the 
# DIET ITEM ANALYSIS below. Would Run that first to adjust this.
ggplot(plsppdiet, aes(fill=fct_relevel(Row, "waterplant","bivalve",
                                       "gastropods",
                                       "non_pred_larv",
                                       "pred_larv",
                                       "malacostraca",
                                       "ostracoda",
                                       "crustacean",
                                       "hydrachnidia",
                                       "cladoceran",
                                       "rotifera",
                                       "eggs",
                                       "annelida",
                                       "nematode"), y=diet, x=species)) + 
  geom_bar(position="stack", stat="identity")+
  labs(x = "Species", y = "% volume of items in stomachs")+
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1.0, by = 0.2)) +
  scale_fill_manual(values = cuscols, labels= diet_labels) +
  annotate("text", x = 1, y = 1.1, label = paste("X"), color = "black", size = 10, family = "Courier") +
  # If using "mtl" data uncomment the following line and commment out the following 2
  #annotate("text", x = c(2,3), y = c(1.1,1.1), label = paste("Y"), color = "black", size = 10, family = "Courier") +
  annotate("text", x = 2, y = 1.1, label = paste("Y"), color = "black", size = 10, family = "Courier") +
  annotate("text", x = 3, y = 1.1, label = paste("XY"), color = "black", size = 10, family = "Courier") +
  theme(aspect.ratio = 0.8)+
  theme_classic() +
  theme(axis.title.x = element_text(size = 22),
        axis.text.x = element_text(size = 24),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 24),
        axis.line = element_line(linewidth=1.3),
        title = element_text (size = 26),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "bottom")

#################### 4. DIET ITEM ANALYSES ###############################
# Section is greatly aided by https://rpubs.com/hafezahmad/948799

# To use PERMANOVA/adonis need to convert data to dissimilarity matrix
# Use vegdist (vegan) to convert to Bray-Curtis dissimilarity
sfddis <- vegdist(dp5[,-c(1:3)], method='bray', na.rm = T)
sfddis <- as.matrix(sfddis)
dim(sfddis)

# Run the pairwise.adonis2 to test pairwise differences in diets
# Stronger alternative to ANOSIM which is not strong
species <- dp5$species
species_dd <- pairwise.adonis2(sfddis ~ species, data = dp5, 
                               permutations = 10000, sqrt.dist = T, add = T)

# Look at pairwise results 
species_dd

#################### 5. PCA OF DIET DATA ###############################
# Here, we need to perform a PCA on the diet data to quantify the overlap
# in diets among the species as we did for the shape.

# For the analyses we want just those cases without NA
diet <- diet[complete.cases(diet),]

# Setting just the numerical values into a matrix
dietmat <- as.matrix(diet[,-c(1:5)])

# Looking at item frequency in the data
itemfrq <- table(unlist(dietmat))

# Make a quick barplot to show that most of the data are zeros
barplot(itemfrq, las = 1, col = grey(5:0/5),
        xlab = "Abundance class", ylab = "Frequency")

# The extremely high proportion of 0s means we should transform 
# data to account for them

# But first look to see how many '0' we have
num0 <- length(which(dietmat == 0))
num0

# And get the proportion of '0' 
prop0 <- num0/(dim(dietmat)[1] * dim(dietmat)[2])
prop0

# so, 73% MTL / 68% OPI of data are zeros
# Patterns in the data would be obscured by this.

# Here, freq distribution of diet abundances wout '0's. 
fd_zcdiet <- apply(dietmat > 0, 1, sum)

# Barplot these quickly to verify
barplot(fd_zcdiet, main = "",
        xlab = "individual fish", ylab = "diet diversity",
        col = "skyblue3", las = 1, 
        ylim = c(0,15), cex.names = 1.5,
        cex.lab = 1.5)

# so, need to take high 0 frequency into account

# Transform data with decostand function from the vegan package in R
# using the option Hellinger which accounts for many 0 in data 
zcdiet <- decostand(dietmat, method = "hellinger")

# Run principal components analysis on these data
zcd_pca <- prcomp(zcdiet, retx = T, center = T, scale = F)
sumpca <- summary(zcd_pca)

pv <- sumpca$importance[2,]
print(pv)

pcs<-as.factor(to("PC", length(pv)))
pcdf<-cbind.data.frame(PCs=pcs, pervar=pv*100)

pcdf <- pcdf[1:20,]

# Plotting scree with ggplot
ggplot(pcdf)+
  geom_col(aes(y=pervar, x=PCs), fill = "firebrick4", col = "black") + 
  labs(x = "PCs", y = "% variance explained")+
  theme(aspect.ratio = 0.8)+
  scale_y_continuous(limits = c(0, 35), breaks = seq(0, 35, by = 5)) +
  scale_x_discrete(drop = T) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 18),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.1),
        axis.ticks.length = unit(.2,"cm"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# This will set the values to display in our overall species PCA plots
pc1_label <- paste0("PC1 (", round(pcdf$pervar[1],1), "%)")
pc2_label <- paste0("PC2 (", round(pcdf$pervar[2],1), "%)")
pc3_label <- paste0("PC3 (", round(pcdf$pervar[3],1), "%)")

# extract scores in the PCA and add them to the diet data
# zcd_pca$x

# diet$pc1 <- zcd_pca$x[,1]
# diet$pc2 <- zcd_pca$x[,2]
# diet$pc3 <- zcd_pca$x[,3]

# Here we use zcd_pca as the object and add species to the end of the data frame
pcscr <- cbind.data.frame(zcd_pca$x, species = species)  

loadings <- as.data.frame(zcd_pca$rotation[, 1:3])  # first 3 PCs
loadings$varname <- rownames(loadings)

loadings <- loadings[abs(loadings$PC1 >= 0.08) | abs(loadings$PC2 >= 0.08),]
loadings2 <- loadings[abs(loadings$PC1 >= 0.08) | abs(loadings$PC3 >= 0.08),]
loadings3 <- loadings[abs(loadings$PC2 >= 0.08) | abs(loadings$PC3 >= 0.08),]


#################### 6. PCA PLOTTING ###############################

# Setting colours for pcs plotting for each species:
sppcol <- c('PKS' = 'darkorange2', 'BLG' = 'dodgerblue4', 'HYB' = 'springgreen4')

# Plot the data (on new shape axes) with the 95% ellipses 
# (for those with enough data) with ggplot2. This will use the "sppcol"
# assigned above and the "pc1_label" and "pc2_label" from above.
ggplot(pcscr, aes(x = PC1, y = PC2, colour=species, shape = species)) + 
  geom_point(size = 5, stroke = 1, alpha = 0.7) + 
  stat_ellipse(aes(fill = species), geom = "polygon", color = NA, alpha = 0.2, 
               level = 0.95, type = "norm") +
  stat_ellipse(aes(color = species), geom = "path", linewidth = 1, 
               level = 0.95, type = "norm") +
  scale_fill_manual(values = sppcol) +
  labs(x = pc1_label, y = pc2_label) + 
  scale_color_manual(values = sppcol) +
  geom_hline(yintercept = 0, linewidth = 1, color="firebrick4", linetype="dotted") +
  geom_vline(xintercept = 0, linewidth = 1, color="firebrick4", linetype="dotted") +
  scale_x_continuous(limits = c(-1.2, 1.2), breaks = seq(-1.0, 1.0, by = 0.25)) +
  scale_y_continuous(limits = c(-1.2, 1.2), breaks = seq(-1.0, 1.0, by = 0.25)) +
  theme(aspect.ratio = 0.70) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "lightgrey"),
    axis.text.x = element_text(size = 22, angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(linewidth=1.3),
    axis.ticks = element_line(linewidth = 1.3),
    axis.ticks.length = unit(.3,"cm"),
    legend.text = element_text(size = 17),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.justification = "center") +
    geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm")), 
               colour = "black", linewidth = 1.2, inherit.aes = FALSE) +
    geom_text(data = loadings, aes(x = PC1, y = PC2, label = varname),
            color = "black", size = 6, vjust = 1.5, inherit.aes = FALSE)


# As above for PC1 vs. PC3
ggplot(pcscr, aes(x = PC1, y = PC3, colour=species, shape = species)) + 
  geom_point(size = 5, stroke = 1, alpha = 0.7) + 
  stat_ellipse(aes(fill = species), geom = "polygon", color = NA, alpha = 0.2, 
               level = 0.95, type = "norm") +
  stat_ellipse(aes(color = species), geom = "path", linewidth = 1, 
               level = 0.95, type = "norm") +
  scale_fill_manual(values = sppcol) +
  labs(x = pc1_label, y = pc3_label) + 
  scale_color_manual(values = sppcol) +
  geom_hline(yintercept = 0, linewidth = 1, color="firebrick4", linetype="dotted") +
  geom_vline(xintercept = 0, linewidth = 1, color="firebrick4", linetype="dotted") +
  scale_x_continuous(limits = c(-1.2, 1.2), breaks = seq(-1.0, 1.0, by = 0.25)) +
  scale_y_continuous(limits = c(-1.0, 1.0), breaks = seq(-1.0, 1.0, by = 0.25)) +
  theme(aspect.ratio = 0.70) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "lightgrey"),
    axis.text.x = element_text(size = 22, angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(linewidth=1.3),
    axis.ticks = element_line(linewidth = 1.3),
    axis.ticks.length = unit(.3,"cm"),
    legend.text = element_text(size = 17),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.justification = "center") +
  geom_segment(data = loadings2, aes(x = 0, y = 0, xend = PC1, yend = PC3),
               arrow = arrow(length = unit(0.25, "cm")), 
               colour = "black", linewidth = 1.2, inherit.aes = FALSE) +
  geom_text(data = loadings2, aes(x = PC1, y = PC3, label = varname),
            color = "black", size = 6, vjust = 1.5, inherit.aes = FALSE)

# As above for PC2 vs. PC3
ggplot(pcscr, aes(x = PC2, y = PC3, colour=species, shape = species)) + 
  geom_point(size = 5, stroke = 1, alpha = 0.7) + 
  stat_ellipse(aes(fill = species), geom = "polygon", color = NA, alpha = 0.2, 
               level = 0.95, type = "norm") +
  stat_ellipse(aes(color = species), geom = "path", linewidth = 1, 
               level = 0.95, type = "norm") +
  scale_fill_manual(values = sppcol) +
  labs(x = pc2_label, y = pc3_label) + 
  scale_color_manual(values = sppcol) +
  geom_hline(yintercept = 0, linewidth = 1, color="firebrick4", linetype="dotted") +
  geom_vline(xintercept = 0, linewidth = 1, color="firebrick4", linetype="dotted") +
  scale_x_continuous(limits = c(-1.2, 1.2), breaks = seq(-1.0, 1.0, by = 0.25)) +
  scale_y_continuous(limits = c(-1.0, 1.0), breaks = seq(-1.0, 1.0, by = 0.25)) +
  theme(aspect.ratio = 0.70) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "lightgrey"),
    axis.text.x = element_text(size = 22, angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(linewidth=1.3),
    axis.ticks = element_line(linewidth = 1.3),
    axis.ticks.length = unit(.3,"cm"),
    legend.text = element_text(size = 17),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.justification = "center") +
  geom_segment(data = loadings3, aes(x = 0, y = 0, xend = PC2, yend = PC3),
               arrow = arrow(length = unit(0.25, "cm")), 
               colour = "black", linewidth = 1.2, inherit.aes = F) +
  geom_text(data = loadings3, aes(x = PC2, y = PC3, label = varname),
            color = "black", size = 6, vjust = 1.5, inherit.aes = F)

############################# 7.ELLIPSE OVERLAP CALCULATIONS ############################
# Using the SIBER package which can estimate the overlap among different ellipses
# generated from various data sets. The results are very comparable with the 
# CAR/sf packages but do not assume normal distribution or even sample sizes.

# First need to make data frame of PCs used and include species in the groups.
# could try to use community to reflect mtl vs opn. but try this for now.
ovdpc_12 <- cbind.data.frame(iso1 = pcscr$PC1, iso2 = pcscr$PC2,
                           group = pcscr$species, community = rep(1,nrow(pcscr)))

# Second, convert data frame to SIBER object
ovdpc12_so <- createSiberObject(ovdpc_12)

# Cmd sets the parameters for the 95% ellipses in the plotSiberObject 
# below.
group.ellipses.args  <- list(n = 1000, p.interval = 0.95, lty = 1, lwd = 2)

# Here, SIBER plots the data as well - which helps check if we get same patterns 
# as biplots (which are prettier).

#par(mfrow=c(1,1))
plotSiberObject(ovdpc12_so, ax.pad = 1, hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args, 
                group.hulls = F, group.hull.args,
                bty = "L", las = 1, iso.order = c(1,2),
                xlab = pc1_label,
                ylab = pc2_label,
                x.limits = c(-1.5,1.5),
                y.limits = c(-1.5,1.5),
                cex = 0.5)

# The graph produced here should look the same. If they are different - there is a 
# problem. For "opi" data, SIBER will generate ellipse for PKS, but it is not likely
# robust

# Calculate the summary statistics for each group 
# Total Area - TA
# Standard Ellipse Area - SEA
# Standard Ellipse Area for small sample sizes - SEAc

# once plotted can look at size of each ellipse and determine 
# amount of overlap
ML12 <- groupMetricsML(ovdpc12_so) 
ML12

# Identify the ellipses
elp1 <-"1.BLG"
elp2 <-"1.HYB"
elp3 <-"1.PKS"

# calculate the maximum likelihood overlap of the two ellipses considered

# First, BLG - HYB
BLG_HYB <- maxLikOverlap(elp1, elp2, ovdpc12_so, p.interval = 0.95, n = 1000)
BLG_HYB

# Calculating overlap between hybrids and bluegill
blghyb_ov <- BLG_HYB[3] / BLG_HYB[1]
print(paste("bluegill hybrid ov: ",blghyb_ov))

# Then, PKS - HYB
PKS_HYB <- maxLikOverlap(elp3, elp2, ovdpc12_so, p.interval = 0.95, n = 1000)
PKS_HYB

# Now between hybrids and pumpkinseed
pkshyb_ov <- PKS_HYB[3] / PKS_HYB[1]
print(paste("pumpkinseed hybrid ov: ",pkshyb_ov))

# Finally, between PKS - BLG
BLG_PKS <- maxLikOverlap(elp1, elp3, ovdpc12_so, p.interval = 0.95, n = 1000)
BLG_PKS

# Calculating overlap between hybrids and bluegill
blgpks_ov <- BLG_PKS[3] / BLG_PKS[1]
pksblg_ov <- BLG_PKS[3] / BLG_PKS[2]

print(paste("bluegill pumpkinseed ov: ", blgpks_ov))
print(paste("pumpkinseed blugill ov: ", pksblg_ov))

# "mtl": Results show greater overlap between between hybrids and pumpkinseed, versus 
# hybrids and bluegill (as in shape). It also shows that the hybrids are much 
# more variable than the bluegill, but comparable to pumpkinseed. The hybrids
# are shifted a bit more to the upper left beyond pumpkinseed - suggestive maybe
# some transgressive segregation in that direction, which may be related to 
# exploiting more water plants and adult insects.

# "opi": Results suggest that the bluegills have some distinct diet items that differ
# considerably from both HYB and PKS. The lack of PKS numbers in our samples reflect 
# the increasing difficulty in collecting actual true PKSs from this region. All observed
# PKSs are well within the HYB ellipse and suggests HYB have to some extent replaced the 
# PKSs in this system as they feed on essentially the same diet items. The overlap 
# calculation may be misleading because it is difficult to get a robust ellipse with only 
# 3 data points. 

############################# 8. ELLIPSE ECCENTRICITY ############################
# Calculating ellipse areas and eccentricities with 95% CIs

# split data by species
sppdata <- split(pcscr, pcscr$species)

# prepare empty data frame for results
elstat <- data.frame(
  species = character(length(sppdata)),
  elA = numeric(length(sppdata)),
  elAlwr = numeric(length(sppdata)),
  elAupr = numeric(length(sppdata)),
  elX = numeric(length(sppdata)),
  elXlwr = numeric(length(sppdata)),
  elXupr = numeric(length(sppdata)),
  stringsAsFactors = FALSE
)

# scaling factor for 95% CI ellipse
c95 <- sqrt(qchisq(0.95, 2))

for (a in seq_along(sppdata)) {
  tmp <- sppdata[[a]]
  
  # covariance matrix of Comp1 and Comp2
  covm <- cov(cbind(tmp$PC1, tmp$PC2))
  evals <- eigen(covm, symmetric = T)$values
  evals <- sort(evals, decreasing = T) # λ1 >= λ2
  
  # Generate semi-axes for the ellipses
  aaxis <- c95 * sqrt(evals[1])
  baxis <- c95 * sqrt(evals[2])
  
  # Calculate the area A of the ellipses
  elA <- pi * aaxis * baxis
  
  # Calculate the eccentricity X
  elX <- if (evals[1] > 0) sqrt(1 - evals[2] / evals[1]) else NA
  
  # This is a bootstrapping function to used to estimate the SEs 
  # for the A and X calculations
  nboot <- 1000
  elit <- matrix(NA, nboot, 2)
  for (b in 1:nboot) {
    xboot1 <- sample(tmp$PC1, replace = T)
    xboot2 <- sample(tmp$PC2, replace = T)
    tcov <- cov(cbind(xboot1, xboot2))
    teval <- sort(eigen(tcov, symmetric = T)$values, decreasing = T)
    ta <- c95 * sqrt(teval[1])
    tb <- c95 * sqrt(teval[2])
    balt <- pi * ta * tb
    becc <- if (teval[1] > 0) sqrt(1 - teval[2] / teval[1]) else NA
    elit[b, ] <- c(balt, becc)
  }
  bse <- apply(elit, 2, sd, na.rm = T)
  
  # Calculate 95% CI ranges for both A and X
  elAci <- c(elA - 1.96 * bse[1], elA + 1.96 * bse[1])
  elXci <- c(elX - 1.96 * bse[2], elX + 1.96 * bse[2])
  
  # Store results (rounded) to use in table
  elstat[a, ] <- list(
    species = names(sppdata)[a],
    elar = round(elA, 5),
    elarl = round(elAci[1], 5),
    elaru = round(elAci[2], 5),
    elec = round(elX, 5),
    elecl = round(elXci[1], 5),
    elecu = round(elXci[2], 5)
  )
}

# display results
elstat

#### END ####

