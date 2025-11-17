# sunfishGM.R
# This script reads in the TPS file for a collection of sunfish (bluegill,
# pumpkinseed, and their hybrids) that have been genetically identified 
# using µsats. The TPS file has been formatted to use the ind. IDs which
# includes their genetic affinity (here called species). We are interested in 
# Assessing their shape variation, testing whether there are any significant 
# differences among and between the species, and finally, assessing the overlap 
# of the shape space among the species. Ultimately, we want to know how the 
# pumpkinseed - bluegill - hybrids overlap to test hypotheses related to 
# hybridisation among the groups.

#### 1. PRELIMINARIRES ####

# This cmd clears and resets the R instance and removes all data and 
# previously loaded variables.
rm(list = ls())

# Next load libraries. These are the list of cmds that are specific 
# to a set of analyses. In this case we have a few libraries to load 
# We can load them all at once using curly brackets.

{
  library(geomorph, warn.conflicts = F, quietly = T)
  library(ggplot2, warn.conflicts = F, quietly = T)
  library(ggfortify, warn.conflicts = F, quietly = T)
  library(ellipse, warn.conflicts = F, quietly = T)
  library(car, warn.conflicts = F, quietly = T)
  library(dplyr, warn.conflicts = F, quietly = T)
  library(vegan, warn.conflicts = F, quietly = T)
  library(lessR, warn.conflicts = F, quietly = T)
  library(emmeans, warn.conflicts = F, quietly = T)
  library(MASS)
  library(rstatix)
  library(effectsize)
  library(RRPP)
  library(parallel)
  library(SIBER)
}

# Setting the core allocation to run parallel computations
options(mc.cores = parallel::detectCores())
cores <- detectCores(logical = T)

# Get working directory
setwd("/Users/denis/Library/CloudStorage/GoogleDrive-denisroy1@gmail.com/My Drive/sunfish/Scripts/GM/")

# Here we load the TPS file we are interested in using. But,
# need to name the file to run and filter later on.
#tpsfile <- "sunfishopi.TPS"
tpsfile <- "sunfishmtl.TPS"

tpsdat <- readland.tps(tpsfile, specID = "ID", readcurves = T, 
                       negNA = FALSE, warnmsg = T)

# Need to apply a filter here too.
#metfile <- "sunfishopi_meta.csv"
metfile <- "sunfishmtl_meta.csv"

# Import the covariate file, containing METADATA for the sunfish.
metdat <- read.csv(metfile, header = T, na.strings = "NA", 
                    stringsAsFactors = T)

# Extract names of the fish we landmarked and store in a vector
idsall <- dimnames(tpsdat)[[3]]
idsall

# Extract and assign site factor to each sample using the IDs
site <- as.factor(substr(idsall,1,3))
site

# Extract and assign species assignment for each individual using the IDs
ifelse(substr(metfile,8,10) == "opi", species <- as.factor(substr(idsall,7,9)),
       species <- as.factor(substr(idsall,9,11)))
species

# Checking to see if the metdat is in same order as tpsdat. 
# If the missass is empty/null/integer(0), then this is good.
misass <- which(idsall %in% as.character(metdat$NEW_TPSID) == "FALSE")
misass

# Build a covar dataframe from the IDs, lengths, and weights.
covar <- cbind.data.frame(site = site, species = species, 
                          lencm = metdat$lencm, weight = metdat$weightg)

# Calculate the summary statistics for the covars to report
# in a table. Mostly interested by species rather than both 
# site and species. But, can try both.
covsumlw <- covar %>%
  group_by(site, species) %>%
  summarise(
    meanl = mean(lencm, na.rm = T),
    meanw = mean(weight, na.rm = T),
    sdl = sd(lencm, na.rm = T),
    sdw = sd(weight, na.rm = T),
    n = n(),
    sel = sdl / sqrt(n),
    sew = sdw / sqrt(n),
    ci_lowl = meanl - qt(0.975, df = n - 1) * sel,
    ci_uppl = meanl + qt(0.975, df = n - 1) * sel,
    ci_loww = meanw - qt(0.975, df = n - 1) * sew,
    ci_uppw = meanw + qt(0.975, df = n - 1) * sew,
    .groups = "drop"
  )

# Show the results so we can put in table later.
covsumlw

# Do the same as above but by species only pooling mtl individuals 
# into the 3 species/groups (PKS, BLG and HYB).
covsum <- covar %>%
  group_by(species) %>%
  summarise(
    meanl = mean(lencm, na.rm = T),
    meanw = mean(weight, na.rm = T),
    sdl = sd(lencm, na.rm = T),
    sdw = sd(weight, na.rm = T),
    n = n(),
    sel = sdl / sqrt(n),
    sew = sdw / sqrt(n),
    ci_lowl = meanl - qt(0.975, df = n - 1) * sel,
    ci_uppl = meanl + qt(0.975, df = n - 1) * sel,
    ci_loww = meanw - qt(0.975, df = n - 1) * sew,
    ci_uppw = meanw + qt(0.975, df = n - 1) * sew,
    .groups = "drop"
  )

# Show these results too to then put in table.
covsum

# Here, the produced table could be show better precision
# So we need to a) select the values we actually want to show in the table,
# and b) make the precision higher with options(pillar.sigfig =) function.
table1 <- covsum %>%
  select(species, meanl, ci_lowl, ci_uppl, meanw, ci_loww, ci_uppw)
options(pillar.sigfig = 6)

# now we can print table of mean locations
table1

################################## 2. PLOT THE DATA ####################################

# Can try to plot lengths by some variables (eg. species)
# to see what sorts of differences we see. 

sppcol <- c("PKS" = "darkorange2", "HYB" = "springgreen4","BLG" = "dodgerblue4")
ifelse(substr(metfile,8,10) == "opi", uplim <- 30, uplim <-90)

# Then we use ggplot to plot the number of samples in each location.
ggplot(covar, aes(x = site, fill = species)) + 
  geom_bar(stat = "count", position = "dodge2", width = 0.5, col = "black") +
  labs(x = "Location", y = "# of individuals") +
  scale_y_continuous(limits = c(0, uplim), breaks = seq(0, uplim, by = 10)) +
  scale_fill_manual(values = sppcol) +  # Corrected to scale_fill_manual
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

# Use filter to pick x,y range
ifelse(substr(metfile,8,10) == "opi", uplim <- 30, uplim <-70)
ifelse(substr(metfile,8,10) == "opi", lolim <- 6, lolim <- 0)

# Also want to plot out the length data to assess the distribution of sizes
# in the different species. Are some species larger on average than the others?
# Can do this using a histogram.
ggplot(covar, aes(x = lencm, fill = species)) + 
  geom_histogram(col = "black",binwidth = 2, boundary = 0) + 
  labs(x = "length (mm)", y = "Frequency") + 
  scale_y_continuous(limits = c(0, uplim), breaks = seq(0, uplim, by = 10)) +
  scale_x_continuous(limits = c(lolim, 24), breaks = seq(lolim, 24, by = 2)) +
  scale_fill_manual(values = sppcol) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3,"cm"), 
        legend.position = "top",
        legend.text = element_text(size = 20),
        legend.title = element_blank())

################################## 4. PROCRUSTES ####################################
# Perform Generalised Procrustes superimposition (GPS) of the data
# This also generates the consensus configuration (CC)
gmdat <- gpagen(tpsdat, curves = NULL, surfaces = NULL, PrinAxes = TRUE,
                Proj = TRUE, approxBE = T, print.progress = TRUE, verbose = T)

# Extract the consensus configuration (CC)
ref <- gmdat$consensus

# Use the below to define the links ourselves.
# Define landmark links based on the consensus configuration.
# Except that this is mostly good for explaoration. For our 
# purposes, we can just define the links in a matrix.
# outshape <- define.links(ref, ptsize = 2)

# Here is where we define the matrix that to use to connect the links in the CC 
# automatically (it will use a matrix to help define the landmarks that should 
# be connected). 
outshape <- cbind(c(1,2,8,11,12,13,14,15,16,17,4,3,6,9),c(2,8,11,12,13,14,15,16,17,4,1,5,7,10))

# Set plotting parameters to visualize the consensus configuration (CC)
plot_params <- list(mean.bg = "dodgerblue3",
                    mean.cex = 2,
                    pt.bg = "darkgray", 
                    pt.cex = .3,
                    link.col = "dodgerblue3",
                    txt.cex = 1.5,
                    txt.col = "springgreen3",
                    txt.pos = 3.5)


# Plot all specimen landmarks with respect to the consensus configuration.
# Help show which LM is most variable and should be checked (if needed).
plotAllSpecimens(gmdat$coords, mean = T, links = outshape, label = T,
                 plot_param = plot_params)

# Plot the concensus configuration of points for reference (e.g., in a paper)
emptycc <- gmdat$consensus
emptycc[,] <- NA

plotRefToTarget(gmdat$consensus, gmdat$consensus, method = "points", mag = 1, links = outshape,
                gridPars = gridPar(tar.pt.bg = "grey80", tar.pt.size = 3,
                                   tar.link.col = "grey80", tar.link.lwd = 3))

text(gmdat$consensus, labels = 1:nrow(gmdat$consensus), pos = 3, offset = 1, cex = 1.4, col = "red")

# Create a dataframe that combines geometric morphometric data with some of 
# the metadata
gdf <- geomorph.data.frame(gmdat,
                           id = idsall,
                           site = site,
                           species = species,
                           length = covar$lencm,
                           weight = covar$weight)

# It would be important to look at how individuals are LM in general
# and whether there may be errors in some of them. Using the outlier 
# plot can help see which ind. are very different in terms of 
# Procrustes distances
outfish <- plotOutliers(gdf$coords)

# From the figure (in plot area) - looks as though there may be 8-9.
# Set to number of outliers to plot = number of red dot in the figure
ifelse(substr(metfile,8,10) == "opi", outn <- 1, outn <- 9)

# Extract and visualise the largest outliers (depends on the figure)
outfish <- as.vector(outfish[1:outn])

# Loop that plots the major outliers against the consensus configuration
for (i in 1:length(outfish)) {
  maj_out <- mshape(gdf$coords[,,outfish[i]])
  
  # Plot the reference to target deformation grid
  plotRefToTarget(ref, maj_out, method = "points", mag = 1, links = outshape,
                  gridPars = gridPar(tar.pt.bg = "darkorange", tar.pt.size = 2,
                                     tar.link.col = "darkorange", tar.link.lwd = 2))
  # Add a title to the plot - identifying the outliers by name
  title(main = paste("Deformation plot for fish ID:", idsall[outfish[i]]), cex.main = 2)
}  

############################# 5. ALLOMERTY CHECK ##############################
# Need to check the data to see if there is an allometric relationship that 
# should be accounted for. 

# To do this we follow the details outlined in Adams et al. 2021, in the
# vignette associated with allometric analyses. Performing size correction 
# of shape data using the procD.lm cmd.

# This will create an object we can then use to test for a significant 
# relationships between shape and size.
sc_gdf <- procD.lm(coords ~ length, iter=9999, data=gdf)
sc_gdf

plot(sc_gdf, type = c("diagnostics", "regression", "PC"), 
     outliers = F,
     predictor = NULL,
     reg.type = c("PredLine"))

# unlike in traditional stats in R, the summary of a procD.lm
# model is a hypothesis test using anova. Here testing for 
# allometric relationship
summary(sc_gdf)

# It would be important to check whether the allometric relationship 
# is different for the different species or sites.

# For this we plot the allometric relationship looking for 
# large deviations in the slopes among groups. We'll try by site:
plotAllometry(sc_gdf, size = gdf$length, logsz = F, method = "PredLine",
              pch = 19, col = as.numeric(gdf$site), cex.axis = 1.5,
              cex.lab = 1.3)

# Then by species;
plotAllometry(sc_gdf, size = gdf$length, logsz = F, method = "PredLine",
              pch = 19, col = as.numeric(gdf$species), cex.axis = 1.5,
              cex.lab = 1.3)

# Here, there is significant relationship between size (ln(length)) and shape
# So, we would need produce a PCA outlining the most important 
# shape differences, while accounting for size variation.

############################# 6. PINCIPLE COMPONENT ANALYSIS ############################
# Perform a Principal components analysis on shape variables, using residuals 
# as preliminary analyses (above) show size relationship with shape.

# Here we use the residuals of the size corrected analyses above 
# 'sc-gdf$residuals'
# which are the shape variables corrected for size differences.
sc_gmpc <- gm.prcomp(sc_gdf$residuals, phy=NULL, align.to.phy=NULL, GLS=F,
                     transform=F)

# Summarise the pca results to see variance explained by new axes 
sumgmpc <- summary(sc_gmpc)

# Set a variable pv to extract the percent variance for each principal 
# component
pv <- as.numeric(sumgmpc$PC.summary[2,])

# Extract the new x and y coordinates of each fish along the new PC axes
# 1, 2, 3,... and put these into a new data frame.
sc_pcs <- as.data.frame(sc_gmpc$x)

# Also add the site and species designations to each individual fish along 
# with its log(length).

# Can also add any other data feature (here weight) we collected too.
sc_pcs$site <- gdf$site
sc_pcs$species <- gdf$species
sc_pcs$loglen <- log(gdf$length)
sc_pcs$weight <- gdf$weight

# Generate a barplot, (screeplot) demonstrating the % variance explained 
# by shape PCs.
pcs <- as.factor(to("PC", 20))
pcdf <- cbind.data.frame(PCs=pcs, pervar=pv[1:20]*100)

# This will set the values to display in our overall species PCA plots
pc1_label <- paste0("PC1 (", round(pcdf$pervar[1],1), "%)")
pc2_label <- paste0("PC2 (", round(pcdf$pervar[2],1), "%)")
pc3_label <- paste0("PC3 (", round(pcdf$pervar[3],1), "%)")

# Set a plot that looks at the percent variance explained by our new 
# most important axes of shape (i.e., PCs)
ggplot(pcdf)+
  geom_col(aes(y=pervar, x=PCs), fill = "firebrick4") + 
  labs(x = "PCs", y = "% variance explained")+
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 5)) +
  theme(aspect.ratio = 0.8)+
  theme_classic() +
  theme(axis.title.x = element_text(size = 22),
        axis.text.x = element_text(size = 18, angle = 90, hjust = 1, vjust = 0.5),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 18),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.1),
        axis.ticks.length = unit(.2,"cm"))

# Plot the transposed data (on new shape axes) with the 95% ellipses 
# (for those with enough data) with ggplot2. This will use the "sppcol"
# assigned above and the "pc1_label" and "pc2_label" from above.
ggplot(sc_pcs, aes(x = Comp1, y = Comp2, colour=species, shape = species)) + 
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
  scale_x_continuous(limits = c(-0.08, 0.08), breaks = seq(-0.08, 0.08, by = 0.02)) +
  scale_y_continuous(limits = c(-0.08, 0.08), breaks = seq(-0.08, 0.08, by = 0.02)) +
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
    legend.position = "top",
    legend.justification = "center"
  )

# Same as above but PC1 vs PC3
ggplot(sc_pcs, aes(x = Comp1, y = Comp3, colour=species, shape = species)) + 
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
  scale_x_continuous(limits = c(-0.08, 0.08), breaks = seq(-0.08, 0.08, by = 0.02)) +
  scale_y_continuous(limits = c(-0.08, 0.08), breaks = seq(-0.08, 0.08, by = 0.02)) +
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
    legend.position = "top",
    legend.justification = "center"
    )

# Same as above but PC2 vs PC3
ggplot(sc_pcs, aes(x = Comp2, y = Comp3, colour=species, shape = species)) + 
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
  scale_x_continuous(limits = c(-0.08, 0.08), breaks = seq(-0.08, 0.08, by = 0.02)) +
  scale_y_continuous(limits = c(-0.08, 0.08), breaks = seq(-0.08, 0.08, by = 0.02)) +
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
    legend.position = "top",
    legend.justification = "center"
  )

############################# 8. SHAPE BASED MANOVA ############################
# MANOVA for shape data to see if there is a significant difference in the 
# size corrected shape variables for the sampled sunfish.
shape_manova <- manova(cbind(Comp1, Comp2) ~ species, data = sc_pcs)

# Summarizing the results in a table
summary(shape_manova, test = "Wilks")
summary(shape_manova, test = "Pillai")

# Calculating the partial variance explained by the variables using eta_square
eta_squared(shape_manova, partial = T)

# Test differences along PC1 using linear model
mod_pc1 <- lm(Comp1 ~ species, data = sc_pcs)

# Make sure a linear model is appropriate by testing normality assumptions
shapiro.test(residuals(mod_pc1))

# All good with residuals - so now summarise model and run 
# anova on to test hypothesis of no differences.
summary(mod_pc1)
anova(mod_pc1)

# Can use emmeans to then summarise the pairwise comparisons
# between each species.
emmeans(mod_pc1, pairwise ~ species)

# Do the same to PC2
# Test differences along PC2
mod_pc2 <- lm(Comp2 ~ species, data = sc_pcs)

# Make sure a linear model is appropriate by testing normality assumptions
shapiro.test(residuals(mod_pc2))

# All good with residuals - so now summarise model and run 
# anova on to test hypothesis of no differences.
summary(mod_pc2)
anova(mod_pc2)

# Here there aren't any significant differences along PC2 
# emmeans can still summarise pairwise differences between each species,
# but these aren't different.
emmeans(mod_pc2, pairwise ~ species)

# So, main differences along PC1.

####################### 9. SPECIES BASED SHAPE ############################
# Next set of cmds create mean shape outline (consensus configuration) and 
# deviation of from this for each species

sp_BLG <- mshape(gmdat$coords[,,which(gdf$species=="BLG")])
plotRefToTarget(ref, sp_BLG, method = "points", mag = 4, links = outshape, 
                gridPars=gridPar(tar.pt.bg = "dodgerblue4", tar.pt.size = 2,
                                 tar.link.col = "dodgerblue4", tar.link.lwd = 2))
# Add a title to the plot
title(main = "mean Bluegill deforation", cex.main = 2, line = -5)

sp_HYB <- mshape(gmdat$coords[,,which(gdf$species=="HYB")])
plotRefToTarget(ref, sp_HYB, method = "points", mag = 4, links = outshape, 
                gridPars=gridPar(tar.pt.bg = "springgreen4", tar.pt.size = 2,
                                 tar.link.col = "springgreen4", tar.link.lwd = 2))
title(main = "mean Hybrid deformation", cex.main = 2, line = -5)


sp_PKS <- mshape(gmdat$coords[,,which(gdf$species=="PKS")])
plotRefToTarget(ref, sp_PKS, method = "points", mag = 4, links = outshape, 
                gridPars=gridPar(tar.pt.bg = "darkorange3", tar.pt.size = 2,
                                 tar.link.col = "darkorange3", tar.link.lwd = 2))
title(main = "mean Pumpkinseed deformation", cex.main = 2, line = -5)

############################# 10.ELLIPSE OVERLAP CALCULATIONS ############################
# Using the SIBER package which can estimate the overlap among different ellipses
# generated from various data sets. The results are very comparable with the CAR/sf
# packages but do not assume normal distribution or even sample sizes.

# First need to make data frame of PCs used and include species in the groups.
# could try to use community to reflect mtl vs opn. biut try this for now.
ovpc12 <- cbind.data.frame(iso1 = sc_pcs$Comp1, iso2 = sc_pcs$Comp2,
                           group = sc_pcs$species, community = rep(1,nrow(sc_pcs)))


# Second, convert data frame to SIBER object
ovpc12_so <- createSiberObject(ovpc12)


# Cmd sets the parameters for the 95% ellipses in the plotSiberObject 
# below.
group.ellipses.args  <- list(n = 1000, p.interval = 0.95, lty = 1, lwd = 2)

# Here, SIBER plots the data as well - which helps check if we get same patterns 
# as biplots (which are prettier).

#par(mfrow=c(1,1))
plotSiberObject(ovpc12_so, ax.pad = 1, hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args, 
                group.hulls = F, group.hull.args,
                bty = "L", las = 1, iso.order = c(1,2),
                xlab = pc1_label,
                ylab = pc2_label,
                x.limits = c(-0.08,0.08),
                y.limits = c(-0.08,0.08),
                cex = 0.5)

# The graph produced here should look the same. If they are different - there is a 
# problem.

# Calculate the summary statistics for each group 
# Total Area - TA
# Standard Ellipse Area - SEA
# Standard Ellipse Area for small sample sizes - SEAc

# once plotted can look at size of each ellipse and determine 
# amount of overlap
ML12 <- groupMetricsML(ovpc12_so) 
ML12

# Identify the ellipses
elp1 <-"1.BLG"
elp2 <-"1.HYB"
elp3 <-"1.PKS"

# calculate the maximum likelihood overlap of the two ellipses considered

# First, BLG - HYB
BLG_HYB <- maxLikOverlap(elp1, elp2, ovpc12_so, p.interval = 0.95, n = 1000)
BLG_HYB

# Calculating overlap between hybrids and bluegill
blghyb_ov <- BLG_HYB[3] / BLG_HYB[1]
print(paste("bluegill hybrid ov: ",blghyb_ov))

# Then, PKS - HYB
PKS_HYB <- maxLikOverlap(elp3, elp2, ovpc12_so, p.interval = 0.95, n = 1000)
PKS_HYB

# Now between hybrids and pumpkinseed
pkshyb_ov <- PKS_HYB[3] / PKS_HYB[1]
print(paste("pumpkinseed hybrid ov: ",pkshyb_ov))

# Finally, between PKS - BLG
BLG_PKS <- maxLikOverlap(elp1, elp3, ovpc12_so, p.interval = 0.95, n = 1000)
BLG_PKS

# Calculating overlap between hybrids and bluegill
blgpks_ov <- BLG_PKS[3] / BLG_PKS[1]
pksblg_ov <- BLG_PKS[3] / BLG_PKS[2]

print(paste("bluegill pumpkinseed ov: ", blgpks_ov))
print(paste("pumpkinseed blugill ov: ", pksblg_ov))

# Results show greater overlap between between hybrids and pumpkinseed, versus 
# hybrids and bluegill. It also shows that the hybrids are much more variable
# than both the bluegill and pumpkinseed, and that the hybrids are shifted 
# a bit more right along PC1 than both parental species - maybe indicating 
# some transgressive segregation in that direction. What is shape along
# PC1 (i.e., what changes teh most along PC1 in terms of shape?). Answering
# this last point would be good.

############################# 11. SHIFTS ALONG PCs ############################
# Look at the landmark (LM) loadings along PC1 (and 2) to see what features 
# change or are most variable along axes. 

# Set the landmark loadings on each PC in a data frame. 
loadings <- cbind.data.frame(LMs = as.factor(rownames(sc_gmpc$rotation)), 
                             sc_gmpc$rotation[,1:2])

# Keep the landmarks in their proper order. Otherwise, R will reorder them 
# numerically (1,10,2,20,...etc). We don't want this.
loadings$LMs <- factor(loadings$LMs, levels = unique(loadings$LMs))

# Loop to look at loadings along PC1 and PC2
for (j in 1:2) {
  
  # Create the ggplot2 bar plot
  lmlp <- ggplot(loadings, aes(x = LMs, y = loadings[,j+1])) +
    geom_bar(stat = "identity", fill = "skyblue2", alpha = 0.6, width = 1, col = "black") +
    labs(y = "Loadings", x = "Landmarks", title = paste("Loadings along PC",j)) +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.2)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 32, hjust = 0.5),
      axis.title.x = element_text(size = 22),
      axis.text.x = element_text(size = 13, angle = 90, vjust = 0.5, hjust = 1),  # Rotates labels
      axis.title.y = element_text(size = 22),
      axis.text.y = element_text(size = 15),
      axis.title = element_text(face = "bold"),
      axis.line = element_line(linewidth=1.3),
      axis.ticks = element_line(linewidth = 1.3),
      axis.ticks.length = unit(.3,"cm"),
      panel.grid.major.y = element_blank(),  # Remove horizontal grid lines
      panel.grid.minor.y = element_blank(),  # Remove minor horizontal grid lines
      panel.grid.major.x = element_line(color = "grey80"))  # Keep vertical grid lines
  print(lmlp)
}

# Shows which landmarks load most on PC1 and 2 and how they most
# likely shift.

############################# 12. ELLIPSE ECCENTRICITY ############################
# Calculating ellipse areas and eccentricities with 95% CIs

# split data by species
sppdata <- split(sc_pcs, sc_pcs$species)

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
  covm <- cov(cbind(tmp$Comp1, tmp$Comp2))
  evals <- eigen(covm, symmetric = T)$values
  evals <- sort(evals, decreasing = T) # λ1 >= λ2
  
  # semi-axes
  aaxis <- c95 * sqrt(evals[1])
  baxis <- c95 * sqrt(evals[2])
  
  # Area
  elA <- pi * aaxis * baxis
  
  # correct eccentricity
  elX <- if (evals[1] > 0) sqrt(1 - evals[2] / evals[1]) else NA
  
  # bootstrap SEs
  nboot <- 1000
  elit <- matrix(NA, nboot, 2)
  for (b in 1:nboot) {
    xboot1 <- sample(tmp$Comp1, replace = T)
    xboot2 <- sample(tmp$Comp2, replace = T)
    tcov <- cov(cbind(xboot1, xboot2))
    teval <- sort(eigen(tcov, symmetric = T)$values, decreasing = T)
    ta <- c95 * sqrt(teval[1])
    tb <- c95 * sqrt(teval[2])
    balt <- pi * ta * tb
    becc <- if (teval[1] > 0) sqrt(1 - teval[2] / teval[1]) else NA
    elit[b, ] <- c(balt, becc)
  }
  bse <- apply(elit, 2, sd, na.rm = T)
  
  # calculate 95% CI ranges
  elAci <- c(elA - 1.96 * bse[1], elA + 1.96 * bse[1])
  elXci <- c(elX - 1.96 * bse[2], elX + 1.96 * bse[2])
  
  # store results (rounded)
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