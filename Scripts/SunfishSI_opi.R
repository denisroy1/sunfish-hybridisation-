# SunfishSI.R
# Script for stable isotope analyses of sunfish and baseline data.
# Data are read with both the fish and baseline data (baseline at 
# the end). Script splits data into fish and baseline, to allow easier 
# manipulations. Baseline data corrected for lipid content as suggested 
# by Post 2007, but not applied to fish (they are not high lipid species).
# We use C signal to estimate alpha (from where the fish get most of their 
# carbon) value which we then use to correct the N signal to convert
# into trophic position (tp).

# Corrected data are plotted along with the mean and 95%CI for the 
# tp and dC of the fish and the littoral("SNA"), and pelagic("MUS") baselines.
# This let's us see where the different species are likely focussing 
# their feeding. 

# The script then tests differences along the C and N axes and also estimates 
# the ML overlap among the species. Finally, it estimates the ellipse areas,
# and their eccentricities as performed in the other scripts.

# written by DR and AG Oct 2025

################# 1. HOUSEKEEPING AND DATA READ IN ############
# Clearing all instances
rm(list=ls())

#Loading needed libraries
{
  library(parallel)
  library(tidyverse)
  library(ggplot2)
  library(plyr)
  library(mvtnorm)
  library(emmeans)
  library(car)
  library(dplyr)
  library(FSA)
  library(SIBER)
  library(rjags)
  library(lsr)
  #library(conover.test)
}

# Setting the core allocation to run parallel computations
options(mc.cores = parallel::detectCores())
cores <- detectCores(logical = T)-1

#Set working directory
setwd('/Users/denis/Library/CloudStorage/GoogleDrive-denisroy1@gmail.com/My Drive/sunfish/Scripts/SI/')

# load the data
sidat <- read.csv("sunfishSI_opi.csv", header = T, stringsAsFactors = T)

# Look to see the data's been read in properly
head(sidat)

# Extract names used for SI (these ought to match TPS file IDs)
ids <- sidat$NEW_TPSID

# Extract and assign species assignment for each individual using the IDs
species <- as.factor(substr(ids,7,9))

# Extract and assign site factor to each sample using the IDs
sidat$site <- as.factor(substr(ids,1,3))

################# 2. SPLITTING DATA FISH/BASELINE ############
# Filter only fish data using the species entries BLG, HYB, and PKS
fsidat <- sidat %>%
  filter(species %in% c("BLG", "HYB", "PKS"))

# Drop unused species levels
fsidat <- droplevels(fsidat)
table(fsidat$species)

# Filter just the baseline data using site information and 
# species. Here use MUS for pelagic and SNA for littoral
# baselines
bsidat <- sidat %>%
  filter(species %in% c("MUS","SNA"))
bsidat <- droplevels(bsidat)
table(bsidat$species)

# Finally, should consider mathematically lipid correcting
# the carbon as per outlined in Post 2007
bsidat$dClc <- (bsidat$dC - 3.32) + (0.99 * bsidat$C.N)

################# 3. PRELIMINARY PLOT #######################
# Get the ranges of the C and N isotope ratio values overall individual
# for plotting later
xr <- range(sidat$dC)
yr <- range(sidat$dN)

# Set colours to use and should be similar to those used in other scripts
mycol <- c('PKS' = 'darkorange2', 
           'BLG' = 'dodgerblue4', 
           'HYB' = 'springgreen4', 
           'SNA' = 'purple', 
           'MUS' = 'red')

# Plot the filtered data with the 95% ellipses
ggplot(fsidat, aes(x = dC, y = dN, colour = species, shape = species)) + 
  geom_point(size = 5, stroke = 1, alpha = 0.7) +
  stat_ellipse(aes(fill = species), geom = "polygon", color = NA, alpha = 0.2, 
               level = 0.95, type = "norm", show.legend = F) +
  stat_ellipse(aes(colour = species), geom = "path", linewidth = 1, 
               level = 0.95, type = "norm") +
  scale_fill_manual(values = mycol) +
  scale_color_manual(values = mycol) +
  labs(x = expression({delta}^13*C~'\u2030'), y = expression({delta}^15*N~'\u2030')) +
  scale_x_continuous(limits = c(-32, -14), breaks = seq(-32, -14, by = 2)) +
  scale_y_continuous(limits = c(5, 16),  breaks = seq(5, 16, by = 2)) +
  theme(aspect.ratio = 0.70) +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "lightgrey"),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(linewidth = 1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3, "cm")) +
  theme(legend.text = element_text(size = 17)) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.justification = "center") +
  guides(fill = "none")   # removes ellipse fill legend

################# 4. BASELINE CORRECTIONS #######################
# The sunfish from Opinicon were all collected in around the same area 
# near QUBS. The sites were very close to one another. So, one baseline 
# correction is good for all the Opinicon fish.

# Here we also don't need a function as we can process all fish at 
# once without having to filter through sites. 

Csidat <- fsidat

lam    <- 2       # baseline trophic level for primary consumers
dCfrac <- 0.5     # d13C diet tissue fractionation for C
dNfrac <- 3.4     # d15N trophic enrichment factor

# Extract the baselines for opi 
bll <- bsidat[bsidat$species == "SNA", c("dClc", "dN")]
blp <- bsidat[bsidat$species == "MUS", c("dClc", "dN")]

# Fraction of the diet that's littoral
all <- (Csidat$dC - dCfrac - blp$dClc) / (bll$dClc - blp$dClc)
all <- pmax(0, pmin(1, all))  
Csidat$all <- all
    
# Fraction of the diet that's pelagic
alp <- (Csidat$dC - dCfrac - bll$dClc) / (blp$dClc - bll$dClc)
alp <- pmax(0, pmin(1, alp))  
Csidat$alp <- alp
    
# Estimating trophic position (tp)
mixdN <- alp * bll$dN + (1 - alp) * blp$dN
tp <- lam + (Csidat$dN - mixdN) / dNfrac
Csidat$tp <- tp
    
# Correcting d13C for littoral source
mixdC <- alp * bll$dClc + (1 - alp) * blp$dClc + dCfrac
dCc <- Csidat$dC - mixdC
Csidat$dCc <- dCc
    
# Formula from above were taken from Post 2002
    

# Here since we only have one value we plot these as reference 
# for littoral vs. pelagic carbon sources 
basestat <- cbind.data.frame(species = bsidat$species, y = c(1.5, 2.0),
                               dCm = bsidat$dClc)


# Plot with baseline points for each group and baseline means + 
# error bars to the plot
ggplot(Csidat, aes(x = dC, y = tp, color = as.factor(species), 
                   fill = as.factor(species), shape = as.factor(species))) + 
  geom_point(size = 5, stroke = 1, alpha = 0.7) +
  stat_ellipse(aes(fill = as.factor(species)), geom = "polygon", 
               color = NA, alpha = 0.3, level = 0.95, type = "norm") +
  stat_ellipse(aes(color = as.factor(species)), geom = "path", 
               linewidth = 1, level = 0.95, type = "norm") +
  # Baseline
  geom_point(data = basestat, inherit.aes = F,
             aes(x = dCm, y = y, color = as.factor(species), 
                 shape = as.factor(species)), 
             size = 5, stroke = 1) +

  scale_fill_manual(values = mycol) +
  scale_color_manual(values = mycol) +
  labs(x = expression({delta}^13*C~'\u2030'), y = "Trophic position") +
  scale_x_continuous(limits = c(-38, -14), breaks = seq(-38, -14, by = 2)) +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 0.5)) +
  theme(aspect.ratio = 0.70) +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "lightgrey"),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(linewidth = 1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3, "cm")) +
  theme(legend.text = element_text(size = 17)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        legend.justification = "center") +
  guides(fill = "none")

################# 5. TESTING FOR dC DIFFERENCES #######################
# Testing the hypothesis of no differences 
# Linear model for Carbon Isotope (Csidat)
lmC <- lm(dC ~ species, data = Csidat)

# Testing normality of residuals for Carbon
shapiro.test(residuals(lmC))

# The data do not seem to be normally distributed and there aren't many 
# easy transformations in a glm that would be applicable
# So, run a KS test for overall differences
ksdC <- kruskal.test(dC ~ species, data = Csidat)
ksdC

# Need to test the pairwise differences noting that the "PKS" only 
# has 3 individuals and so not likely very robust. Also correcting 
# the tests for multiple comparisons
pwksdC <- dunnTest(dC ~ species, method = "by", two.sided = T, data = Csidat)
pwksdC

# Using emmeans to get the central tendencies of the dC data
dCstats <- emmeans(lmC, specs = "species")
dCstats

# Use the above to create a df to plot the differences
# in a stripchart
sumdC <- as.data.frame(dCstats)

# Visualise the results of the model
ggplot(Csidat, aes(x = species, y = dC, colour = species)) +
  geom_jitter(size = 5, alpha = 0.7, stroke = 1) +
  geom_point(data = sumdC, aes(x = species, y = emmean), colour = "black", 
             fill = "black", alpha = 0.7, shape = 21, size = 6, inherit.aes = F) +
  geom_errorbar(data = sumdC, aes(x = species, ymin = lower.CL, ymax = upper.CL),
                colour = "black", width = 0.15, linewidth = 1.1, inherit.aes = F) +
  scale_y_continuous(limits = c(-32, -14), breaks = seq(-32, -14, by = 2)) +
  scale_color_manual(values = mycol) +
  labs(x = "Species", y = expression({delta}^13*C~'\u2030')) +
  annotate("text", x = 1, y = -14, label = paste("X"), color = "black", size = 10, family = "Courier") +
  annotate("text", x = 2, y = -14, label = paste("Y"), color = "black", size = 10, family = "Courier") +
  annotate("text", x = 3, y = -14, label = paste("XY"), color = "black", size = 10, family = "Courier") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(linewidth = 1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3, "cm")) +
  theme(legend.text = element_text(size = 17)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        legend.justification = "center") +
  guides(fill = "none")

################# 6. TESTING FOR dN DIFFERENCES #######################
# Testing the hypothesis of no differences 
# Linear model for trophic position differences (Csidat)
# However _ important to note that differences are not reals and 
# anything within 3.4 tps from one another are typically considered 
# the same. We can, however, call these differences in trophic signals
# with the trophic level.

# First try linear model
lmN <- lm(tp ~ species, data = Csidat)

# Testing normality of residuals for Carbon
shapiro.test(residuals(lmN))

# Visualise the residuals' distribution
hist(residuals(lmN))

# Summarise the model and see coefficients
summary(lmN)

# Use anova to test hypothesis for Carbon Isotope (dC)
aN <- anova(lmN)
aN

# Pairwise comparisons for Carbon with adjustment for multiple tests
dNstats <- emmeans(lmN, pairwise ~ species)

# Show results
dNstats

# Use the above to also create a df to plot the differences
# in a stripchart
sumdN <- as.data.frame(dNstats$emmeans)

# Visualise the results of the model
ggplot(Csidat, aes(x = species, y = tp, colour = species)) +
  geom_jitter(size = 5, alpha = 0.7, stroke = 1) +
  geom_point(data = sumdN, aes(x = species, y = emmean), colour = "black", 
             fill = "black", alpha = 0.7, shape = 21, size = 6, inherit.aes = F) +
  geom_errorbar(data = sumdN, aes(x = species, ymin = lower.CL, ymax = upper.CL),
                colour = "black", width = 0.15, linewidth = 1.1, inherit.aes = F) +
  scale_y_continuous(limits = c(2.5, 4), breaks = seq(2.5, 4, by = 0.5)) +
  scale_color_manual(values = mycol) +
  labs(x = "Species", y = expression({delta}^15*N~'\u2030')) +
  annotate("text", x = 1, y = 4, label = paste("X"), color = "black", size = 10, family = "Courier") +
  annotate("text", x = 2, y = 4, label = paste("XY"), color = "black", size = 10, family = "Courier") +
  annotate("text", x = 3, y = 4, label = paste("Y"), color = "black", size = 10, family = "Courier") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(linewidth = 1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3, "cm")) +
  theme(legend.text = element_text(size = 17)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        legend.justification = "center") +
  guides(fill = "none")

############################# 7.ELLIPSE OVERLAP CALCULATIONS ############################
# Using the SIBER package which can estimate the overlap among different ellipses
# generated from si data. 

# Make data frame of Csidat and include species in the groups.
# could try to use community to reflect mtl vs opn. but try this for now.
ovcsi <- cbind.data.frame(iso1 = Csidat$dC, iso2 = Csidat$tp,
                             group = Csidat$species, community = rep(1,nrow(Csidat)))

# Second, convert data frame to SIBER object
ovcsi_so <- createSiberObject(ovcsi)

# Cmd sets the parameters for the 95% ellipses in the plotSiberObject 
# below.
group.ellipses.args  <- list(n = 1000, p.interval = 0.95, lty = 1, lwd = 2)

# Here, SIBER plots the data as well - which helps check if we get same patterns 
# as biplots (which are prettier).

#par(mfrow=c(1,1))
plotSiberObject(ovcsi_so, ax.pad = 1, hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args, 
                group.hulls = F, group.hull.args,
                bty = "L", las = 1, iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = "Trophic position",
                x.limits = c(-32,-14),
                y.limits = c(2,5),
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
MLCN <- groupMetricsML(ovcsi_so) 
MLCN

# Identify the ellipses
elp1 <-"1.BLG"
elp2 <-"1.HYB"
elp3 <-"1.PKS"

# calculate the maximum likelihood overlap of the two ellipses considered

# First, BLG - HYB
BLG_HYB <- maxLikOverlap(elp1, elp2, ovcsi_so, p.interval = 0.95, n = 1000)
BLG_HYB

# Calculating overlap between hybrids and bluegill
blghyb_ov <- BLG_HYB[3] / BLG_HYB[1]
print(paste("bluegill hybrid ov: ",blghyb_ov))

# Then, PKS - HYB
PKS_HYB <- maxLikOverlap(elp3, elp2, ovcsi_so, p.interval = 0.95, n = 1000)
PKS_HYB

# Now between hybrids and pumpkinseed
pkshyb_ov <- PKS_HYB[3] / PKS_HYB[1]
print(paste("pumpkinseed hybrid ov: ",pkshyb_ov))

# Finally, between PKS - BLG
BLG_PKS <- maxLikOverlap(elp1, elp3, ovcsi_so, p.interval = 0.95, n = 1000)
BLG_PKS

# Calculating overlap between hybrids and bluegill
blgpks_ov <- BLG_PKS[3] / BLG_PKS[1]
pksblg_ov <- BLG_PKS[3] / BLG_PKS[2]

print(paste("bluegill pumpkinseed ov: ", blgpks_ov))
print(paste("pumpkinseed blugill ov: ", pksblg_ov))

# "opi": Results show some overlap in all 3 species here too. However, the
# BLG (in black) this time is more than the HYB, and the HYB tends to be 
# more littoral than the BLG. PKS is not super reliable as we only have 3 ind.
# for this species. Some BLG seem more pelagic, but both HYB and BLG 
# seem centered on SNA vs MUS in relation to baseline 

# Test results along the dC and tp axes tend to confirm that HYB and BLG 
# are different and that both are not different from PKS which has ind. 
# in both other species SI space.

############################# 8. ELLIPSE ECCENTRICITY ############################
# Calculating ellipse areas and eccentricities with 95% CIs

# split data by species
sppdata <- split(Csidat, Csidat$species)

# prepare empty data frame for results
elstat <- data.frame(
  species = character(length(sppdata)),
  elA = numeric(length(sppdata)),
  elAlwr = numeric(length(sppdata)),
  elAupr = numeric(length(sppdata)),
  elX = numeric(length(sppdata)),
  elXlwr = numeric(length(sppdata)),
  elXupr = numeric(length(sppdata)),
  stringsAsFactors = F
)

# scaling factor for 95% CI ellipse
c95 <- sqrt(qchisq(0.95, 2))

# Looping through the spp
for (a in seq_along(sppdata)) {
  tmp <- sppdata[[a]]
  
  # covariance matrix of dC and tp
  covm <- cov(cbind(tmp$dC, tmp$tp))
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
    xboot1 <- sample(tmp$dC, replace = T)
    xboot2 <- sample(tmp$tp, replace = T)
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

