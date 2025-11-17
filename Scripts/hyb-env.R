# hyb-env.R
# This script was written to be test whether the presence of hybrids
# can be linked to different environmental variables associated to shoreline
# development. Data are taken from different sites sampled around Montreal's
# west Island and tests whether the hybrids are best explained by different 
# variables or by a null model. The differences between the models are 
# evaluated by maximum likelihood and AICc. A table is produced of the results.
# Results are also plotted but could be done better

# Written by ag and dr Fall 2025.

### Reset all instances in R
rm(list=ls())

#load packages
{
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(MuMIn)
  library(emmeans)
  library(visreg)
  }

#set working directory
setwd("/Users/denis/Library/CloudStorage/GoogleDrive-denisroy1@gmail.com/My Drive/sunfish/Scripts/MLmodel/")

#load in the data
sfdat <- read.csv("envcomp.csv", header = T, stringsAsFactors = T, na.strings = "NA")

################################## 1. NULL MODEL ############################################
#Create the null model (no effect of variables on the prob of hybrids)
mod1 <- glm(cbind(hybrids, non_hybrids) ~ 1, family = binomial(link = "logit"), data = sfdat)

#Get the summary of the model
summary(mod1)

# get the value of the coefficient on the non-transformed scale
emmeans(mod1, specs = ~1, type = "response")

#perform anova on null model NO other model is added
anova(mod1, test = "LRT") 

################################## 2. HMS2 MODEL ############################################
# Create an alt. model with HMS score explain hybrid proportions
mhms <- glm(cbind(hybrids, non_hybrids) ~ HMS2, family = binomial(link = "logit"), data = sfdat)

#Get the summary of the model
summary(mhms)

#perform anova on null model vs. mhms
anova(mod1, mhms, test = "LRT")

# get the value of the coefficient on the non-transformed scale
emmeans(mhms, specs = "HMS2", type = "response")

################################## 3. SHORELINESLOPE MODEL #################################
#Create an alt. model with shoreline slope affecting hybridization
mss <- glm(cbind(hybrids, non_hybrids) ~ shoreslope, family = binomial(link = "logit"), data = sfdat)

#Get the summary of the model
summary(mss)

#perform anova on null model vs. mss
anova(mod1, mss, test = "LRT")

# get the value of the coefficient on the non-transformed scale
emmeans(mss, specs = "shoreslope", type = "response")

################################## 4. TURBIDITY MODEL #################################
#Create an alt. model with turbidity affecting hybridization
mt <- glm(cbind(hybrids, non_hybrids) ~ turbidity, family = binomial(link = "logit"), data = sfdat)

#Get the summary of the model
summary(mt)

#perform anova on null model vs. mt
anova(mod1, mt, test = "LRT")

# get the value of the coefficient on the non-transformed scale
emmeans(mt, specs = "turbidity", type = "response")

################################## 5. CONDUCTIVITY MODEL #################################
#Create an alt. model with conductivity affecting hybridization
mc <- glm(cbind(hybrids, non_hybrids) ~ conductivity, family = binomial(link = "logit"), data = sfdat)

#Get the summary of the model
summary(mc)

#perform anova on null model vs. mc
anova(mod1, mc, test = "LRT")

# get the value of the coefficient on the non-transformed scale
emmeans(mc, specs = "conductivity", type = "response")

############################### 6. MODELS COMPARED TO NULL #################################
# Anova comparing all models to the null. We cannot use the LRT to compare the different
# models as they are not subsets of each other or nested. The best we can do here 
# is to use the LRT to test the models relative to hte null.
anova(mod1, mhms, test = "LRT")
anova(mod1, mss, test = "LRT")
anova(mod1, mt, test = "LRT")
anova(mod1, mc, test = "LRT")

############################## 7. MODEL TABLE ##############################################
# Create a list of the models tested above and define each model
models <- list(
  "Null" = mod1,
  "HMS2 (Shoreline Dev)" = mhms,
  "Shoreline Slope" = mss,
  "Conductivity" = mc,
  "Turbidity" = mt
)

# Set an empty data frame that will hold the final model table.
modelcomp <- data.frame(Model = names(models), logL = NA, K = NA, AICc = NA)

# Use a loop to compute model selection statistics
for (i in 1:length(models)) {
  model <- models[[i]]
  modelcomp$logL[i] <- as.numeric(logLik(model))
  modelcomp$K[i] <- length(coef(model))
  modelcomp$AICc[i] <- AICc(model)
}

# Calculate the AICc and AICc weights by hand
modelcomp$dAICc <- modelcomp$AICc - min(modelcomp$AICc)
modelcomp$Akaikew <- exp(-0.5 * modelcomp$dAICc) / sum(exp(-0.5 * modelcomp$dAICc))

# Likelihood Ratio Tests against null model
lrtres <- data.frame(Model = names(models)[-1], LRT = NA, pvalue = NA)

for (i in 2:length(models)) {  # Skip null model
  lrt <- anova(models[["Null"]], models[[i]], test = "Chisq")
  lrtres$LRT[i - 1] <- lrt$Deviance[2]
  lrtres$pvalue[i - 1] <- lrt$`Pr(>Chi)`[2]
}

# Merge LRT results with model comparison table
modelcom <- merge(modelcomp, lrtres, by = "Model", all.x = T)

# Sort by AICc
modelcomp <- modelcomp[order(modelcomp$AICc), ]

# Create and combine for the final table
print(modelcomp)

############################# 8. PLOTTING #############################################
# Use the correct dataset
df <- sfdat  

# Inverse logit function
inv_logit <- function(x) exp(x) / (1 + exp(x))

# Predicted probabilities and CI for the null model
df$predicted_prob_null <- predict(mod1, type = "response")
logit_ci <- confint(mod1)                  # Confidence intervals on logit scale
ci_prob <- inv_logit(logit_ci)             # Transform to probability scale
df$lower_ci <- ci_prob[1]                  # Lower bound
df$upper_ci <- ci_prob[2]                  # Upper bound

# Observed proportion of hybrids (optional for plotting)
df$observed_hybrid_prop <- df$hybrids / (df$hybrids + df$non_hybrids)

# Reorder shoreline categories
df$HMS1 <- factor(df$HMS1, levels = c("predominantly unmodified", 
                                      "obviously modified", 
                                      "significantly modified", 
                                      "severely modified"))

# Null model plot with CIs
plot_null_model_with_CI <- ggplot(df, aes(x = HMS2, y = predicted_prob_null)) +
  geom_point(aes(color = HMS1, shape = HMS1), size = 4, alpha = 0.7) +
  geom_line(aes(group = 1), color = "black", linetype = "solid", size = 0.5) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.5, alpha = 0.6, color = "black") +
  geom_text(aes(label = site), vjust = -1, hjust = 0.5, size = 3) +
  labs(
    title = "Predicted Proportion of Hybrids by Shoreline Development Score (Null Model)",
    x = "Shoreline Development Score",
    y = "Predicted Proportion of Hybrids"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("darkgreen", "#bfa500", "#ff6600", "#d62728")) +
  scale_shape_manual(values = c(16, 17, 18, 15)) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "grey55"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 11), 
    axis.title.y = element_text(size = 16)
  ) +
  scale_x_continuous(breaks = seq(0, 80, by = 10), limits = range(df$HMS2, na.rm = TRUE)) +
  scale_y_continuous(breaks = seq(0, 0.1, by = 0.05), limits = c(0, 0.1), expand = c(0, 0))

# Display the plot
print(plot_null_model_with_CI)


######################## 8b. (HMS Score plot) ######################################

# Define inverse logit function
inv_logit <- function(x) exp(x) / (1 + exp(x))

# Step 1: Get 95% CI for coefficients on logit scale
confint(mhms)  # Optional - just to see the CI of the model's intercept and slope

# Step 2: Create new data for predictions
newdata <- data.frame(HMS2 = seq(min(df$HMS2, na.rm = TRUE),
                                 max(df$HMS2, na.rm = TRUE),
                                 length.out = 100))

# Step 3: Predict on logit scale, with standard errors
preds <- predict(mhms, newdata = newdata, type = "link", se.fit = TRUE)

# Step 4: Convert to probability scale using inverse logit
newdata$fit <- inv_logit(preds$fit)
newdata$lower_ci <- inv_logit(preds$fit - 1.96 * preds$se.fit)
newdata$upper_ci <- inv_logit(preds$fit + 1.96 * preds$se.fit)

df$observed_se <- sqrt((df$hybrat * (1 - df$hybrat)) / (df$hybrids + df$non_hybrids))
z_value <- qnorm(0.975)
df$lower_ci2 <- df$hybrat - z_value * df$observed_se
df$upper_ci2 <- df$hybrat + z_value * df$observed_se

#plot the shoreline development plot
HMS2_plot <- ggplot(df, aes(x = HMS2, y = hybrat)) +
  geom_ribbon(data = newdata, aes(x = HMS2, ymin = lower_ci, ymax = upper_ci), 
              inherit.aes = FALSE, fill = "grey", alpha = 0.6) +
  geom_point(aes(shape = HMS1, color = HMS1), size = 4, alpha = 0.7) +
  geom_errorbar(aes(ymin = lower_ci2, ymax = upper_ci2), width = 0.5,alpha = 0.6, color = "black") +
  geom_line(data = newdata, aes(x = HMS2, y = fit), inherit.aes = FALSE, 
            color = "black", size = 0.5) +
  geom_text(aes(label = site), vjust = -1, hjust = 0.5, size = 3) +
  labs(
    x = "Shoreline Development Score",
    y = "Observed Proportion of Hybrids"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("darkgreen", "#bfa500", "#ff6600", "#d62728")) +
  scale_shape_manual(values = c(16, 17, 18, 15)) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "grey55"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 11), 
    axis.title.y = element_text(size = 16)
  ) +
  scale_x_continuous(breaks = seq(0, 70, by = 10), limits = c(0, 80), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(-0.025, 0.225, by = 0.05), limits = c(-0.025, 0.225))

# Display the plot
plot(HMS2_plot)

######################## 8c. Conductivity ###############################
# Define inverse logit function
inv_logit <- function(x) exp(x) / (1 + exp(x))

# Step 1: Get 95% CI for coefficients (optional)
confint(mc)

# Step 2: Create new data for predictions
newdata_cond <- data.frame(conductivity = seq(min(df$conductivity, na.rm = TRUE),
                                              max(df$conductivity, na.rm = TRUE),
                                              length.out = 100))

# Step 3: Predict on logit scale, with standard errors
preds_cond <- predict(mc, newdata = newdata_cond, type = "link", se.fit = TRUE)

# Step 4: Convert to probability scale using inverse logit
newdata_cond$fit <- inv_logit(preds_cond$fit)
newdata_cond$lower_ci <- inv_logit(preds_cond$fit - 1.96 * preds_cond$se.fit)
newdata_cond$upper_ci <- inv_logit(preds_cond$fit + 1.96 * preds_cond$se.fit)

# Step 5: Calculate 95% CI for observed proportions
df$observed_hybrid_prop <- df$hybrids / (df$hybrids + df$non_hybrids)
df$observed_se <- sqrt((df$observed_hybrid_prop * (1 - df$observed_hybrid_prop)) / 
                         (df$hybrids + df$non_hybrids))
z_value_cond <- qnorm(0.975)
df$lower_ci_cond <- df$observed_hybrid_prop - z_value_cond * df$observed_se
df$upper_ci_cond <- df$observed_hybrid_prop + z_value_cond * df$observed_se

# Step 6: Plot
conductivity_plot <- ggplot(df, aes(x = conductivity, y = observed_hybrid_prop)) +
  geom_ribbon(data = newdata_cond, aes(x = conductivity, ymin = lower_ci, ymax = upper_ci), 
              inherit.aes = FALSE, fill = "grey", alpha = 0.6) +
  geom_point(aes(shape = HMS1, color = HMS1), size = 4, alpha = 0.7) +
  geom_errorbar(aes(ymin = lower_ci_cond, ymax = upper_ci_cond), width = 0.5,alpha = 0.6, color = "black") +
  geom_line(data = newdata_cond, aes(x = conductivity, y = fit), inherit.aes = FALSE, 
            color = "black", size = 0.5) +
  geom_text(aes(label = site), vjust = -1, hjust = 0.5, size = 3) +
  labs(
    x = "Conductivity (?S/cm)",
    y = "Observed Proportion of Hybrids"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("darkgreen", "#bfa500", "#ff6600", "#d62728")) +
  scale_shape_manual(values = c(16, 17, 18, 15)) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "grey55"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 11), 
    axis.title.y = element_text(size = 16)
  ) +
  scale_x_continuous(breaks = seq(0, 450, by = 100), limits = c(0, 450), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-0.05, 0.250, by = 0.05), limits = c(-0.05, 0.250), expand = c(0,0))

# Display the plot
print(conductivity_plot)

######################## 8d. Turbidity ###############################
# Define inverse logit function
inv_logit <- function(x) exp(x) / (1 + exp(x))

# Step 1: Get 95% CI for coefficients (optional)
confint(mt)

# Step 2: Create new data for predictions
newdata_turb <- data.frame(turbidity = seq(min(df$turbidity, na.rm = TRUE),
                                           max(df$turbidity, na.rm = TRUE),
                                           length.out = 100))

# Step 3: Predict on logit scale with standard errors
preds_turb <- predict(mt, newdata = newdata_turb, type = "link", se.fit = TRUE)

# Step 4: Convert to probability scale
newdata_turb$fit <- inv_logit(preds_turb$fit)
newdata_turb$lower_ci <- inv_logit(preds_turb$fit - 1.96 * preds_turb$se.fit)
newdata_turb$upper_ci <- inv_logit(preds_turb$fit + 1.96 * preds_turb$se.fit)

# Step 5: Calculate 95% CI for observed proportions
df$observed_hybrid_prop <- df$hybrids / (df$hybrids + df$non_hybrids)
df$observed_se <- sqrt((df$observed_hybrid_prop * (1 - df$observed_hybrid_prop)) / 
                         (df$hybrids + df$non_hybrids))
z_value_turb <- qnorm(0.975)
df$lower_ci_turb <- df$observed_hybrid_prop - z_value_turb * df$observed_se
df$upper_ci_turb <- df$observed_hybrid_prop + z_value_turb * df$observed_se

# Step 6: Plot
turbidity_plot <- ggplot(df, aes(x = turbidity, y = observed_hybrid_prop)) +
  geom_ribbon(data = newdata_turb, aes(x = turbidity, ymin = lower_ci, ymax = upper_ci), 
              inherit.aes = FALSE, fill = "grey", alpha = 0.6) +
  geom_point(aes(shape = HMS1, color = HMS1), size = 4, alpha = 0.7) +
  geom_errorbar(aes(ymin = lower_ci_turb, ymax = upper_ci_turb), width = 0.5, alpha = 0.6, color = "black") +
  geom_line(data = newdata_turb, aes(x = turbidity, y = fit), inherit.aes = FALSE, 
            color = "black", size = 0.5) +
  geom_text(aes(label = site), vjust = -1, hjust = 0.5, size = 3) +
  labs(
    x = "Turbidity (FNU)",
    y = "Observed Proportion of Hybrids"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("darkgreen", "#bfa500", "#ff6600", "#d62728")) +
  scale_shape_manual(values = c(16, 17, 18, 15)) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "grey55"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 11), 
    axis.title.y = element_text(size = 16)
  ) +
  scale_x_continuous(breaks = seq(0, 12, by = 2), limits = c(0, 12), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-0.025, 0.275, by = 0.05), limits = c(-0.025, 0.275))

# Display the plot
print(turbidity_plot)

######################## 8d. shoreslope ###############################
# Define inverse logit function
inv_logit <- function(x) exp(x) / (1 + exp(x))

# Step 1: Get 95% CI for coefficients (optional)
confint(mss)

# Step 2: Create new data for predictions
newdata_slope <- data.frame(shoreslope = seq(min(df$shoreslope, na.rm = TRUE),
                                        max(df$shoreslope, na.rm = TRUE),
                                        length.out = 100))

# Step 3: Predict on logit scale with standard errors
preds_slope <- predict(mss, newdata = newdata_slope, type = "link", se.fit = TRUE)

# Step 4: Convert to probability scale
newdata_slope$fit <- inv_logit(preds_slope$fit)
newdata_slope$lower_ci <- inv_logit(preds_slope$fit - 1.96 * preds_slope$se.fit)
newdata_slope$upper_ci <- inv_logit(preds_slope$fit + 1.96 * preds_slope$se.fit)

# Step 5: Calculate 95% CI for observed proportions
df$observed_hybrid_prop <- df$hybrids / (df$hybrids + df$non_hybrids)
df$observed_se <- sqrt((df$observed_hybrid_prop * (1 - df$observed_hybrid_prop)) / 
                         (df$hybrids + df$non_hybrids))
z_value_slope <- qnorm(0.975)
df$lower_ci_slope <- df$observed_hybrid_prop - z_value_slope * df$observed_se
df$upper_ci_slope <- df$observed_hybrid_prop + z_value_slope * df$observed_se

# Step 6: Plot
slope_plot <- ggplot(df, aes(x = shoreslope, y = observed_hybrid_prop)) +
  geom_ribbon(data = newdata_slope, aes(x = shoreslope, ymin = lower_ci, ymax = upper_ci), 
              inherit.aes = FALSE, fill = "grey", alpha = 0.6) +
  geom_point(aes(shape = HMS1, color = HMS1), size = 4, alpha = 0.7) +
  geom_errorbar(aes(ymin = lower_ci_slope, ymax = upper_ci_slope), width = 0.4, alpha = 0.6, color = "black") +
  geom_line(data = newdata_slope, aes(x = shoreslope, y = fit), inherit.aes = FALSE, 
            color = "black", size = 0.5) +
  geom_text(aes(label = site), vjust = -1, hjust = 0.5, size = 3) +
  labs(
    x = "Shoreline Slope (Degrees)",
    y = "Observed Proportion of Hybrids"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("darkgreen", "#bfa500", "#ff6600", "#d62728")) +
  scale_shape_manual(values = c(16, 17, 18, 15)) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "grey55"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 11), 
    axis.title.y = element_text(size = 16)
  ) +
  scale_x_continuous(breaks = seq(0, 20, by = 4), limits = c(0, 20), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-0.025, 0.225, by = 0.05), limits = c(-0.025, 0.225), expand = c(0,0))

  # Display the plot
print(slope_plot)


############################ 9.COMBINE PLOTS #############################################
# Remove y-axis labels from all but the conductivity plot
plot_null_model_with_CI <- plot_null_model_with_CI + theme(axis.title.y = element_blank())
HMS2_plot <- HMS2_plot + theme(axis.title.y = element_blank())
slope_plot <- slope_plot + theme(axis.title.y = element_blank())
turbidity_plot <- turbidity_plot + theme(axis.title.y = element_blank())


# Combine plots with turbidity
plot_null_model_with_CI / HMS2_plot / conductivity_plot / slope_plot / turbidity_plot + plot_layout(guides = "collect") & theme(legend.position = "bottom")

# Combine plots without the turbidity
plot_null_model_with_CI / HMS2_plot / conductivity_plot / slope_plot + plot_layout(guides = "collect") & theme(legend.position = "bottom")

##### END #####

