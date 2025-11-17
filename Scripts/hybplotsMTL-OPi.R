# This script is to make the plot for the assignment probabilities of
# Hybrid and Parentals sunfish in Montreal and Opinicon

# It takes in data that were generated from the Geneclass/Hybridlab analyses
# and makes stacked horizontal bar charts out of them. This hels visualise 
# assignment probabilities of individuals to derived classes of parental 
# fish, and relatively simple hybrid classes.

# Written by AG and DR Fall 2025

# reset all instance in r
rm(list=ls())

# load packages
library(ggplot2)
library(dplyr)
 
################## Montreal Dataset #################
#Create the dataset
data <- data.frame(
  Category = rep(c("PKS", "BLG ass. to PKS", "BLG", "PKS ass. to BLG"), each = 3),
  Group = factor(rep(c("Excluded (P??? 0.05)", "Assigned (0.05 < P < 0.90)", "Assigned (P ??? 0.90)"), times = 4), 
                 levels = c("Assigned (P ??? 0.90)", "Assigned (0.05 < P < 0.90)", "Excluded (P??? 0.05)")),
  Probability = c(
    0.010204082, 0.836734694, 0.153061224,
    0.959183673, 0.040816327, 0,
    0.015384615, 0.692307692, 0.292307692,
    0.953846154, 0.030769231, 0.015384615
  )
)

# Define category order, must be in this order for plot to make sense
category_order <- c("PKS", "BLG ass. to PKS", "BLG", "PKS ass. to BLG")

# Define custom colors
custom_colors <- c("Assigned (P ??? 0.90)" = "black", "Assigned (0.05 < P < 0.90)" = "white", "Excluded (P??? 0.05)" = "darkgrey")

# legend labels
legend_labels <- c(
  "Assigned (P ??? 0.90)" = "Assigned (P \u2265 0.90)", 
  "Assigned (0.05 < P < 0.90)" = "Assigned (0.05 < P < 0.90)", 
  "Excluded (P ??? 0.05)" = "Excluded (P \u2264 0.05)"
)

# Plot
p <- ggplot(data, aes(x = factor(Category, levels = rev(category_order)), y = Probability, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = custom_colors) +
  coord_flip() +
  labs(x = "Category", y = "Probability", fill = "Group") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Print the plot
print(p)


# FIX aesthetics
p <- ggplot(data, aes(x = factor(Category, levels = rev(category_order)), y = Probability, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = 0.5) +
  scale_fill_manual(values = custom_colors) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(scale = 1), breaks = seq(0, 1, by = 0.25)) +  # Ensure percentage format on y-axis
  labs(x = "Potential source population", y = "Proportion of individuals", fill = NULL) +  # Remove "Group" from legend
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 14),  # Increase legend text size
    axis.title = element_text(size = 14),  # Increase axis title size
    axis.text = element_text(color = "black", size = 12)  # Make the ticks darker and larger
  ) +
  guides(fill = guide_legend(
    title = NULL,  # Remove legend title
    override.aes = list(label = c(
      "Excluded (P ??? 0.05)", 
      "Assigned (0.05 < P < 0.90)", 
      "Assigned (P ??? 0.90)"  # Use proper Unicode for "???"
    ))
  ))

# Print the plot
print(p)

############################ Opinicon Dataset ########################
#Create the dataset
Opidata <- data.frame(
  Category = rep(c("PKS", "BLG ass. to PKS", "BLG", "PKS ass. to BLG"), each = 3),
  Group = factor(rep(c("Excluded (P??? 0.05)", "Assigned (0.05 < P < 0.90)", "Assigned (P ??? 0.90)"), times = 4), 
                 levels = c("Assigned (P ??? 0.90)", "Assigned (0.05 < P < 0.90)", "Excluded (P??? 0.05)")),
  Probability = c(
    0, 1.0, 0,
    1.0, 0, 0,
    0, 0.86486, 0.13514,
    1.0, 0, 0
  )
)

# Define category order, must be in this order
OPI_category_order <- c("PKS", "BLG ass. to PKS", "BLG", "PKS ass. to BLG")

# Define custom colors
OPI_custom_colors <- c("Assigned (P ??? 0.90)" = "black", "Assigned (0.05 < P < 0.90)" = "white", "Excluded (P??? 0.05)" = "darkgrey")

# legend labels
OPI_legend_labels <- c(
  "Assigned (P ??? 0.90)" = "Assigned (P \u2265 0.90)", 
  "Assigned (0.05 < P < 0.90)" = "Assigned (0.05 < P < 0.90)", 
  "Excluded (P ??? 0.05)" = "Excluded (P \u2264 0.05)"
)

# Plot
Opiplot <- ggplot(Opidata, aes(x = factor(Category, levels = rev(OPI_category_order)), y = Probability, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = OPI_custom_colors) +
  coord_flip() +
  labs(x = "Category", y = "Probability", fill = "Group") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Print the plot
print(Opiplot)

# FIX aesthetics
Opiplot <- ggplot(Opidata, aes(x = factor(Category, levels = rev(OPI_category_order)), y = Probability, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = 0.5) +
  scale_fill_manual(values = OPI_custom_colors) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(scale = 1), breaks = seq(0, 1, by = 0.25)) +  # Ensure percentage format on y-axis
  labs(x = "Potential source population", y = "Proportion of individuals", fill = NULL) +  # Remove "Group" from legend
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 14),  # Increase legend text size
    axis.title = element_text(size = 14),  # Increase axis title size
    axis.text = element_text(color = "black", size = 12)  # Make the ticks darker and larger
  ) +
  guides(fill = guide_legend(
    title = NULL,  # Remove legend title
    override.aes = list(label = c(
      "Excluded (P ??? 0.05)", 
      "Assigned (0.05 < P < 0.90)", 
      "Assigned (P ??? 0.90)"  # Use proper Unicode for "???"
    ))
  ))

# Print the plot
print(Opiplot)

######################## Combine/display plots ###############################
library(patchwork)

# Combine plots & keep only one legend
p + Opiplot + plot_layout(guides = "collect") & theme(legend.position = "bottom")


######################## Hybrids in Montreal ###########################################
# Create the dataset
mtlhybdata <- data.frame(
  Category = rep(c("Bluegill", "Pumpkinseed", "F1", "F1XBLG", "F1XPKS", "F2", "F2XBLG", "F2XPKS"), each = 3),
  Group = factor(rep(c("Excluded (P??? 0.05)", "Assigned (0.05 < P < 0.90)", "Assigned (P ??? 0.90)"), times = 8), 
                 levels = c("Assigned (P ??? 0.90)", "Assigned (0.05 < P < 0.90)", "Excluded (P??? 0.05)")),
  Probability = c(
    0.368421053, 0.631578947, 0,
    0.526315789, 0.473684211, 0,
    0.052631579, 0.789473684, 0.157894737,
    0.157894737, 0.789473684, 0.052631579,
    0.052631579, 0.842105263, 0.105263158,
    0.052631579, 0.894736842, 0.052631579,
    0.210526316, 0.789473684, 0,
    0.052631579, 0.894736842, 0.052631579
  )
)

# Define category order
hybcategory_order <- c("Bluegill", "Pumpkinseed", "F1", "F1XBLG", "F1XPKS", "F2", "F2XBLG", "F2XPKS")

# Define custom colors
hybcustom_colors <- c("Assigned (P ??? 0.90)" = "black", "Assigned (0.05 < P < 0.90)" = "white", "Excluded (P??? 0.05)" = "darkgrey")

# Legend labels
hyblegend_labels <- c(
  "Assigned (P ??? 0.90)" = "Assigned (P \u2265 0.90)", 
  "Assigned (0.05 < P < 0.90)" = "Assigned (0.05 < P < 0.90)", 
  "Excluded (P ??? 0.05)" = "Excluded (P \u2264 0.05)"
)

# Load necessary package
library(ggplot2)

# Plot
hybmontp <- ggplot(mtlhybdata, aes(x = factor(Category, levels = rev(hybcategory_order)), y = Probability, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = hybcustom_colors) +
  coord_flip() +
  labs(x = "Category", y = "Probability", fill = "Group") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Print the plot
print(hybmontp)

# Fix aesthetics
hybmontp <- ggplot(mtlhybdata, aes(x = factor(Category, levels = rev(hybcategory_order)), y = Probability, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = 0.5) +
  scale_fill_manual(values = hybcustom_colors) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(scale = 1), breaks = seq(0, 1, by = 0.25)) +
  labs(x = "Potential source population", y = "Proportion of individuals", fill = NULL) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(color = "black", size = 12)
  ) +
  guides(fill = guide_legend(
    title = NULL,
    override.aes = list(label = c(
      "Excluded (P ??? 0.05)", 
      "Assigned (0.05 < P < 0.90)", 
      "Assigned (P ??? 0.90)"
    ))
  ))

# Print the refined plot
print(hybmontp)


######################## Hybrids in Opinicon ###########################################
# Create the dataset
opihybdata <- data.frame(
  Category = rep(c("Bluegill", "Pumpkinseed", "F1", "F1XBLG", "F1XPKS", "F2", "F2XBLG", "F2XPKS"), each = 3),
  Group = factor(rep(c("Excluded (P??? 0.05)", "Assigned (0.05 < P < 0.90)", "Assigned (P ??? 0.90)"), times = 8), 
                 levels = c("Assigned (P ??? 0.90)", "Assigned (0.05 < P < 0.90)", "Excluded (P??? 0.05)")),
  Probability = c(
    0.85, 0.15, 0,
    1.00, 0, 0,
    0.60, 0.40, 0,
    0.60, 0.40, 0,
    0.75, 0.25, 0,
    0.65, 0.35, 0,
    0.65, 0.35, 0,
    0.75, 0.25, 0
  ) # Removed the trailing comma
)


# Define category order
ohybcategory_order <- c("Bluegill", "Pumpkinseed", "F1", "F1XBLG", "F1XPKS", "F2", "F2XBLG", "F2XPKS")

# Define custom colors
ohybcustom_colors <- c("Assigned (P ??? 0.90)" = "black", "Assigned (0.05 < P < 0.90)" = "white", "Excluded (P??? 0.05)" = "darkgrey")

# Legend labels
ohyblegend_labels <- c(
  "Assigned (P ??? 0.90)" = "Assigned (P \u2265 0.90)", 
  "Assigned (0.05 < P < 0.90)" = "Assigned (0.05 < P < 0.90)", 
  "Excluded (P ??? 0.05)" = "Excluded (P \u2264 0.05)"
)

# Load necessary package
library(ggplot2)

# Plot
hybopip <- ggplot(opihybdata, aes(x = factor(Category, levels = rev(ohybcategory_order)), y = Probability, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = ohybcustom_colors) +
  coord_flip() +
  labs(x = "Category", y = "Probability", fill = "Group") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Print the plot
print(hybopip)

# Fix aesthetics
hybopip <- ggplot(opihybdata, aes(x = factor(Category, levels = rev(ohybcategory_order)), y = Probability, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = 0.5) +
  scale_fill_manual(values = ohybcustom_colors) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(scale = 1), breaks = seq(0, 1, by = 0.25)) +
  labs(x = "Potential source population", y = "Proportion of individuals", fill = NULL) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(color = "black", size = 12)
  ) +
  guides(fill = guide_legend(
    title = NULL,
    override.aes = list(label = c(
      "Excluded (P ??? 0.05)", 
      "Assigned (0.05 < P < 0.90)", 
      "Assigned (P ??? 0.90)"
    ))
  ))

# Print the refined plot
print(hybopip)

################ Combine Plots ####################

library(patchwork)

# Combine plots & keep only one legend
p + Opiplot + hybmontp +hybopip + plot_layout(guides = "collect") & theme(legend.position = "bottom")

