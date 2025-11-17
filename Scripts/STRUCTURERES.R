#K=2 Evanno plots script

#reset all instances in R
rm(list=ls())

library(dplyr)
library(ggplot2)

# --- Step 1: Read CSV ---
df <- read.csv("Evannoplotdata.csv")

# --- Step 2: Calculate Delta K ---
diff1 <- diff(df$mean_LnP)
diff2 <- diff(diff1)
deltaK <- abs(diff2) / df$sd_LnP[2:(nrow(df)-1)]

plot_df <- data.frame(
  K = 2:(nrow(df)-1),
  deltaK = deltaK
)
 plot_df
 
#######Evanno plot K=1-5###############################
plot_1to5 <- ggplot(plot_df %>% filter(K <= 5), aes(x = K, y = deltaK)) +
  geom_point(size = 5, color = "blue", alpha = 0.8) +
  geom_line(color = "blue", linewidth = 1.2) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "red", linewidth = 1) +
  scale_x_continuous(limits = c(1, 5), breaks = 1:5) +
  scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, by = 100)) +
  labs(
    x = "K",
    y = expression(Delta * K)
  ) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "lightgrey"),
    axis.text.x = element_text(size = 22),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(linewidth = 1.3),
    axis.ticks = element_line(linewidth = 1.3),
    axis.ticks.length = unit(.3, "cm")
  )

#print
plot_1to5

#export the plot at a certain dpi 
ggsave("k1-5_500dpi.jpg",
       plot = plot_1to5,   # specify the plot object
       width = 8,           # width in inches
       height = 6,          # height in inches
       units = "in",
       dpi = 500
)

####################now lnp########################################
lnp_plot <- ggplot(df %>% filter(K <= 5), aes(x = K, y = mean_LnP)) +
  geom_point(size = 5, color = "blue", alpha = 0.8) +
  geom_line(color = "blue", linewidth = 1.2) +
  geom_errorbar(
    aes(ymin = mean_LnP - sd_LnP, ymax = mean_LnP + sd_LnP),
    width = 0.2,
    color = "forestgreen",
    linewidth = 1
  ) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "red", linewidth = 1) +
  scale_x_continuous(limits = c(1, 5), breaks = 1:5) +
  scale_y_continuous(limits = c(-8000, -6500), breaks = seq(-8000, -6500, by = 200)) +
  labs(
    x = "K",
    y = "Mean LnP(K) ± SD"
  ) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "lightgrey"),
    axis.text.x = element_text(size = 22),  # horizontal labels
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(linewidth = 1.3),
    axis.ticks = element_line(linewidth = 1.3),
    axis.ticks.length = unit(.3, "cm")
  )

#print
lnp_plot

#export the plot at a certain dpi 
ggsave("lnpplot_500dpi.jpg",
       plot = lnp_plot,   # specify the plot object
       width = 8,           # width in inches
       height = 6,          # height in inches
       units = "in",
       dpi = 500
)


#Now let's add the structure plot

library(tidyr)

# Step 1: Read CSV
data <- read.csv("sunfish_indivq.csv")

# Step 2: Ensure ancestry columns are numeric
data$blue <- as.numeric(data$blue)
data$orange <- as.numeric(data$orange)

# Step 3: Create a combined group: Region + Status
data <- data %>%
  mutate(Group = paste(Region, status, sep = "_")) %>%
  arrange(Region, status) %>%
  mutate(IndIndex = 1:nrow(data))

# Step 4: Convert to long format
data_long <- data %>%
  pivot_longer(cols = c("blue", "orange"),
               names_to = "Ancestry", values_to = "Proportion") %>%
  mutate(Ancestry = factor(Ancestry, levels = c("orange", "blue")))  # invert stacking

# Step 5: Find positions for black bars
group_positions <- data %>%
  group_by(Group) %>%
  summarize(end = max(IndIndex)) %>%
  pull(end)

###first draft
Strucplot <- ggplot(data_long, aes(x = IndIndex, y = Proportion, fill = Ancestry)) +
  geom_col(width = 1, size = 0.1) +
  scale_fill_manual(values = c("orange" = "darkorange2", "blue" = "dodgerblue4")) +
  geom_vline(xintercept = group_positions + 0.5,
             color = "black", size = 1.1) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 16),  
    axis.title.x = element_text(size = 16),  
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "Individual", y = "Ancestry Proportion", fill = "Ancestry")

Strucplot




#try to add labels
group_labels <- data %>%
  group_by(status, Region) %>%
  summarise(
    start = min(IndIndex),
    end = max(IndIndex)
  ) %>%
  mutate(midpoint = (start + end) / 2)


group_labels$custom_name <- c(
  "BLG", 
  "BLG", 
  "HYB",
  "HYB", 
  "PKS",
  ""
)

# Example: add first line at x = 0.5
all_vlines <- c(0.5, group_positions + 0.5)

Strucplot2 <- ggplot(data_long, aes(x = IndIndex, y = Proportion, fill = Ancestry)) +
  geom_col(width = 1, color = NA) +
  scale_fill_manual(values = c("orange" = "darkorange2", "blue" = "dodgerblue4")) +
  geom_vline(xintercept = all_vlines + 0.5, color = "black", size = 0.5) +
  geom_text(
    data = group_labels,
    aes(x = midpoint, y = -0.07, label = custom_name),
    inherit.aes = FALSE,    # important
    size = 4.5,
    vjust = 1              # text below the bars
  ) +
  scale_y_continuous(limits = c(-0.15, 1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "Individual", y = "Ancestry Proportion")

Strucplot2

#export the plot at a certain dpi 
ggsave("Strucplot2_500dpi.jpg",
       plot = Strucplot2,   # specify the plot object
       width = 11,           # width in inches
       height = 2.5,          # height in inches
       units = "in",
       dpi = 500
)

