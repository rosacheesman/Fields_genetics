setwd("~/Dropbox/PROMENTA/Choice of educational fields/results")

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)

# Load the data
dat <- read_excel("meta_edu_control_201024.xlsx")

# Set custom order for fields
field_order <- c("Natural sciences, mathematics and statistics", "Social sciences, journalism and information",
                 "Arts and humanities", "Information and Communication Technologies (ICTs)", "Education",
                 "Agriculture, forestry, fisheries and veterinary", "Engineering, manufacturing and construction",
                 "Services", "Business, administration and law", "Health and welfare")

# Set custom order for Version with "Unadjusted" first and "EA-adjusted" second
dat$Version <- factor(dat$Version, levels = c('Unadjusted', 'EA-adjusted'))

# Plot with ggplot
ggplot(dat, aes(x = factor(Field, levels = field_order), y = SNP_h2, color = Version, fill = Version)) + 
  geom_point(position = position_dodge(width = 0.5), size = 3) + 
  geom_errorbar(aes(ymin = SNP_h2 - 1.96 * SE, ymax = SNP_h2 + 1.96 * SE), 
                position = position_dodge(width = 0.5), 
                width = 0) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  labs(x = " ", y = "SNP-based heritability", color = NULL, fill = NULL) + 
  scale_color_manual(values = c("Unadjusted" = "dark green", "EA-adjusted" = "light green")) + 
  scale_fill_manual(values = c("Unadjusted" = "dark green", "EA-adjusted" = "light green")) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 70, hjust = 1, size = 11)) + 
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 11)) + 
  guides(
    color = guide_legend(override.aes = list(shape = 21, size = 4, stroke = 1)),  # Set appropriate legend symbol for color
    fill = guide_legend(override.aes = list(shape = 21, size = 4, stroke = 1))   # Same for fill
  )
