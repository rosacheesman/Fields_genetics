library(dplyr)
library(data.table)
library(factoextra)
library(FactoMineR)
library(readxl)
library(Matrix)
library(ggplot2)
library(readxl)
library(tidyr)
library(dplyr)

data <- read_excel("~/Dropbox/PROMENTA/Choice of educational fields/results/lifelines/figdat_120324.xlsx")
data<-na.omit(data)


# Convert Effect and SE to numeric
data$Effect <- as.numeric(data$Effect)
data$SE <- as.numeric(data$SE)



# Make Eff_type a factor
data$Eff_type <- factor(data$Eff_type, levels = c("Pop", "Direct", "Parental"))
levels(data$Eff_type)[levels(data$Eff_type) == "Pop"] <- "Population"

table(data$Field)
field_order <- c("Natural sciences, mathematics and statistics","Social sciences, journalism and information",
                 "Arts and humanities","Information and Communication Technologies (ICTs)","Education",
                 "Agriculture, forestry, fisheries and veterinary","Engineering, manufacturing and construction",
                 "Services","Business, administration and law","Health and welfare")


# Convert Field to factor with custom order
data$Field <- factor(data$Field, levels = field_order)



# Create the scatter plot
# Calculate significant points
data$Significance <- ifelse((data$Effect - 2.807 * data$SE) > 0 | (data$Effect + 2.807 * data$SE) < 0, "*", "")

# Plot with significance indication
ggplot(data, aes(x = Field, y = Effect, color = Eff_type, fill = Eff_type, label = Significance)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = Effect - 2.807 * SE, ymax = Effect + 2.807 * SE), 
                position = position_dodge(width = 0.5), 
                width = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_text(position = position_dodge(width = 0.5), vjust = -0.5, hjust = -0.5, size = 4) +  # Adjust size
  labs(x = "Field", y = "Field-specific PGI association (beta)", 
       color = "Genetic effect", fill = "Genetic effect") +
  scale_color_manual(values = c("red", "darkblue", "lightblue")) +
  scale_fill_manual(values = c("red", "darkblue", "lightblue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 65, hjust = 1,size=11))  

theme(axis.text=element_text(size=12),
      axis.title=element_text(size=11,face="bold"))
