library(dplyr)
library(ggplot2)
library(ggrepel)  
library(data.table)

load("Fields_wo_EA_LDSCStand.RData")

# .......................................................................................................
#check names and order
# .......................................................................................................
# perform PCA

# Apply eigen decomposition
eigen_result <- eigen(cormatrix)

# Calculate loadings for PC1 and PC2
loadings_pc1 <- eigen_result$vectors[, 1] * sqrt(eigen_result$values[1])
loadings_pc2 <- eigen_result$vectors[, 2] * sqrt(eigen_result$values[2])

# Create a dataframe for plotting
loadings_df <- data.frame(
  Variable = rownames(cormatrix),
  PC1 = loadings_pc1,
  PC2 = loadings_pc2
)

# .......................................................................................................
# plot PC1 and 2 coords

# Create circle 
theta <- seq(0, 2*pi, length.out = 100)
circle <- data.frame(x = cos(theta), y = sin(theta))

ggplot() +
  geom_path(data = circle, aes(x = x, y = y), color = "gray") +
  geom_point(data = loadings_df, aes(x = PC1, y = PC2), color = "black") +
  geom_text_repel(data = loadings_df, aes(x = PC1, y = PC2, label = Variable), 
                 color = "black", size = 3) +
  geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  coord_fixed() +
  xlim(-1.1, 1.1) +
  ylim(-1.1, 1.1) +
  labs(x = "Principal Component 1", y = "Principal Component 2", title = "") +
  theme(
    panel.background = element_blank(),         
    plot.background = element_rect(fill = "white", color = NA),  
    panel.border = element_blank(),              
    panel.grid.major = element_blank(),          
    panel.grid.minor = element_blank(),          
    axis.line = element_line(color = "black")    
  ) +
  guides(color = "none")

# .......................................................................................................

# correlation matrix plot

get_upper_tri <- function(cormatrix){
  cormatrix[lower.tri(cormatrix)]<- NA
  return(cormatrix)
}
reorder_cormat <- function(cormatrix){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormatrix)/2)
  hc <- hclust(dd)
  cormatrix <-cormatrix[hc$order, hc$order]
}

# Reorder the correlation matrix
cormatrix <- reorder_cormat(cormatrix)
upper_tri <- get_upper_tri(cormatrix)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="LDSC\nGenetic\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 9, hjust = 1))+
  coord_fixed()
# Print the heatmap
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))


# .......................................................................................................
# Perform parallel analysis to figure out number of PCs to use

try(source("Parallel_Anallysis_paLDSC_JF.R"))

paLDSC(S_Stand = Fields_wo_EA_LDSCStand$S_Stand, V_Stand = Fields_wo_EA_LDSCStand$V_Stand, r = 100, p = .95, diag = F,
       fa = F, fm = "minres", save.pdf = T)

