library(dplyr)
library(ggplot2)
library(readxl)
library(dplyr)
library(data.table)
library(FactoMineR)
library(factoextra)

load("Fields_wo_EA_LDSCoutput.RData")

# make cor matrix
cormat<-cov2cor(Fields_wo_EA_LDSCoutput$S)
names <- c("edu_EA4sub","arts_EA4sub","social_EA4sub","business_EA4sub","natural_sci_EA4sub","ict_EA4sub","engineering_EA4sub","agri_EA4sub","health_EA4sub","services_EA4sub",
           "edu_EA_adj","arts_EA_adj","social_EA_adj","business_EA_adj","natural_sci_EA_adj","ict_EA_adj","engineering_EA_adj","agri_EA_adj","health_EA_adj","services_EA_adj")
colnames(cormat)<-names
rownames(cormat)<-names
cormat<-round(cormat,2)
# subset to gwas by sub variables
eacor<-cormat[,c("edu_EA4sub","arts_EA4sub","social_EA4sub","business_EA4sub","natural_sci_EA4sub","ict_EA4sub","engineering_EA4sub","agri_EA4sub","health_EA4sub","services_EA4sub")]
eacor<-eacor[c("edu_EA4sub","arts_EA4sub","social_EA4sub","business_EA4sub","natural_sci_EA4sub","ict_EA4sub","engineering_EA4sub","agri_EA4sub","health_EA4sub","services_EA4sub"),]

colnames(eacor)<-c("Education","Arts and humanities","Social sciences, journalism and information","Business, administration and law","Natural sciences, mathematics and statistics","Information and Communication Technologies (ICTs)","Engineering, manufacturing and construction","Agriculture, forestry, fisheries and veterinary","Health and welfare","Services")
rownames(eacor)<-c("Education","Arts and humanities","Social sciences, journalism and information","Business, administration and law","Natural sciences, mathematics and statistics","Information and Communication Technologies (ICTs)","Engineering, manufacturing and construction","Agriculture, forestry, fisheries and veterinary","Health and welfare","Services")
cormatrix<-eacor
# .......................................................................................................

# perform PCA
# and plot PC1 and 2 coords


res.pca <- PCA(cormatrix, scale.unit = TRUE)

res.pca$var$coord <- -res.pca$var$coord
res.pca$ind$coord <- -res.pca$ind$coord

print(res.pca$eig)
print(res.pca$var$coord)

fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("black", "black", "black"),
             repel = TRUE     # Avoid text overlapping
)+labs(x = "Principal Component 1", y = "Principal Component 2", title = "") +
  theme(
    panel.background = element_blank(),          # Removes the gray background
    plot.background = element_rect(fill = "white", color = NA),  # Sets the overall plot background to white
    panel.border = element_blank(),              # Removes the border around the plot area
    panel.grid.major = element_blank(),          # Removes major grid lines
    panel.grid.minor = element_blank(),          # Removes minor grid lines
    axis.line = element_line(color = "black")    # Optionally adds axis lines for clarity
  ) +
  guides(color = "none")  # Remove color legends
print(res.pca$var)

# pca_result <- princomp(cormatrix, cor=T)
# summary(pca_result)
# pca_result$loadings

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
load("Fields_wo_EA_LDSCStand.RData")

paLDSC(S_Stand = Fields_wo_EA_LDSCStand$S_Stand, V_Stand = Fields_wo_EA_LDSCStand$V_Stand, r = 100, p = .95, diag = F,
       fa = F, fm = "minres", save.pdf = T)

