
setwd("~/summary_stats/meta")
load("Fields_META_LDSCoutput.RData")
# .......................................................................................................
broad_trait_names<-c("edu","arts","social","business","natural_sci","ict","engineering","agri","health","services",
                     "edu_EA_adj","arts_EA_adj","social_EA_adj","business_EA_adj","natural_sci_EA_adj","ict_EA_adj","engineering_EA_adj","agri_EA_adj","health_EA_adj","services_EA_adj")

covmat<-LDSCoutput2$S
cormat<-cov2cor(LDSCoutput2$S)
names <- broad_trait_names

colnames(covmat)<-names
rownames(covmat)<-names        
colnames(cormat)<-names
rownames(cormat)<-names
cormat<-round(cormat,2)


# .......................................................................................................
# if looking at EA-adjusted
eacor<-cormat[,c("edu_EA_adj","arts_EA_adj","social_EA_adj","business_EA_adj","natural_sci_EA_adj","ict_EA_adj","engineering_EA_adj","agri_EA_adj","health_EA_adj","services_EA_adj")]
eacor<-eacor[c("edu_EA_adj","arts_EA_adj","social_EA_adj","business_EA_adj","natural_sci_EA_adj","ict_EA_adj","engineering_EA_adj","agri_EA_adj","health_EA_adj","services_EA_adj"),]
cormat<-eacor


# .......................................................................................................
# actually give better names
names<- c("Education","Arts and humanities","Social sciences, journalism and information","Business, administration and law","Natural sciences, mathematics and statistics","Information and Communication Technologies (ICTs)","Engineering, manufacturing and construction","Agriculture, forestry, fisheries and veterinary","Health and welfare","Services"
)
colnames(cormat)<-names
rownames(cormat)<-names
# .......................................................................................................

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
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

# ...........................................................................................................
# .......................................................................................................
# PCA
res.pca <- PCA(cormat, scale.unit = TRUE, ncp = 5, graph = TRUE)

print(res.pca$eig)
fviz_eig(res.pca)
print(res.pca$var$coord)
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
