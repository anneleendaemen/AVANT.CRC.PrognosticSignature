source("Load_Libraries.R")

### CODE FOR FIGURES 4A, 4B, 4C
### CODE FOR SUPPLEMENTARY FIGURE S14A, S14B, S14C


load("../data/tsneRun_CD45_immune_allcells_excludedmarkers_noPBMC_dataset.RData")

clusters <- data.frame("Cluster"=c(0:25),"celltype"=rep(NA,26))
clusters$celltype[clusters$Cluster==0] <- "Outliers"
clusters$celltype[clusters$Cluster==1] <- "CD16+ NK cells" # "Mature NK cells" -> CD16+ NK cells
clusters$celltype[clusters$Cluster==2] <-  "Monocytes"# "Myeloid_5"
clusters$celltype[clusters$Cluster==3] <- "GammaDelta_Tcells" 
clusters$celltype[clusters$Cluster==4] <- "CD4_Tcells"
clusters$celltype[clusters$Cluster==5] <- "Outliers" # "Bcells_Naive" -> Outliers
clusters$celltype[clusters$Cluster==6] <- "Outliers" # "CD103a_DCs" -> Outliers
clusters$celltype[clusters$Cluster==7] <- "Outliers"
clusters$celltype[clusters$Cluster==8] <- "GammaDelta_Tcells"
clusters$celltype[clusters$Cluster==9] <- "Granulocytes"# "Myeloid_7"
clusters$celltype[clusters$Cluster==10] <- "pDCs"
clusters$celltype[clusters$Cluster==11] <- "CD8_Tcells"
clusters$celltype[clusters$Cluster==12] <- "Outliers"
clusters$celltype[clusters$Cluster==13] <- "Bcells"  # *** includes the middle junk cluster as well. Figure out some way to clean this up later. 
clusters$celltype[clusters$Cluster==14] <- "Granulocytes"# "Myeloid_6"
clusters$celltype[clusters$Cluster==15] <- "Monocytes"# "Myeloid_4"
clusters$celltype[clusters$Cluster==16] <- "Myeloid" #"Myeloid_1"
clusters$celltype[clusters$Cluster==17] <-   "Outliers" # "NKcells" --> Outliers
clusters$celltype[clusters$Cluster==18] <- "Outliers"
clusters$celltype[clusters$Cluster==19] <- "Outliers" # "NKcells" --> Outliers
clusters$celltype[clusters$Cluster==20] <- "NK cells" # "Intraepithelial lymphocytes" --> NK cells
clusters$celltype[clusters$Cluster==21] <- "CD4_Treg"
clusters$celltype[clusters$Cluster==22] <- "Outliers"
clusters$celltype[clusters$Cluster==23] <- "Myeloid" # "Myeloid_3" 
clusters$celltype[clusters$Cluster==24] <- "Myeloid" # "Myeloid_2"
clusters$celltype[clusters$Cluster==25] <- "Outliers"

# Add cell type information 
df$Celltype <- clusters$celltype[match(df$dbclusters25,clusters$Cluster)]


########################## Figure 4A: tSNE CYTOF PLOT ##########################

# coordinates for cell type annotation
coord <- data.frame(
  x = c(11, 5, -4, -19, -22, -10, -13, -9, -9, 3, -3, 8 , -10), 
  y = c(11, -10, -6, 6, -1, 2.5, 14, -17, -13, 10, 8 , -1, -9), 
  label=c("B-cells","CD8 T-cells","Monocytes","Myeloid","Mono\ncytes","CD4 T-cells", "pDCs", "NK cells", "CD16+\nNK-cells", "Granulocytes", "Granulocytes", "gd\nT-cells", "Treg"), 
  size = c(7,7,5,7, 5, 8, 5, 3, 3, 2, 2, 2, 2)
  )

mycols <- setNames(hue_pal()(12), levels(factor(df$Celltype)))

# Colored for all clusters
p <- ggplot(df) + 
geom_point(aes(x=tSNE.1, y=tSNE.2,colour=factor(Celltype, exclude="Outliers")), size=1, alpha=1) +  
scale_color_manual(values=mycols) + 
theme_bw() + 
ylim(-22,22) + 
xlim(-22,22)
# Plot with legend
p + annotate("text", x = coord[,"x"], y = coord[,"y"], label = coord[,"label"], size=coord[,"size"], fontface='bold', hjust=0.4) + 
theme(legend.position="bottom", 
      axis.title=element_text(size=16, face="bold",hjust=0.5), 
      axis.text=element_text(size=14)) + 
guides(color = guide_legend(nrow=3, ncol=9)) 
# plot without legend

pdf("../figures/Figure4A_tSNE_OutliersExcluded.pdf", width=6, height=6, onefile=FALSE)
print(p + annotate("text", x = coord[,"x"], y = coord[,"y"], label = coord[,"label"], size=coord[,"size"], fontface='bold', hjust=0.4) + 
theme(legend.position="none", 
      axis.title=element_text(size=16, face="bold",hjust=0.5), 
      axis.text=element_text(size=14)))
dev.off()


###################### FIGURE 4B: BI-AXAL PLOT OF CD8 VS. GZMB IN CD16+ NK CELLS #############

colfunc <- colorRampPalette(c("white", "lightblue", "green", "yellow", "red"))

j <- "CD16+ NK cells"
  v2 <- ggplot(subset(df, Celltype==j), aes(Granzyme_B, CD8)) +
  stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
  scale_fill_gradientn(colours=colfunc(400)) + # gives the colour plot
  geom_density2d(colour="black", bins=10) + # draws the lines inside
  theme_bw() + 
  xlim(1,8) + 
  ylim(1,8) +
  geom_vline(xintercept=3, lty=2) +
  geom_hline(yintercept=3, lty=2) +
  ggtitle(paste(j)) + 
  theme(axis.text = element_text(size=14), 
    axis.title = element_text(size=16, face="bold"), 
    plot.title = element_text(size=20, face='bold', hjust=0.5)) 

pdf(paste0("../figures/Figure4B_Biaxial_",j,".pdf"), width=6, height=6, onefile=FALSE)
print(v2)
dev.off()


###################### FIGURE 4C: BI-AXAL PLOT OF CD8 VS. GZMB IN pDCs #############

colfunc <- colorRampPalette(c("white", "lightblue", "green", "yellow", "red"))

j <- "pDCs"
  v2 <- ggplot(subset(df, Celltype==j), aes(Granzyme_B, CD8)) +
  stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
  scale_fill_gradientn(colours=colfunc(400)) + # gives the colour plot
  geom_density2d(colour="black", bins=10) + # draws the lines inside
  theme_bw() + 
  xlim(1,8) + 
  ylim(1,8) +
  geom_vline(xintercept=3, lty=2) +
  geom_hline(yintercept=3, lty=2) +
  ggtitle(paste(j)) + 
  theme(axis.text = element_text(size=14), 
    axis.title = element_text(size=16, face="bold"), 
    plot.title = element_text(size=20, face='bold', hjust=0.5)) 

pdf(paste0("../figures/Figure4C_Biaxial_",j,".pdf"), width=6, height=6, onefile=FALSE)
print(v2)
dev.off()


########################## Figure S14A: HEATMAP OF IMMUNE MARKERS IN THE CYTOF DATASET ##########################

dfsummary_mod <- df[-c(1,2,37,42:45)] %>% group_by(Celltype) %>% summarise_each(funs(mean)) %>% data.frame(., check.names=FALSE) 
rownames(dfsummary_mod) <- dfsummary_mod[,c(1)]
# add original cluster information
dfsummary_mod$Cluster <- clusters$Cluster[match(dfsummary_mod$Celltype,clusters$celltype)]

# Remove outliers, CD103a_DCs, Bcells_Naive and NKcells
test <- dfsummary_mod[!dfsummary_mod$Celltype %in% c("Outliers","CD103a_DCs", "Bcells_Naive", "NKcells"),]
# Rename Intraepithelial lymphocytes to NK cells
rownames(test)[rownames(test)=="Intraepithelial lymphocytes"] <- "NK cells"
# Rename Mature NK cells to CD16+ NK cells
rownames(test)[rownames(test)=="Mature NK cells"] <- "CD16+ NK cells"
test$Celltype <- rownames(test)

# Heatmap
temp <- as.matrix(test[c(2:36,39)])
temp <- scale(temp)

pdf("../figures/FigureS14A_CytofHeatmap_CellTypeMarkers_ScaledPerMarker.pdf", width=8, height=8, onefile=FALSE)
aheatmap(temp, distfun = "pearson", hclustfun = "ward", scale="none", cexRow=8, cexCol=4,
			Rowv=FALSE, Colv=FALSE)
dev.off()


########################## Figure S14B: PREVALENCE OF IMMUNE POPULATIONS ##########################

indices.rm <- which(df$Celltype %in% c("Outliers", "CD103a_DCs", "Bcells_Naive", "NKcells"))
df<- df[-indices.rm, ]
df$Celltype[which(df$Celltype %in% "Intraepithelial lymphocytes")] <- "NK cells"
df$Celltype[which(df$Celltype %in% "Mature NK cells")] <- "CD16+ NK cells"

df$Position <- with(df, factor(Celltype, levels=names(sort(table(Celltype), decreasing=TRUE))))
mycols <- setNames(hue_pal()(15), levels(factor(df$Celltype)))

p <- ggplot(df, aes(fct_infreq(Position, ordered=TRUE), fill=Celltype)) + 
geom_bar(aes(y = (..count..)/sum(..count..))) +   
geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))), stat = "count", hjust = -0.1) + 
scale_y_continuous(labels = percent, limits=c(0,0.4)) + 
coord_flip() + scale_color_manual(values=mycols) +
theme_bw() + 
labs(x = "", y="Percentage") + 
theme(axis.text = element_text(size=13, face='bold'), 
      axis.title = element_text(size=16, face='bold'), 
      legend.position='none')
      
pdf("../figures/FigureS14B_PrevalenceImmunePopulations.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()


###################### FIGURE S14C: BI-AXAL PLOT OF CD8 VS. GZMB IN CD8 T CELLS #############

colfunc <- colorRampPalette(c("white", "lightblue", "green", "yellow", "red"))

j <- "CD8_Tcells"
  v2 <- ggplot(subset(df, Celltype==j), aes(Granzyme_B, CD8)) +
  stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
  scale_fill_gradientn(colours=colfunc(400)) + # gives the colour plot
  geom_density2d(colour="black", bins=10) + # draws the lines inside
  theme_bw() + 
  xlim(1,8) + 
  ylim(1,8) +
  geom_vline(xintercept=3, lty=2) +
  geom_hline(yintercept=3, lty=2) +
  ggtitle(paste(j)) + 
  theme(axis.text = element_text(size=14), 
    axis.title = element_text(size=16, face="bold"), 
    plot.title = element_text(size=20, face='bold', hjust=0.5)) 

pdf(paste0("../figures/FigureS14C_Biaxial_",j,".pdf"), width=6, height=6, onefile=FALSE)
print(v2)
dev.off()
