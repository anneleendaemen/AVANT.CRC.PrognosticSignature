source("Functions.R")
source("DataLoad_AVANT.R")

####### CODE FOR FIGURES 1A, 1C, 1D, 1E, 1F, 2A, 2C, 2E, 2G, 2H, 3A, 3C, 3D, 3F, 3G, 3I, 3K, 3L, 4D
####### CODE FOR SUPPLEMENTARY FIGURES S1A, S1B, S3A, S3D, S3G, S3I, S5A, S6A-F, S8A-H, S9A-D, S11, S12A


####### FIGURE 1A: DISEASE-FREE SURVIVAL BY TREATMENT ARM ########

x <- 'DFS'
m <- dfs
fit<- survfit(dfs ~ Elasticnet.score.quartiles, data = pData(avant))
# Drawing survival curves
p <- ggsurvplot(fit, 
   title = paste("AVANT"),
   font.main = c(14, "bold", "black"),
   font.x = c(14, "bold", "black"),
   font.y = c(14, "bold", "black"),
   font.tickslab = 14, 
   palette =  c(1:4),
   xlab = "Time (months)", 
   ylab = x,
   conf.int = FALSE, 
   pval=TRUE, 
   # linetype = "strata", change line type based on groups - solid/dashed
   risk.table = TRUE, # Add risk table
   risk.table.col = "strata", # Change risk table color by groups
   risk.table.title = "Numbers at risk",
   risk.table.fontsize = 3, 
   risk.table.y.text = FALSE,
   risk.table.y.text.col = TRUE, 
   risk.table.height = 0.3,
   break.time.by = 10, 
   legend="bottom", 
   legend.title = "Signature Quartiles",
   legend.labs = c('lowest', 'low', 'high', 'highest'),
   ggtheme = theme_bw()) # Change ggplot2 theme
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=12))

pdf("../figures/Figure1A_AVANT_SignatureQuartiles_DFS.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()


####### FIGURE 1C: GENE-GENE CORRELATION OF AVANT SIGNATURE ########

i <- 0.2 # Alpha value
m <- paste0('alpha_',i)
myd <- data.frame(read_excel("../data/AVANT_Elasticnet_GZMBsearch_prognosticgenes.xlsx", sheet=as.character(m)))
mycor <- cor(t(exprs(avant)[myd$gene,]))# compute correlation matrix 

##  Getting order from corrplot dendrogram
mydist <- hclust(dist(1-mycor, method='euclidean'), method='ward.D2')  # This calculates hclust distances exactly as the corrplot function
dend1 <- as.dendrogram(mydist)
dend1 <- color_branches(dend1, k = 4) # COLOR BRANCHES ACCORDING TO GROUP LABEL
     
# Reorder the dendrogram 
dend2 <- dend1 %>% rotate(2:20,1) # swaps the branches of the original dendrogram

# Obtain group labels for each cluster  
groups <- cutree(dend2, k=4) 
myd <- cbind(myd, groups)
myd <- myd[order(myd$groups),]
myd <- myd[order(match(myd$gene, labels(dend2))),]

# Replot the reordered corrplot                                                                         
myd1 <- myd[order(-match(myd$gene, labels(dend2))),]                                                                           
mycor <- cor(t(exprs(avant)[myd1$gene,]))# compute new correlation matrix                                                                                                          

pdf("../figures/Figure1C_AVANT_Signature_GeneGeneCorrelation.pdf", width=6, height=6, onefile=FALSE)
corrplot(mycor, order='original', tl.cex=0.6, tl.col="black") #mar=c(0,0,1,0)) # reordered corrplot      
dev.off()


###### FIGURE 1D, 1E and 1F: DISEASE-FREE SURVIVAL BY AVANT SUB-SIGNATURE ########

fit<- survfit(dfs ~ Immune_signature.Levels, data = pData(avant))
# Drawing survival curves
p <- ggsurvplot(fit, 
   title = paste("AVANT: Immune signature"),
   font.main = c(14, "bold", "black"),
   font.x = c(14, "bold", "black"),
   font.y = c(14, "bold", "black"),
   font.tickslab = 14, 
   palette =  c(1:4),
   xlab = "Time (months)", 
   ylab = x,
   conf.int = FALSE, 
   pval=TRUE, 
   # linetype = "strata", change line type based on groups - solid/dashed
   risk.table = TRUE, # Add risk table
   risk.table.col = "strata", # Change risk table color by groups
   risk.table.title = "Numbers at risk",
   risk.table.fontsize = 3, 
   risk.table.y.text = FALSE,
   risk.table.y.text.col = TRUE, 
   risk.table.height = 0.3,
   #break.time.by = 10, 
   legend="bottom", 
   legend.title = "Immune signature",
   legend.labs = c('Low', 'High'),
   ggtheme = theme_bw()) # Change ggplot2 theme
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=12))

pdf("../figures/Figure1D_AVANT_ImmuneSignature_DFS.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()

x <- 'DFS'
m <- dfs
fit<- survfit(dfs ~ Stromal_signature.Levels, data = pData(avant))
# Drawing survival curves
p <- ggsurvplot(fit, 
   title = paste("AVANT: Stromal signature"),
   font.main = c(14, "bold", "black"),
   font.x = c(14, "bold", "black"),
   font.y = c(14, "bold", "black"),
   font.tickslab = 14, 
   palette =  c(1:4),
   xlab = "Time (months)", 
   ylab = x,
   conf.int = FALSE, 
   pval=TRUE, 
   # linetype = "strata", change line type based on groups - solid/dashed
   risk.table = TRUE, # Add risk table
   risk.table.col = "strata", # Change risk table color by groups
   risk.table.title = "Numbers at risk",
   risk.table.fontsize = 3, 
   risk.table.y.text = FALSE,
   risk.table.y.text.col = TRUE, 
   risk.table.height = 0.3,
   #break.time.by = 10, 
   legend="bottom", 
   legend.title = "Stromal signature",
   legend.labs = c('Low', 'High'),
   ggtheme = theme_bw()) # Change ggplot2 theme
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=12))

pdf("../figures/Figure1E_AVANT_StromalSignature_DFS.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()

fit<- survfit(dfs ~ Proliferative_signature.Levels, data = pData(avant))
# Drawing survival curves
p <- ggsurvplot(fit, 
   title = paste("AVANT: Proliferative signature"),
   font.main = c(14, "bold", "black"),
   font.x = c(14, "bold", "black"),
   font.y = c(14, "bold", "black"),
   font.tickslab = 14, 
   palette =  c(1:4),
   xlab = "Time (months)", 
   ylab = x,
   conf.int = FALSE, 
   pval=TRUE, 
   # linetype = "strata", change line type based on groups - solid/dashed
   risk.table = TRUE, # Add risk table
   risk.table.col = "strata", # Change risk table color by groups
   risk.table.title = "Numbers at risk",
   risk.table.fontsize = 3, 
   risk.table.y.text = FALSE,
   risk.table.y.text.col = TRUE, 
   risk.table.height = 0.3,
   #break.time.by = 10, 
   legend="bottom", 
   legend.title = "Proliferative signature",
   legend.labs = c('Low', 'High'),
   ggtheme = theme_bw()) # Change ggplot2 theme
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=12))

pdf("../figures/Figure1F_AVANT_ProliferativeSignature_DFS.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()


########## FIGURE 2A: DFS BY GZMB #####################

x <- 'DFS'
m <- dfs
fit<- survfit(m ~ GZMB.Levels, data = pData(avant))
# Drawing survival curves
p <- ggsurvplot(fit, 
   title = paste("AVANT: GZMB"),
   font.main = c(14, "bold", "black"),
   font.x = c(14, "bold", "black"),
   font.y = c(14, "bold", "black"),
   font.tickslab = 14, 
   palette =  c(1:4),
   xlab = "Time (months)", 
   ylab = x,
   conf.int = FALSE, 
   pval=TRUE, 
   # linetype = "strata", change line type based on groups - solid/dashed
   risk.table = TRUE, # Add risk table
   risk.table.col = "strata", # Change risk table color by groups
   risk.table.title = "Numbers at risk",
   risk.table.fontsize = 3, 
   risk.table.y.text = FALSE,
   risk.table.y.text.col = TRUE, 
   risk.table.height = 0.3,
   #break.time.by = 10, 
   legend="bottom", 
   legend.title = "GZMB",
   legend.labs = levels(pData(avant)$GZMB.Levels),
   ggtheme = theme_bw()) # Change ggplot2 theme
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5))
# p$table <- p$table + theme(axis.line = element_blank()) to change numbers at risk table

pdf("../figures/Figure2A_GZMBexpression_AVANT.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()


########## FIGURE 2C: DFS BY T-EFFECTOR SIGNATURE #####################

i <- "Teffector_signature"
j <- signatures[i]
genes <- signatures[[i]]
k <- paste0(names(j),".Levels")
x <- 'DFS'
m <- dfs
fit <- survfit(m ~ pData(avant)[,k], data = pData(avant))
# Drawing survival curves
p <- ggsurvplot(fit, 
   title = paste("AVANT: T-effector signature"),
   font.main = c(14, "bold", "black"),
   font.x = c(14, "bold", "black"),
   font.y = c(14, "bold", "black"),
   font.tickslab = 14, 
   palette =  c(1:2),
   xlab = "Time (months)", 
   ylab = x,
   conf.int = FALSE, 
   pval=TRUE, 
   # linetype = "strata", change line type based on groups - solid/dashed
   risk.table = TRUE, # Add risk table
   risk.table.col = "strata", # Change risk table color by groups
   risk.table.title = "Numbers at risk",
   risk.table.fontsize = 5, 
   risk.table.y.text = FALSE,
   risk.table.y.text.col = TRUE, 
   risk.table.height = 0.3,
   #break.time.by = 10, 
   legend="bottom", 
   legend.title = names(j),
   legend.labs = c('Low','High'),
   ggtheme = theme_bw()) # Change ggplot2 theme
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=12))
# p$table <- p$table + theme(axis.line = element_blank()) to change numbers at risk table

pdf(paste0("../figures/Figure2C_",i,"_AVANT.pdf"), width=6, height=6, onefile=FALSE)
print(p)
dev.off()


####### FIGURE 2E: FOREST PLOT FOR T-EFFECTOR AND GZMB FOR DFS #############

x <- 'DFS'
m <- dfs

  # Forest plot
  myforest <- data.frame(matrix(ncol=0,nrow=0))
  cols <- c("Teffector_signature_noGZMB","GZMB")
  for(i in cols){
  # Individual gene signatures prognostic value 
  temp <- with(pData(avant), summary(coxph(m ~ get(i)))) 
  myforest[i,"HR"] <- temp$coefficients[[2]]
  myforest[i,"CI_lower"] <- temp$conf.int[[3]]
  myforest[i,"CI_upper"] <- temp$conf.int[[4]]
  myforest[i,"Cox_Pvalue"] <- signif(temp$coefficients[[5]],2)
  }
  # Added value of Teff to GZMB
   myanova <- anova(
    with(pData(avant), coxph(m ~ GZMB)), 
    with(pData(avant), coxph(m ~ GZMB + Teffector_signature_noGZMB))
    )
  myforest[i,"Teffector signature added\nbenefit to GZMB \n(P-value)"] <- signif(myanova$`P(>|Chi|)`[2],2)
  # Added value of elastic net to individual signatures
    myanova <- anova(
    with(pData(avant), coxph(m ~ Teffector_signature_noGZMB)), 
    with(pData(avant), coxph(m ~ Teffector_signature_noGZMB + GZMB))
    )
  myforest[i,"GZMB added benefit \nto Teffector signature \n(P-value)"] <- signif(myanova$`P(>|Chi|)`[2],2)
  myforest$Group <- rownames(myforest)
   myforest2 <- myforest[2, , drop=FALSE]
  rownames(myforest2)[1] <- ""
  myforest2$Group <- ""

  # Forest plot with p-values
  p <- ggplot(myforest, aes(x = Group, y = HR, ymin = CI_lower, ymax = CI_upper), colour = "black") +
              geom_hline(aes(yintercept = 1), colour = 'grey', linetype = 2) +
              geom_pointrange(size = 0.7, shape = 15) +
              scale_y_log10(limits = c(0.3, 2), breaks = c(0.5,1,1.5)) +
              scale_x_discrete(limits = rev(myforest$Group)) + 
              theme_bw() + 
              ggtitle(paste(x,'\nAVANT Dataset')) + 
              coord_flip() + 
              theme(legend.position="none", 
                title=element_text(face="bold", hjust=0.5),
                plot.title=element_text(face="bold", hjust=0.5),
                axis.text = element_text(size=12,face="bold"))  + 
              annotate("text", x = 2:1, y=1.5, label=paste("P=", myforest$Cox_Pvalue))
tbl <- gridExtra::tableGrob(myforest2[,c(7,5,6)], rows=NULL, theme=ttheme_minimal(colhead=list(fg_params = list(parse=TRUE))))
  # Plot chart and table into one object
 pdf("../figures/Figure2E_ForestPlot_AVANT_DFS.pdf", width=7, height=6, onefile=FALSE)
  gridExtra::grid.arrange(arrangeGrob(p, ncol=1,nrow=1), arrangeGrob(tbl,ncol=1,nrow=1, widths=c(2)),
               as.table=TRUE,
               heights=c(4,1.5))
dev.off()


########### FIGURE 2G and 2H: DFS BY SIGNATURE AND GZMB LEVELS #############

cols <- c("Proliferative_signature.Levels","Stromal_signature.Levels","Immune_signature.Levels")
for(i in cols[1:2]){
  j <- paste0(gsub("_signature.Levels","",i), "_GZMB")
  pData(avant)[,j] <- paste0(gsub("_signature.Levels"," signature",i),"=" ,pData(avant)[,i],"; GZMB=", pData(avant)$GZMB.Levels)
  # Plot KM curves for each signature 
  x <- 'DFS'
  # Make Surv objects
  dfs <- with( 
  	pData(avant),
    Surv(
      time  = TTMDFS,
      event = CSDFS
    )
  )
m <- dfs
fit <- survfit(m ~ get(j), data = pData(avant))
p <- ggsurvplot(fit, 
   title = paste("AVANT:",gsub("_signature.Levels"," signature",i)),
   font.main = c(14, "bold", "black"),
   font.x = c(14, "bold", "black"),
   font.y = c(14, "bold", "black"),
   font.tickslab = 14, 
   palette =  c(1:4),
   xlab = "Time (months)", 
   ylab = x,
   conf.int = FALSE, 
   pval=TRUE, 
   # linetype = "strata", change line type based on groups - solid/dashed
   risk.table = FALSE, # Add risk table
   risk.table.col = "strata", # Change risk table color by groups
   risk.table.title = NULL,
   risk.table.fontsize = 4, 
   risk.table.y.text = NULL,
   risk.table.y.text.col = TRUE, 
   risk.table.height = 0.3,
   # break.time.by = 10, 
   legend="bottom", 
   legend.labs = levels(factor(pData(avant)[[j]])),
   ggtheme = theme_bw()) # Change ggplot2 theme
p$plot <- p$plot + guides(colour=guide_legend(nrow = 4,ncol=1)) + 
          theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=12),legend.title = element_text(size=12)) 
# p$table <- p$table + theme(axis.line = element_blank()) to change numbers at risk table

if (i=="Proliferative_signature.Levels") pdf("../figures/Figure2G_GZMBvsProliferative_AVANT.pdf", width=6, height=6, onefile=FALSE) else if (i=="Stromal_signature.Levels") pdf("../figures/Figure2H_GZMBvsStromal_AVANT.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()
}


###### FIGURE 3A: GZMB EXPRESSION BY CMS SUBTYPE #########

i <- 'GZMB'
  xlabs <- paste(levels(factor(pData(avant)$CMS_Subtypes)),"\n(N=",table(pData(avant)$CMS_Subtypes),")",sep="")
p <- ggplot(data=pData(avant)[pData(avant)$CMS_Subtypes !="unclassified",], aes_string(x = 'CMS_Subtypes', y= as.character(i))) +
  geom_boxplot(aes(fill=factor(CMS_Subtypes))) + ylim(-1.5,5) +
  theme_bw() + 
  labs(title=paste("AVANT:",i), x = '', y = "Expression") +
  theme(plot.title = element_text(lineheight=1.5, face="bold", hjust=0.5, size=20), 
        axis.title = element_text(lineheight=1.5, face="bold", hjust=0.5,size=16), 
        axis.text = element_text(lineheight=1.5, face="bold", hjust=0.5, colour = "black", size=16),
        legend.position="none") + 
  scale_x_discrete(labels=xlabs) 
# pairwise t-test
m <- as.data.frame(with(pData(avant)[pData(avant)$CMS_Subtypes !="unclassified",], pairwise.t.test(get(i), CMS_Subtypes))[["p.value"]]) %>% signif(.,2)
m[is.na(m)] <- "-"

pdf("../figures/Figure3A_GZMBbyCMS_AVANT.pdf", width=6, height=6, onefile=FALSE)
print(p + annotation_custom(tableGrob(m), xmin=3, xmax=3.5, ymin=3.5, ymax=4.5))
dev.off()


###### FIGURE 3C: GZMB EXPRESSION BY MSI #########

i <- "GZMB"
  xlabs <- paste(levels(factor(pData(avant)$MSI)),"\n(N=",table(pData(avant)$MSI),")",sep="")
p <- ggplot(data=pData(avant), aes_string(x = 'MSI', y= as.character(i))) +
  geom_boxplot(aes(fill=factor(MSI))) + 
  ylim(-1.5,5) +
  theme_bw() + 
  labs(title=paste("AVANT:",i), x = '', y = "Expression") +
  theme(plot.title = element_text(lineheight=1.5, face="bold", hjust=0.5, size=20), 
        axis.title = element_text(lineheight=1.5, face="bold", hjust=0.5,size=16), 
        axis.text = element_text(lineheight=1.5, face="bold", hjust=0.5, colour = "black", size=16),
        legend.position="none") + 
  scale_x_discrete(labels=xlabs) + stat_compare_means(label.x.npc="center")

pdf("../figures/Figure3C_GZMBbyMSI_AVANT.pdf", width=6, height=6, onefile=FALSE)
  print(p)
dev.off()


###### FIGURE 3D: T-EFFECTOR SIGNATURE BY CMS SUBTYPE #########

i <- 'Teffector_signature'
  xlabs <- paste(levels(factor(pData(avant)$CMS_Subtypes)),"\n(N=",table(pData(avant)$CMS_Subtypes),")",sep="")
p <- ggplot(data=pData(avant)[pData(avant)$CMS_Subtypes !="unclassified",], aes_string(x = 'CMS_Subtypes', y= as.character(i))) +
  geom_boxplot(aes(fill=factor(CMS_Subtypes))) + ylim(-1.5,5) +
  theme_bw() + 
  labs(title=paste("AVANT:",i), x = '', y = "Expression") +
  theme(plot.title = element_text(lineheight=1.5, face="bold", hjust=0.5, size=20), 
        axis.title = element_text(lineheight=1.5, face="bold", hjust=0.5,size=16), 
        axis.text = element_text(lineheight=1.5, face="bold", hjust=0.5, colour = "black", size=16),
        legend.position="none") + 
  scale_x_discrete(labels=xlabs) 
# pairwise t-test
m <- as.data.frame(with(pData(avant)[pData(avant)$CMS_Subtypes !="unclassified",], pairwise.t.test(get(i), CMS_Subtypes))[["p.value"]]) %>% signif(.,2)
m[is.na(m)] <- "-"

pdf("../figures/Figure3D_TeffectorbyCMS_AVANT.pdf", width=6, height=6, onefile=FALSE)
print(p + annotation_custom(tableGrob(m), xmin=3, xmax=3.5, ymin=3.5, ymax=4.5))
dev.off()


###### FIGURE 3F: T-EFFECTOR SIGNATURE BY MSI #########

i <- "Teffector_signature"
  xlabs <- paste(levels(factor(pData(avant)$MSI)),"\n(N=",table(pData(avant)$MSI),")",sep="")
p <- ggplot(data=pData(avant), aes_string(x = 'MSI', y= as.character(i))) +
  geom_boxplot(aes(fill=factor(MSI))) + 
  ylim(-1.5,5) +
  theme_bw() + 
  labs(title=paste("AVANT:",i), x = '', y = "Expression") +
  theme(plot.title = element_text(lineheight=1.5, face="bold", hjust=0.5, size=20), 
        axis.title = element_text(lineheight=1.5, face="bold", hjust=0.5,size=16), 
        axis.text = element_text(lineheight=1.5, face="bold", hjust=0.5, colour = "black", size=16),
        legend.position="none") + 
  scale_x_discrete(labels=xlabs) + stat_compare_means(label.x.npc="center")

pdf("../figures/Figure3F_TeffectorbyMSI_AVANT.pdf", width=6, height=6, onefile=FALSE)
  print(p)
dev.off()


###### FIGURE 3G: CORRELATION OF GZMA/B AND T-EFFECTOR #########

# Teffector signature vs GGZMA/GZMB in CMS
d <- data.frame(matrix(nrow=4, ncol=3))
names(d) <-c('CMS.subtypes','GZMA','GZMB')
d <-ddply(pData(avant)[pData(avant)$CMS_Subtypes !='unclassified',], "CMS_Subtypes", summarise, corr=cor(GZMB, Teffector_signature_modified))
names(d) <-c('CMS.subtypes','GZMB')
d$GZMA <- ddply(pData(avant)[pData(avant)$CMS_Subtypes !='unclassified',], "CMS_Subtypes", summarise, corr=cor(GZMA, Teffector_signature_modified))[,2]
rownames(d) <- d$CMS.subtypes

# Teffector signature vs GGZMA/GZMB in MSI
d1 <- data.frame(matrix(nrow=4, ncol=3))
names(d1) <-c('MSI','GZMA','GZMB')
d1 <-ddply(pData(avant), "MSI", summarise, corr=cor(GZMB, Teffector_signature_modified))
names(d1) <-c('MSI','GZMB')
d1$GZMA <- ddply(pData(avant), "MSI", summarise, corr=cor(GZMA, Teffector_signature_modified))[,2]
rownames(d1) <- d1$MSI
col1 <- colorRampPalette(c("red","white","blue"))
x <- rbind(d[-1],d1[-1])

pdf("../figures/Figure3G_GZMABvsModifiedTeff_AVANT.pdf", width=6, height=6, onefile=FALSE)
corrplot(t(x[c(2,1)]), is.corr=FALSE, cl.lim=c(0,1), col=col1(20), title = NA, mar=c(0,0,1,0), tl.cex=1.4)
dev.off()


###### FIGURE 3I: KLRK1 EXPRESSION BY CMS SUBTYPE ##########

i <- 'KLRK1'
  xlabs <- paste(levels(factor(pData(avant)$CMS_Subtypes)),"\n(N=",table(pData(avant)$CMS_Subtypes),")",sep="")
p <- ggplot(data=pData(avant)[pData(avant)$CMS_Subtypes !="unclassified",], aes_string(x = 'CMS_Subtypes', y= as.character(i))) +
  geom_boxplot(aes(fill=factor(CMS_Subtypes))) + ylim(-1.5,5) +
  theme_bw() + 
  labs(title=paste("AVANT:",i), x = '', y = "Expression") +
  theme(plot.title = element_text(lineheight=1.5, face="bold", hjust=0.5, size=20), 
        axis.title = element_text(lineheight=1.5, face="bold", hjust=0.5,size=16), 
        axis.text = element_text(lineheight=1.5, face="bold", hjust=0.5, colour = "black", size=16),
        legend.position="none") + 
  scale_x_discrete(labels=xlabs) 
# pairwise t-test
m <- as.data.frame(with(pData(avant)[pData(avant)$CMS_Subtypes !="unclassified",], pairwise.t.test(get(i), CMS_Subtypes))[["p.value"]]) %>% signif(.,2)
m[is.na(m)] <- "-"

pdf("../figures/Figure3I_KLRK1byCMS_AVANT.pdf", width=6, height=6, onefile=FALSE)
print(p + annotation_custom(tableGrob(m), xmin=3, xmax=3.5, ymin=3.5, ymax=4.5))
dev.off()


####### FIGURE 3K: KLRK1 EXPRESSION BY MSI STATUS #########

i <- "KLRK1"
xlabs <- paste(levels(factor(pData(avant)$MSI)),"\n(N=",table(pData(avant)$MSI),")",sep="")
p <- ggplot(data=pData(avant), aes_string(x = 'MSI', y= as.character(i))) +
  geom_boxplot(aes(fill=factor(MSI))) + 
  ylim(-1.5,5) +
  theme_bw() + 
  labs(title=paste("AVANT:",i), x = '', y = "Expression") +
  theme(plot.title = element_text(lineheight=1.5, face="bold", hjust=0.5, size=20), 
        axis.title = element_text(lineheight=1.5, face="bold", hjust=0.5,size=16), 
        axis.text = element_text(lineheight=1.5, face="bold", hjust=0.5, colour = "black", size=16),
        legend.position="none") + 
  scale_x_discrete(labels=xlabs) + stat_compare_means(label.x.npc="center")

pdf("../figures/Figure3K_KLRK1byMSI_AVANT.pdf", width=6, height=6, onefile=FALSE)
  print(p)
dev.off()


########### FIGURE 3L: DFS BY GZMB IN CMS2, T-EFFECTOR LOW TUMORS #############

### Does the association of GZMB with prognosis hold true in the subset of CMS2 tumors with low Teffector? (use of Teffector signature without GZMB, mean cutoff) [mean better cutoff than median in this case, given the large range in Teffector score in CMS1 tumors only, and low levels in all of CMS2/3/4]
x <- 'DFS'
i <- 'CMS2'
  dsub <- subset(pData(avant), CMS_Subtypes==i)
  dsub <- subset(dsub, Teffector_signature_noGZMB.MeanLevels=="Low")
  # Make Surv objects
# Make survival curves 
dfs <- with( 
  dsub,
  Surv(
    time  = TTMDFS,
    event = CSDFS
  )
)
m <- dfs

fit<- survfit(m ~ GZMB.Levels, data = dsub)
p <- ggsurvplot(fit, 
   title = paste0("GZMB\nAVANT dataset (",i,", n=",nrow(dsub),")"),
   font.main = c(14, "bold", "black"),
   font.x = c(14, "bold", "black"),
   font.y = c(14, "bold", "black"),
   font.tickslab = 14, 
   palette = c(1:2),
   xlab = "Time (months)", 
   ylab = x,
   conf.int = FALSE, 
   pval=TRUE, 
   risk.table = FALSE, # Add risk table
   risk.table.col = "strata", # Change risk table color by groups
   risk.table.title = "Numbers at risk",
   risk.table.fontsize = 4, 
   risk.table.y.text = FALSE,
   risk.table.y.text.col = TRUE, 
   risk.table.height = 0.3,
   #break.time.by = 10, 
   legend="bottom", 
   legend.title = "GZMB",
   legend.labs = c('Low','High'),
   ggtheme = theme_bw()) # Change ggplot2 theme
p$plot <- p$plot + theme(
  plot.title = element_text(hjust = 0.5, size=18),
  legend.text=element_text(size=14), 
  legend.title = element_text(size=14))
#   legend.margin=margin(t = 0, unit='cm'),
#   plot.margin=unit(c(1,1,1,1),"cm"))

pdf("../figures/Figure3L_GZMBexpression_CMS2andTeffLow_AVANT.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()


####### FIGURE 4D: pDC SIGNATURE BY CMS SUBTYPE #########

pData(avant)$pDC_signature <- (avant$TCF4 + avant$NRP1)/2
i <- 'pDC_signature'
  xlabs <- paste(levels(factor(pData(avant)$CMS_Subtypes)),"\n(N=",table(pData(avant)$CMS_Subtypes),")",sep="")
p <- ggplot(data=pData(avant)[pData(avant)$CMS_Subtypes !="unclassified",], aes_string(x = 'CMS_Subtypes', y= as.character(i))) +
  geom_boxplot(aes(fill=factor(CMS_Subtypes))) + ylim(-1.5,5) +
  theme_bw() + 
  labs(title=paste0("AVANT: ",i), x = '', y = "Expression") +
  theme(plot.title = element_text(lineheight=1.5, face="bold", hjust=0.5, size=20), 
        axis.title = element_text(lineheight=1.5, face="bold", hjust=0.5,size=16), 
        axis.text = element_text(lineheight=1.5, face="bold", hjust=0.5, colour = "black", size=16),
        legend.position="none") + 
  scale_x_discrete(labels=xlabs) 
# pairwise t-test
m <- as.data.frame(with(pData(avant)[pData(avant)$CMS_Subtypes !="unclassified",], pairwise.t.test(get(i), CMS_Subtypes))[["p.value"]]) %>% signif(.,2)
m[is.na(m)] <- "-"

pdf("../figures/Figure4D_pDCsignaturebyCMS_AVANT.pdf", width=6, height=6, onefile=FALSE)
print(p + annotation_custom(tableGrob(m), xmin=3, xmax=3.5, ymin=3.5, ymax=4.5))
dev.off()


########## FIGURE S1A: OS BY TREATMENT ARM #####################

x <- 'OS'
m <- os
fit<- survfit(m ~ RND, data = pData(avant))
# Drawing survival curves
p <- ggsurvplot(fit, 
   title = paste("Association of treatment with OS"),
   font.main = c(14, "bold", "black"),
   font.x = c(14, "bold", "black"),
   font.y = c(14, "bold", "black"),
   font.tickslab = 14, 
   palette =  c(1:4),
   xlab = "Time (months)", 
   ylab = x,
   conf.int = FALSE, 
   pval=TRUE, 
   # linetype = "strata", change line type based on groups - solid/dashed
   risk.table = TRUE, # Add risk table
   risk.table.col = "strata", # Change risk table color by groups
   risk.table.title = "Numbers at risk",
   risk.table.fontsize = 3, 
   risk.table.y.text = FALSE,
   risk.table.y.text.col = TRUE, 
   risk.table.height = 0.3,
   #break.time.by = 10, 
   legend="bottom", 
   legend.title = "Treatment arms",
   legend.labs = levels(pData(avant)$RND),
   ggtheme = theme_bw()) # Change ggplot2 theme
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5))
# p$table <- p$table + theme(axis.line = element_blank()) to change numbers at risk table

pdf("../figures/FigureS1A_Treatment_AVANT_OS.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()


########## FIGURE S1B: DFS BY TREATMENT ARM #####################

dfs <- with(
  pData(avant),
  Surv(
	time  = TTMDFS,
    event = CSDFS
  )
)

x <- 'DFS'
m <- dfs
fit<- survfit(m ~ RND, data = pData(avant))
# Drawing survival curves
p <- ggsurvplot(fit, 
   title = paste("Association of treatment with DFS"),
   font.main = c(14, "bold", "black"),
   font.x = c(14, "bold", "black"),
   font.y = c(14, "bold", "black"),
   font.tickslab = 14, 
   palette =  c(1:4),
   xlab = "Time (months)", 
   ylab = x,
   conf.int = FALSE, 
   pval=TRUE, 
   # linetype = "strata", change line type based on groups - solid/dashed
   risk.table = TRUE, # Add risk table
   risk.table.col = "strata", # Change risk table color by groups
   risk.table.title = "Numbers at risk",
   risk.table.fontsize = 3, 
   risk.table.y.text = FALSE,
   risk.table.y.text.col = TRUE, 
   risk.table.height = 0.3,
   #break.time.by = 10, 
   legend="bottom", 
   legend.title = "Treatment arms",
   legend.labs = levels(pData(avant)$RND),
   ggtheme = theme_bw()) # Change ggplot2 theme
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5))
# p$table <- p$table + theme(axis.line = element_blank()) to change numbers at risk table

pdf("../figures/FigureS1B_Treatment_AVANT_DFS.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()


####### FIGURE S3A: PREVALENCE OF KRAS MUTATION BY CMS SUBTYPE IN AVANT ########

table(CMS.subtypes, avant$KRAS)
chisq.test(CMS.subtypes, avant$KRAS) # 0.0001223
KRAS <- factor(avant$KRAS, levels=c("MUT", "WT"))

m <- 100*prop.table(xtabs(~ KRAS + CMS.subtypes), 2)

pdf("../figures/FigureS3A_AVANT_KRASvsCMS.pdf", width=6, height=6, onefile=FALSE)
par(mar=c(5.1, 5.1, 4.1, 7.1), xpd=TRUE)
barplot(m, col=c('red','grey'), main='KRAS mutation: AVANT', cex.names = 1.5, cex.axis=1.5, ylab='Percentage', width=2, cex.lab=1.5, las=2)
legend('topright', inset=c(-0.3,0), rownames(m), fill=c('red','grey'))
dev.off()


####### FIGURE S3D: PREVALENCE OF BRAF MUTATION BY CMS SUBTYPE IN AVANT ########

table(CMS.subtypes, avant$BRAF_T1799A)
chisq.test(CMS.subtypes, avant$BRAF_T1799A) # 1.23e-15

BRAF <- factor(avant$BRAF_T1799A, levels=c("MUT", "WT"))

m <- 100*prop.table(xtabs(~ BRAF + CMS.subtypes), 2)

pdf("../figures/FigureS3D_AVANT_BRAFvsCMS.pdf", width=6, height=6, onefile=FALSE)
par(mar=c(5.1, 5.1, 4.1, 7.1), xpd=TRUE)
barplot(m, col=c('red','grey'), main='BRAF mutation: AVANT', cex.names = 1.5, cex.axis=1.5, ylab='Percentage', width=2, cex.lab=1.5, las=2)
legend('topright', inset=c(-0.3,0), rownames(m), fill=c('red','grey'))
dev.off()


####### FIGURE S3G: PREVALENCE OF MSI STATUS BY CMS SUBTYPE IN AVANT ########

table(CMS.subtypes, avant$MSI)
chisq.test(CMS.subtypes, avant$MSI) # p<2.2e-16

MSI <- factor(avant$MSI, levels=c("MSI-H", "MSS"))

m <- 100*prop.table(xtabs(~ MSI + CMS.subtypes), 2)

pdf("../figures/FigureS3G_AVANT_MSIvsCMS.pdf", width=6, height=6, onefile=FALSE)
par(mar=c(5.1, 5.1, 4.1, 7.1), xpd=TRUE)
barplot(m, col=c('red','grey'), main='MSI status: AVANT', cex.names = 1.5, cex.axis=1.5, ylab='Percentage', width=2, cex.lab=1.5, las=2)
legend('topright', inset=c(-0.3,0), rownames(m), fill=c('red','grey'))
dev.off()


####### FIGURE S3I: PREVALENCE OF SIDEDNESS BY CMS SUBTYPE IN AVANT ########

side <- as.character(avant$TYPECOLO)
side[avant$TYPECOLO %in% c("LEFT COLECTOMY", "SIGMOIDECTOMY")] <- "LEFT"
side[avant$TYPECOLO %in% c("RIGHT COLECTOMY")] <- "RIGHT" 
side[avant$TYPECOLO %in% c("OTHER", "TOTAL COLECTOMY", "TRANSVERSE COLECTOMY", "LOWER ANTERIOR RESECTION")] <- NA
table(CMS.subtypes, side)
prop.table(table(CMS.subtypes, side),1)
chisq.test(CMS.subtypes, side) # p<2.2e-16

side <- factor(side, levels=c("LEFT", "RIGHT"))

m <- 100*prop.table(xtabs(~ side + CMS.subtypes), 2)

pdf("../figures/FigureS3I_AVANT_SidednessvsCMS.pdf", width=6, height=6, onefile=FALSE)
par(mar=c(5.1, 5.1, 4.1, 7.1), xpd=TRUE)
barplot(m, col=c('red','grey'), main='Sidedness: AVANT', cex.names = 1.5, cex.axis=1.5, ylab='Percentage', width=2, cex.lab=1.5, las=2)
legend('topright', inset=c(-0.3,0), rownames(m), fill=c('red','grey'))
dev.off()


####### FIGURE S5A: ELASTIC NET GENES IN FUNCTION OF ALPHA PARAMETER ########

test <- list()
mya <- c(0.1,0.2,0.3, 0.4,0.5,0.6,0.7, 0.8,0.9) # Set of Alpha values for the loop below

for(i in mya){
  m <- paste0('alpha_',i)
# Read in the gene signatures for each alpha value

  myd <- data.frame(read_excel('../data/AVANT_Elasticnet_GZMBsearch_prognosticgenes.xlsx', sheet=as.character(m)))
  test[[paste(m)]] <-  paste(myd$gene)
}

test2 <- matrix(ncol=length(test), nrow = length(unique(unlist(test, use.names=FALSE)))) 
rownames(test2) <- unique(unlist(test, use.names=FALSE))
colnames(test2) <- names(test)

for(i in names(test)){
  x <- test[[i]]
  if(x %in% rownames(test2)){
  test2[rownames(test2) %in% x,i] <- 1
  # within(test2[rownames(test2) %in% x, i] <- 1)
}
}
test2[is.na(test2)] <- 0

order.rows <- c("SMAD3", "SPP1", "WBSCR17", "TLL1", "CD28", "REG4", "TCF12", "SMAD9", "GMNN", "CXCL1", "DTX1",
			"CRYAB", "RGS2", "E2F1", "KDM5D", "KDM1A", "IGFBP1", "DTX2", "HEYL", "CDCA5", "ABCC9", "PCSK1",
			"FLT1", "CD36", "CDK6", "TNF", "TAP1", "SHISA5", "IGFBP3", "MAPK3", "EPHA4", "FCRL5", "FOS",
			"CDKN1B", "DNMT1", "CENPM", "GZMB", "RPS6KA1", "CDK2", "SP2", "CXCR4", "AXIN2", "RHOA", "SERPINB13",
			"RPL23", "TP73", "ANGPT1", "SCD5", "POU5F1", "APCDD1", "HBEGF", "PDCD1", "MLH1", "RHOJ", "ZEB2",
			"KDM4B", "DLL4", "LIMCH1", "JMJD1C", "MYBL2", "IL7", "MYB", "HMGCS2", "PRND", "ALK", "TNFRSF9",
			"WNT5B", "AKT3", "KCNMA1", "SHCBP1", "MS4A12", "ARID1A")
test2 <- test2[order.rows, ]

# Plot heatmap
pdf("../figures/FigureS5A_ElasticNetGenes.pdf", width=6, height=8, onefile=FALSE)
aheatmap(test2,Rowv=NA, Colv=NA,color=c("grey", "black"), scale='none', distfun = "man", hclustfun = "average", cexRow=1.2, cexCol=1.2, legend=FALSE)
dev.off()


####### FIGURE S6A, S6B, S6C: DISEASE-FREE SURVIVAL BY TREATMENT ARM ########

# PLOT FOR TREATMENT ARM: FOLFOX-4
avant.sub <- avant[,avant$RND %in% "FOLFOX-4"]

dfs <- with(
  pData(avant.sub),
  Surv(
	time  = TTMDFS,
    event = CSDFS
  )
)

x <- 'DFS'
m <- dfs
fit<- survfit(dfs ~ Elasticnet.score.quartiles, data = pData(avant.sub))
# Drawing survival curves
p <- ggsurvplot(fit, 
   title = paste("AVANT: FOLFOX-4"),
   font.main = c(14, "bold", "black"),
   font.x = c(14, "bold", "black"),
   font.y = c(14, "bold", "black"),
   font.tickslab = 14, 
   palette =  c(1:4),
   xlab = "Time (months)", 
   ylab = x,
   conf.int = FALSE, 
   pval=TRUE, 
   # linetype = "strata", change line type based on groups - solid/dashed
   risk.table = TRUE, # Add risk table
   risk.table.col = "strata", # Change risk table color by groups
   risk.table.title = "Numbers at risk",
   risk.table.fontsize = 3, 
   risk.table.y.text = FALSE,
   risk.table.y.text.col = TRUE, 
   risk.table.height = 0.3,
   break.time.by = 10, 
   legend="bottom", 
   legend.title = "Signature Quartiles",
   legend.labs = c('lowest', 'low', 'high', 'highest'),
   ggtheme = theme_bw()) # Change ggplot2 theme
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=12))

pdf("../figures/FigureS6A_AVANT_FOLFOX4_SignatureQuartiles_DFS.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()

# PLOT FOR TREATMENT ARM: FOLFOX-4+Bevacizumab
avant.sub <- avant[,avant$RND %in% "FOLFOX-4+Bevacizumab"]

dfs <- with(
  pData(avant.sub),
  Surv(
	time  = TTMDFS,
    event = CSDFS
  )
)

x <- 'DFS'
m <- dfs
fit<- survfit(dfs ~ Elasticnet.score.quartiles, data = pData(avant.sub))
# Drawing survival curves
p <- ggsurvplot(fit, 
   title = paste("AVANT: FOLFOX-4 + Bevacizumab"),
   font.main = c(14, "bold", "black"),
   font.x = c(14, "bold", "black"),
   font.y = c(14, "bold", "black"),
   font.tickslab = 14, 
   palette =  c(1:4),
   xlab = "Time (months)", 
   ylab = x,
   conf.int = FALSE, 
   pval=TRUE, 
   # linetype = "strata", change line type based on groups - solid/dashed
   risk.table = TRUE, # Add risk table
   risk.table.col = "strata", # Change risk table color by groups
   risk.table.title = "Numbers at risk",
   risk.table.fontsize = 3, 
   risk.table.y.text = FALSE,
   risk.table.y.text.col = TRUE, 
   risk.table.height = 0.3,
   break.time.by = 10, 
   legend="bottom", 
   legend.title = "Signature Quartiles",
   legend.labs = c('lowest', 'low', 'high', 'highest'),
   ggtheme = theme_bw()) # Change ggplot2 theme
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=12))

pdf("../figures/FigureS6B_AVANT_FOLFOX4_Bevacizumab_SignatureQuartiles_DFS.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()


# PLOT FOR TREATMENT ARM: XELOX+Bevacizumab
avant.sub <- avant[,avant$RND %in% "XELOX+Bevacizumab"]

dfs <- with(
  pData(avant.sub),
  Surv(
	time  = TTMDFS,
    event = CSDFS
  )
)

x <- 'DFS'
m <- dfs
fit<- survfit(dfs ~ Elasticnet.score.quartiles, data = pData(avant.sub))
# Drawing survival curves
p <- ggsurvplot(fit, 
   title = paste("AVANT: XELOX + Bevacizumab"),
   font.main = c(14, "bold", "black"),
   font.x = c(14, "bold", "black"),
   font.y = c(14, "bold", "black"),
   font.tickslab = 14, 
   palette =  c(1:4),
   xlab = "Time (months)", 
   ylab = x,
   conf.int = FALSE, 
   pval=TRUE, 
   # linetype = "strata", change line type based on groups - solid/dashed
   risk.table = TRUE, # Add risk table
   risk.table.col = "strata", # Change risk table color by groups
   risk.table.title = "Numbers at risk",
   risk.table.fontsize = 3, 
   risk.table.y.text = FALSE,
   risk.table.y.text.col = TRUE, 
   risk.table.height = 0.3,
   break.time.by = 10, 
   legend="bottom", 
   legend.title = "Signature Quartiles",
   legend.labs = c('lowest', 'low', 'high', 'highest'),
   ggtheme = theme_bw()) # Change ggplot2 theme
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=12))

pdf("../figures/FigureS6C_AVANT_XELOX_Bevacizumab_SignatureQuartiles_DFS.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()


####### FIGURE S6D, S6E, S6F: OVERALL SURVIVAL BY TREATMENT ARM ########

# PLOT FOR TREATMENT ARM: FOLFOX-4
avant.sub <- avant[,avant$RND %in% "FOLFOX-4"]

os <- with(
  pData(avant.sub),
  Surv(
	time  = TTMDIED,
    event = CSDIED
  )
)

x <- 'OS'
m <- os
fit<- survfit(os ~ Elasticnet.score.quartiles, data = pData(avant.sub))
# Drawing survival curves
p <- ggsurvplot(fit, 
   title = paste("AVANT: FOLFOX-4"),
   font.main = c(14, "bold", "black"),
   font.x = c(14, "bold", "black"),
   font.y = c(14, "bold", "black"),
   font.tickslab = 14, 
   palette =  c(1:4),
   xlab = "Time (months)", 
   ylab = x,
   conf.int = FALSE, 
   pval=TRUE, 
   # linetype = "strata", change line type based on groups - solid/dashed
   risk.table = TRUE, # Add risk table
   risk.table.col = "strata", # Change risk table color by groups
   risk.table.title = "Numbers at risk",
   risk.table.fontsize = 3, 
   risk.table.y.text = FALSE,
   risk.table.y.text.col = TRUE, 
   risk.table.height = 0.3,
   break.time.by = 10, 
   legend="bottom", 
   legend.title = "Signature Quartiles",
   legend.labs = c('lowest', 'low', 'high', 'highest'),
   ggtheme = theme_bw()) # Change ggplot2 theme
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=12))

pdf("../figures/FigureS6D_AVANT_FOLFOX4_SignatureQuartiles_OS.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()

# PLOT FOR TREATMENT ARM: FOLFOX-4+Bevacizumab
avant.sub <- avant[,avant$RND %in% "FOLFOX-4+Bevacizumab"]

os <- with(
  pData(avant.sub),
  Surv(
	time  = TTMDIED,
    event = CSDIED
  )
)

x <- 'OS'
m <- os
fit<- survfit(os ~ Elasticnet.score.quartiles, data = pData(avant.sub))
# Drawing survival curves
p <- ggsurvplot(fit, 
   title = paste("AVANT: FOLFOX-4 + Bevacizumab"),
   font.main = c(14, "bold", "black"),
   font.x = c(14, "bold", "black"),
   font.y = c(14, "bold", "black"),
   font.tickslab = 14, 
   palette =  c(1:4),
   xlab = "Time (months)", 
   ylab = x,
   conf.int = FALSE, 
   pval=TRUE, 
   # linetype = "strata", change line type based on groups - solid/dashed
   risk.table = TRUE, # Add risk table
   risk.table.col = "strata", # Change risk table color by groups
   risk.table.title = "Numbers at risk",
   risk.table.fontsize = 3, 
   risk.table.y.text = FALSE,
   risk.table.y.text.col = TRUE, 
   risk.table.height = 0.3,
   break.time.by = 10, 
   legend="bottom", 
   legend.title = "Signature Quartiles",
   legend.labs = c('lowest', 'low', 'high', 'highest'),
   ggtheme = theme_bw()) # Change ggplot2 theme
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=12))

pdf("../figures/FigureS6E_AVANT_FOLFOX4_Bevacizumab_SignatureQuartiles_OS.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()

# PLOT FOR TREATMENT ARM: XELOX+Bevacizumab
avant.sub <- avant[,avant$RND %in% "XELOX+Bevacizumab"]

os <- with(
  pData(avant.sub),
  Surv(
	time  = TTMDIED,
    event = CSDIED
  )
)

x <- 'OS'
m <- os
fit<- survfit(os ~ Elasticnet.score.quartiles, data = pData(avant.sub))
# Drawing survival curves
p <- ggsurvplot(fit, 
   title = paste("AVANT: XELOX + Bevacizumab"),
   font.main = c(14, "bold", "black"),
   font.x = c(14, "bold", "black"),
   font.y = c(14, "bold", "black"),
   font.tickslab = 14, 
   palette =  c(1:4),
   xlab = "Time (months)", 
   ylab = x,
   conf.int = FALSE, 
   pval=TRUE, 
   # linetype = "strata", change line type based on groups - solid/dashed
   risk.table = TRUE, # Add risk table
   risk.table.col = "strata", # Change risk table color by groups
   risk.table.title = "Numbers at risk",
   risk.table.fontsize = 3, 
   risk.table.y.text = FALSE,
   risk.table.y.text.col = TRUE, 
   risk.table.height = 0.3,
   break.time.by = 10, 
   legend="bottom", 
   legend.title = "Signature Quartiles",
   legend.labs = c('lowest', 'low', 'high', 'highest'),
   ggtheme = theme_bw()) # Change ggplot2 theme
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=12))

pdf("../figures/FigureS6F_AVANT_XELOX_Bevacizumab_SignatureQuartiles_OS.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()


###### FIGURE S8A, S8B, S8C, S8D: AVANT SUB-SIGNATURES BY AVANT SIGNATURE QUARTILES #######

xlabs <- paste(levels(factor(pData(avant)$Elasticnet.score.quartiles)),"\n(N=",table(pData(avant)$Elasticnet.score.quartiles),")",sep="")
xlabs <- gsub(" quartile", "", xlabs)

# Boxplots
for(i in names(mygenes)[-4]){
p <- ggplot(data=pData(avant), aes_string(x = 'Elasticnet.score.quartiles', y= i)) + 
  geom_boxplot(aes(fill=factor(avant$Elasticnet.score.quartiles))) + theme_bw() + 
  labs(title=paste(as.character(i)), x = NULL, y = 'Pathway Z-score') + 
      theme(plot.title = element_text(lineheight=1.2, face="bold"), legend.position="none") + 
      guides(fill = guide_legend(title = "AVANT signature quartiles", title.position = "top")) + theme(axis.text=element_text(size=20, face="bold"),axis.title=element_text(size=20,face="bold"), plot.title = element_text(size=20,face="bold", hjust=0.5)) + scale_y_continuous(limits = c(-2, 4)) + scale_x_discrete(labels= xlabs) + scale_fill_manual(values=c("black","red","green","blue"))

m <- as.data.frame(with(pData(avant), pairwise.t.test(get(i), Elasticnet.score.quartiles))[["p.value"]]) %>% signif(.,2)
m[is.na(m)] <- "-"

if (i=="Proliferative_signature") pdf(paste0("../figures/FigureS8A_",i,"_byElasticNetQuartiles_AVANT.pdf"), width=6, height=6, onefile=FALSE) else if (i=="Immune_signature") pdf(paste0("../figures/FigureS8B_",i,"_byElasticNetQuartiles_AVANT.pdf"), width=6, height=6, onefile=FALSE) else if (i=="Stromal_signature_TGFb_enriched") pdf(paste0("../figures/FigureS8C_",i,"_byElasticNetQuartiles_AVANT.pdf"), width=6, height=6, onefile=FALSE) else if (i=="Stromal_signature_ECM_endocrine") pdf(paste0("../figures/FigureS8D_",i,"_byElasticNetQuartiles_AVANT.pdf"), width=6, height=6, onefile=FALSE)
print(p + annotation_custom(tableGrob(m), xmin=2, xmax=3, ymin=3, ymax=4))
dev.off()
}


###### FIGURE S8E, S8F, S8G, S8H: AVANT SUB-SIGNATURES BY CMS SUBTYPE #######

xlabs <- paste(levels(factor(pData(avant)$CMS_Subtypes)),"\n(N=",table(pData(avant)$CMS_Subtypes),")",sep="")
# Boxplots
for(i in names(mygenes)[-4]){
  test <- pData(avant)[pData(avant)$CMS_Subtypes !='unclassified',]
p <- ggplot(data=test, aes_string(x = 'CMS_Subtypes', y= i)) + 
  geom_boxplot(aes(fill=factor(test$CMS_Subtypes))) + theme_bw() + 
  labs(title=paste(as.character(i)), x = NULL, y = 'Pathway Z-score') + 
      theme(plot.title = element_text(lineheight=1.2, face="bold"), legend.position="none") + 
      guides(fill = guide_legend(title = "CMS subtypes", title.position = "top")) + theme(axis.text=element_text(size=20, face="bold"),axis.title=element_text(size=20,face="bold"), plot.title = element_text(size=20,face="bold", hjust=0.5)) + scale_y_continuous(limits = c(-2, 3)) + scale_x_discrete(labels= xlabs) + scale_fill_manual(values=c("orange","blue","magenta","darkgreen"))

m <- as.data.frame(with(test, pairwise.t.test(get(i), CMS_Subtypes))[["p.value"]]) %>% signif(.,2)
m[is.na(m)] <- "-"

if (i=="Proliferative_signature") pdf(paste0("../figures/FigureS8E_",i,"_byCMS_AVANT.pdf"), width=6, height=6, onefile=FALSE) else if (i=="Immune_signature") pdf(paste0("../figures/FigureS8F_",i,"_byCMS_AVANT.pdf"), width=6, height=6, onefile=FALSE) else if (i=="Stromal_signature_TGFb_enriched") pdf(paste0("../figures/FigureS8G_",i,"_byCMS_AVANT.pdf"), width=6, height=6, onefile=FALSE) else if (i=="Stromal_signature_ECM_endocrine") pdf(paste0("../figures/FigureS8H_",i,"_byCMS_AVANT.pdf"), width=6, height=6, onefile=FALSE)
print(p + annotation_custom(tableGrob(m), xmin=2.5, xmax=3.5, ymin=2, ymax=3))
dev.off()
}


###### FIGURE S9A, S9B, S9C: AVANT SUB-SIGNATURES BY MSI #######

xlabs <- paste(levels(factor(pData(avant)$MSI)),"\n(N=",table(pData(avant)$MSI),")",sep="")

# Boxplots
for(i in names(mygenes)[-c(2,3)]){
p <- ggplot(data=pData(avant), aes_string(x = 'MSI', y= as.character(i))) +
  geom_boxplot(aes(fill=factor(MSI))) + ylim(-2,5) +
  theme_bw() + 
  labs(title=paste(as.character(i)), x = NULL, y = 'Pathway Z-score') + 
  theme(plot.title = element_text(lineheight=1.5, face="bold", hjust=0.5, size=20), 
        axis.title = element_text(lineheight=1.5, face="bold", hjust=0.5,size=16), 
        axis.text = element_text(lineheight=1.5, face="bold", hjust=0.5, colour = "black", size=16),
        legend.position="none") + scale_x_discrete(labels=xlabs) + scale_fill_manual(values=c("green","red"))

m <- with(pData(avant), t.test(get(i)[MSI == "MSI-H"], get(i)[MSI =="MSS"]))$p.value
df1 <- data.frame(a = c(1,1,2,2,2), b = c(4,4.5,4.5,4.5,4)) # defines the significance line on ggplot

if (i=="Proliferative_signature") pdf(paste0("../figures/FigureS9A_",i,"_byMSI_AVANT.pdf"), width=6, height=6, onefile=FALSE) else if (i=="Immune_signature") pdf(paste0("../figures/FigureS9B_",i,"_byMSI_AVANT.pdf"), width=6, height=6, onefile=FALSE) else if (i=="Stromal_signature") pdf(paste0("../figures/FigureS9C_",i,"_byMSI_AVANT.pdf"), width=6, height=6, onefile=FALSE)
print(p + geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 4.8, label = paste("p =", signif(m,2)), size = 5))
dev.off()
}


###### FIGURE S9D: CORRELATION OF AVANT SUB-SIGNATURES #######

my_data <- pData(avant)[,c("Proliferative_signature", "Immune_signature", "Stromal_signature")]
colnames(my_data) <- c("Prolife\nrative","Immune", "Stromal")

pdf("../figures/FigureS9D_CorrelationSignatures_AVANT.pdf", width=5, height=5)
pairs(my_data, lower.panel=panel.smooth, upper.panel=panel.cor)
dev.off()


###### FIGURE S11: COMPARISON OF SIGNATURE CONTRIBUTION TO DFS #########
signatures.names <- names(signatures)
signatures.names[1] <- "OncotypeDx_weighted"

dfs <- with(
  pData(avant),
  Surv(
	time  = TTMDFS,
    event = CSDFS
  )
)

x <- 'DFS'
m <- dfs
    myforest <- data.frame(matrix(ncol=0,nrow=0))
  cols <- c("Elasticnet.score",signatures.names[-c(7,8)]) %>% paste0(., "_std")
  for(i in cols){
  # Individual gene signatures prognostic value 
  temp <- with(pData(avant), summary(coxph(m ~ get(i)))) 
  myforest[i,"HR"] <- temp$coefficients[[2]]
  myforest[i,"CI_lower"] <- temp$conf.int[[3]]
  myforest[i,"CI_upper"] <- temp$conf.int[[4]]
  myforest[i,"Cox_Pvalue"] <- signif(temp$coefficients[[5]],2)
  }
  cols <- signatures.names[-c(7,8)] %>% paste0(., "_std")
  for(i in cols){
  # Added value of individual signatures to elastic net
   myanova <- anova(
    with(pData(avant), coxph(m ~ Elasticnet.score_std)), 
    with(pData(avant), coxph(m ~ Elasticnet.score_std + get(i)))
    )
  #myforest[i,"ANOVA_Chisquare"] <- myanova$Chisq[2]
  myforest[i,"Published Signature added benefit \n to AVANT signature \n (P-value)"] <- signif(myanova$`P(>|Chi|)`[2],2)
  # Added value of elastic net to individual signatures
    myanova <- anova(
    with(pData(avant), coxph(m ~ get(i))), 
    with(pData(avant), coxph(m ~ get(i) + Elasticnet.score_std))
    )
  myforest[i,"AVANT signature added benefit \n to published signatures \n (P-value)"] <- signif(myanova$`P(>|Chi|)`[2],2)
  }
  rownames(myforest) <- gsub("_std", "", rownames(myforest))
  myforest$Group <- rownames(myforest)
  myforest <- myforest[-c(4,5), ] # remove Teffector signature modified
  myforest$Group <- c("AVANT signature", "OncotypeDx", "Teffector signature", "CAF Isella et al", "F-TBRS Calon et al")
    myforest2 <- data.frame(NA, NA, NA, NA, 0.010, 1.3e-4, "CMS subtypes")
  colnames(myforest2) <- colnames(myforest)
  myforest2 <- rbind(myforest, myforest2)
  rownames(myforest2)[nrow(myforest2)] <- "CMS_subtypes"

# Forest plot with p-values
p <- ggplot(myforest, aes(x = Group, y = HR, ymin = CI_lower, ymax = CI_upper), colour = "black") +
            geom_hline(aes(yintercept = 1), colour = 'grey', linetype = 2) +
            geom_pointrange(size = 0.7, shape = 15) +
            scale_y_log10(limits = c(0.3, 3), breaks = c(0.3, 1, 3)) +
            scale_x_discrete(limits = rev(myforest$Group)) + 
            theme_bw() + 
            ggtitle(paste(x,'\nAVANT Dataset')) + 
            coord_flip() + 
            theme(legend.position="none", 
              title=element_text(face="bold", hjust=0.5), 
              plot.title=element_text(face="bold", hjust=0.5),
              axis.text = element_text(size=12,face="bold"))  + 
            annotate("text", x = 5:1, y=2.5, label=paste("P=", myforest$Cox_Pvalue))
tt <- gridExtra::ttheme_minimal(colhead=list(fg_params = list(parse=TRUE)))
tbl <- gridExtra::tableGrob(myforest2[-1,c(7,5,6)], rows=NULL, theme=ttheme_minimal(colhead=list(fg_params = list(parse=TRUE))))

pdf("../figures/FigureS11_ForestPlot_AVANT_OncotypeDx4genes.pdf", width=7.2, height=6.5)
gridExtra::grid.arrange(p, tbl,
             ncol=1,nrow=2,
             as.table=TRUE,
             heights=c(5,3))
dev.off()


######## FIGURE S12A RESULTS ########
x <- 'DFS'
m <- dfs
print(with(pData(avant), coxph(m ~ Elasticnet.score + 
                                    AGE_CAT + 
                                    CEABL + 
                                    ECOGSCAT + 
                                    SEX +  
                                    STRATA)))

