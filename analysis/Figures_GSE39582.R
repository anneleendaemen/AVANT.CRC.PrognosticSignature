source("Functions.R")
source("DataLoad_GSE39582.R")

####### CODE FOR FIGURES 1B, 1G, 2B, 2D, 2F, 2I, 2J, 3B, 3E, 3H, 3J, 4E
####### CODE FOR SUPPLEMENTARY FIGURES S3C, S3F, S7A, S7B, S10A-C, S12B, S13


############## FIGURE 1B: RECURRENCE-FREE SURVIVAL BY SIGNATURE QUANTILE #############

x <- 'RFS'
m <- rfs
fit<- survfit(rfs ~ glmnetU.quartiles, data = pData(eset))
# Drawing survival curves
p <- ggsurvplot(fit, 
   title = paste("GSE39582"),
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
   legend.title = "Prediction score Quartiles",
   legend.labs = c('lowest', 'low','high','highest'),
   ggtheme = theme_bw()) # Change ggplot2 theme
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=12))

pdf("../figures/Figure1B_AVANTsignature_GSE39582_UnweightedScores.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()


############## FIGURE 1G: FOREST PLOT FOR RFS #############

signatures.names <- names(signatures)
signatures.names[1] <- "OncotypeDx_weighted"

x <- 'RFS'
m <- rfs
    myforest <- data.frame(matrix(ncol=0,nrow=0))
  cols <- c("glmnet.unweighted",signatures.names[-c(7,8)]) %>% paste0(., "_std")
  for(i in cols){
  # Individual gene signatures prognostic value 
  temp <- with(pData(eset), summary(coxph(m ~ get(i)))) 
  myforest[i,"HR"] <- temp$coefficients[[2]]
  myforest[i,"CI_lower"] <- temp$conf.int[[3]]
  myforest[i,"CI_upper"] <- temp$conf.int[[4]]
  myforest[i,"Cox_Pvalue"] <- signif(temp$coefficients[[5]],2)
  }
  cols <- signatures.names[-c(7,8)] %>% paste0(., "_std")
  for(i in cols){
  # Added value of individual signatures to elastic net
   myanova <- anova(
    with(pData(eset), coxph(m ~ glmnet.unweighted_std)), 
    with(pData(eset), coxph(m ~ glmnet.unweighted_std + get(i)))
    )
  myforest[i,"Published signature added benefit \n to AVANT signature \n (P-value)"] <- signif(myanova$`P(>|Chi|)`[2],2)
  # Added value of elastic net to individual signatures
    myanova <- anova(
    with(pData(eset), coxph(m ~ get(i))), 
    with(pData(eset), coxph(m ~ get(i) + glmnet.unweighted_std))
    )
  myforest[i,"AVANT signature added benefit \n to published signatures \n (P-value)"] <- signif(myanova$`P(>|Chi|)`[2],2)
  }
  rownames(myforest) <- gsub("_std", "", rownames(myforest))
  myforest$Group <- rownames(myforest)
  myforest <- myforest[-c(4,5),]
  myforest$Group <- c("AVANT signature", "OncotypeDx", "Teffector signature", "CAF Isella et al", "F-TBRS Calon et al")
  myforest2 <- data.frame(NA, NA, NA, NA, 0.033, 1.1e-4, "CMS subtypes")
  colnames(myforest2) <- colnames(myforest)
  myforest2 <- rbind(myforest, myforest2)
  rownames(myforest2)[nrow(myforest2)] <- "CMS_subtypes"
# Forest plot with p-values
p <- ggplot(myforest, aes(x = Group, y = HR, ymin = CI_lower, ymax = CI_upper), colour = "black") +
            geom_hline(aes(yintercept = 1), colour = 'grey', linetype = 2) +
            geom_pointrange(size = 0.7, shape = 15) +
            scale_y_log10(limits = c(0.3, 3), breaks = c(0.3, 1, 3)) +
            scale_x_discrete(limits = rev(myforest$Group)) + theme_bw() + 
            ggtitle(paste(x, "\n Validation Dataset: GSE39582")) + 
            coord_flip() + 
            theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5), axis.text = element_text(size=12,face="bold"))  + annotate("text", x = 5:1, y=2.5, label=paste("P=", myforest$Cox_Pvalue))
tt <- gridExtra::ttheme_minimal(colhead=list(fg_params = list(parse=TRUE)))
tbl <- gridExtra::tableGrob(myforest2[-1,c(7,5,6)], rows=NULL, theme=tt)
# Plot chart and table into one object

pdf("../figures/Figure1G_ForestPlot_GSE39582_UnweightedScores_OncotypeDx7genes_WEIGHTED.pdf", width=7.2, height=6.5)
gridExtra::grid.arrange(p, tbl,
             ncol=1,nrow=2,
             as.table=TRUE,
             heights=c(5,3))
dev.off()


##################### FIGURE 2B: RFS BY GZMB EXPRESSION #####################

x <- 'RFS'
m <- rfs
fit<- survfit(m ~ GZMB.Levels, data = pData(eset))
# Drawing survival curves
p <- ggsurvplot(fit, 
   title = paste("GSE39582: GZMB"),
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
   risk.table.fontsize = 4, 
   risk.table.y.text = FALSE,
   risk.table.y.text.col = TRUE, 
   risk.table.height = 0.3,
   #break.time.by = 10, 
   legend="bottom", 
   legend.title = "GZMB",
   legend.labs = levels(pData(eset)$GZMB.Levels),
   ggtheme = theme_bw()) # Change ggplot2 theme
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5))

pdf("../figures/Figure2B_GZMBexpression_GSE39582.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()


##################### FIGURE 2D: RFS BY T-EFFECTOR SIGNATURE #####################

x <- 'RFS'
m <- rfs
fit <- survfit(m ~ Teffector_signature.levels, data = pData(eset))
# Drawing survival curves
p <- ggsurvplot(fit, 
   title = paste("GSE39582: T-effector signature"),
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
   risk.table.fontsize = 4, 
   risk.table.y.text = FALSE,
   risk.table.y.text.col = TRUE, 
   risk.table.height = 0.3,
   #break.time.by = 10, 
   legend="bottom", 
   legend.title = "T-effector signature",
   legend.labs = levels(pData(eset)$Teffector_signature.levels),
   ggtheme = theme_bw()) # Change ggplot2 theme
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5))

pdf("../figures/Figure2D_Teffector_signature_GSE39582.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()


############## FIGURE 2F: FOREST PLOT FOR T EFFECTOR AND GZMB #############

x <- 'RFS'
m <- rfs
  # Forest plot
  myforest <- data.frame(matrix(ncol=0,nrow=0))
  cols <- c("Teffector_signature_noGZMB","GZMB")
  for(i in cols){
  # Individual gene signatures prognostic value 
  temp <- with(pData(eset), summary(coxph(m ~ get(i)))) 
  myforest[i,"HR"] <- temp$coefficients[[2]]
  myforest[i,"CI_lower"] <- temp$conf.int[[3]]
  myforest[i,"CI_upper"] <- temp$conf.int[[4]]
  myforest[i,"Cox_Pvalue"] <- signif(temp$coefficients[[5]],2)
  }
  # Added value of Teff to GZMB
   myanova <- anova(
    with(pData(eset), coxph(m ~ GZMB)), 
    with(pData(eset), coxph(m ~ GZMB + Teffector_signature_noGZMB))
    )
  myforest[i,"Teffector signature added\nbenefit to GZMB \n(P-value)"] <- signif(myanova$`P(>|Chi|)`[2],2)
  # Added value of elastic net to individual signatures
    myanova <- anova(
    with(pData(eset), coxph(m ~ Teffector_signature_noGZMB)), 
    with(pData(eset), coxph(m ~ Teffector_signature_noGZMB + GZMB))
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
              ggtitle(paste(x,'\nGSE39582 Dataset')) + 
              coord_flip() + 
              theme(legend.position="none", 
                title=element_text(face="bold", hjust=0.5),
                plot.title=element_text(face="bold", hjust=0.5),
                axis.text = element_text(size=12,face="bold"))  + 
              annotate("text", x = 2:1, y=1.5, label=paste("P=", myforest$Cox_Pvalue))
tbl <- gridExtra::tableGrob(myforest2[,c(7,5,6)], rows=NULL, theme=ttheme_minimal(colhead=list(fg_params = list(parse=TRUE))))
  # Plot chart and table into one object

pdf("../figures/Figure2F_ForestPlot_GSE39582_RFS.pdf", width=7, height=6, onefile=FALSE)
  gridExtra::grid.arrange(arrangeGrob(p, ncol=1,nrow=1), arrangeGrob(tbl,ncol=1,nrow=1, widths=c(2)),
               as.table=TRUE,
               heights=c(4,1.5))
dev.off()


############ FIGURE 2I AND 2J: RFS BY SIGNATURE AND GZMB LEVELS IN GSE39582 ######################

cols <- c("Proliferative_signature.levels","Stromal_signature.levels","Immune_signature.levels")

for(i in cols[1:2]){
  j <- paste0(gsub("_signature.levels","",i), "_GZMB")
  pData(eset)[,j] <- paste0(gsub("_signature.levels"," signature",i),"=" ,pData(eset)[,i],"; GZMB=", pData(eset)$GZMB.Levels)

# Plot KM curves for each signature 
x <- 'RFS'
m <- rfs

fit <- survfit(m ~ get(j), data = pData(eset))
p <- ggsurvplot(fit, 
   title = paste("GSE39582:",gsub("_signature.levels"," signature",i)),
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
   #legend.title = "EN signature & GZMB",
   legend.labs = levels(factor(pData(eset)[[j]])),
   ggtheme = theme_bw()) # Change ggplot2 theme
p$plot <- p$plot + guides(colour=guide_legend(nrow = 4,ncol=1)) + theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size=12), legend.title=element_text(size=12)) 

if (i=="Proliferative_signature.levels") pdf("../figures/Figure2I_GZMBvsProliferative_GSE39582.pdf", width=6, height=6, onefile=FALSE) else if (i=="Stromal_signature.levels") pdf("../figures/Figure2J_GZMBvsStromal_GSE39582.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()
}


########################### FIGURE 3B: GZMB EXPRESSION BY CMS SUBTYPE ###########################

i <- 'GZMB'
xlabs <- paste(levels(factor(pData(eset)$CMS)),"\n(N=",table(pData(eset)$CMS),")",sep="")
p <- ggplot(data=pData(eset), aes_string(x = 'CMS', y= as.character(i))) +
  geom_boxplot(aes(fill=factor(CMS))) + ylim(-1.5,5) +
  theme_bw() + 
  labs(title=paste0("GSE39582: ",i), x = '', y = "Expression") +
  theme(plot.title = element_text(lineheight=1.5, face="bold", hjust=0.5, size=20), 
        axis.title = element_text(lineheight=1.5, face="bold", hjust=0.5,size=16), 
        axis.text = element_text(lineheight=1.5, face="bold", hjust=0.5, colour = "black", size=14),
        legend.position="none") + 
  scale_x_discrete(labels=xlabs) 

m <- as.data.frame(with(pData(eset), pairwise.t.test(get(i), CMS))[["p.value"]]) %>% signif(.,2)
m[is.na(m)] <- "-"

pdf(paste0("../figures/Figure3B_",i,"_byCMS_GSE39582.pdf"), width=6, height=6, onefile=FALSE)
print(p + annotation_custom(tableGrob(m), xmin=3, xmax=3.5, ymin=3.5, ymax=4.5))
dev.off()


########################### FIGURE 3E: T-EFFECTOR SIGNATURE BY CMS SUBTYPE ###########################

i <- 'Teffector_signature'
j <- "T effector signature"
xlabs <- paste(levels(factor(pData(eset)$CMS)),"\n(N=",table(pData(eset)$CMS),")",sep="")
p <- ggplot(data=pData(eset), aes_string(x = 'CMS', y= as.character(i))) +
  geom_boxplot(aes(fill=factor(CMS))) + ylim(-1.5,5) +
  theme_bw() + 
  labs(title=paste0("GSE39582: ",j), x = '', y = "Expression") +
  theme(plot.title = element_text(lineheight=1.5, face="bold", hjust=0.5, size=20), 
        axis.title = element_text(lineheight=1.5, face="bold", hjust=0.5,size=16), 
        axis.text = element_text(lineheight=1.5, face="bold", hjust=0.5, colour = "black", size=14),
        legend.position="none") + 
  scale_x_discrete(labels=xlabs) 

m <- as.data.frame(with(pData(eset), pairwise.t.test(get(i), CMS))[["p.value"]]) %>% signif(.,2)
m[is.na(m)] <- "-"

pdf(paste0("../figures/Figure3E_",i,"_byCMS_GSE39582.pdf"), width=6, height=6, onefile=FALSE)
print(p + annotation_custom(tableGrob(m), xmin=3, xmax=3.5, ymin=3.5, ymax=4.5))
dev.off()


########################### FIGURE 3H: CORRELATION OF GZMA/B AND T-EFFECTOR ###########################

# Teffector signature vs GZMA/GZMB in CMS
pvalues <-ddply(pData(eset), "CMS", summarise, corr=cor.test(GZMB, Teffector_signature_modified, method="spearman")$p.val)
d <-ddply(pData(eset), "CMS", summarise, corr=cor(GZMB, Teffector_signature_modified, method="spearman"))
names(d) <-c('CMS_subtypes','GZMB')
eset$GZMA[which(eset$GZMA>21)] <- NA
d$GZMA <- ddply(pData(eset), "CMS", summarise, corr=cor(GZMA, Teffector_signature_modified, use="pairwise.complete.obs"))[,2]
rownames(d) <- d$CMS_subtypes
d <- d[,c(1,3,2)]
col1 <- colorRampPalette(c("red","white","blue"))

pdf("../figures/Figure3H_GZMABvsModifiedTeff_GSE39582.pdf", width=6, height=6, onefile=FALSE)
corrplot(t(d[c(2,3)]), is.corr=FALSE, cl.lim=c(0,1), col=col1(20), title = NA,mar=c(0,0,1,0), tl.cex=1.4)
dev.off()


########################### FIGURE 3J: NK SIGNATURE BY CMS SUBTYPE ###########################

i <- 'NK_signature'
j <- "NK signature"
xlabs <- paste(levels(factor(pData(eset)$CMS)),"\n(N=",table(pData(eset)$CMS),")",sep="")
p <- ggplot(data=pData(eset), aes_string(x = 'CMS', y= as.character(i))) +
  geom_boxplot(aes(fill=factor(CMS))) + ylim(-1.5,5) +
  theme_bw() + 
  labs(title=paste0("GSE39582: ",j), x = '', y = "Expression") +
  theme(plot.title = element_text(lineheight=1.5, face="bold", hjust=0.5, size=20), 
        axis.title = element_text(lineheight=1.5, face="bold", hjust=0.5,size=16), 
        axis.text = element_text(lineheight=1.5, face="bold", hjust=0.5, colour = "black", size=14),
        legend.position="none") + 
  scale_x_discrete(labels=xlabs) 

m <- as.data.frame(with(pData(eset), pairwise.t.test(get(i), CMS))[["p.value"]]) %>% signif(.,2)
m[is.na(m)] <- "-"

pdf(paste0("../figures/Figure3J_",i,"_byCMS_GSE39582.pdf"), width=6, height=6, onefile=FALSE)
print(p + annotation_custom(tableGrob(m), xmin=3, xmax=3.5, ymin=3.5, ymax=4.5))
dev.off()


########################### FIGURE 4E: pDC SIGNATURE BY CMS SUBTYPE ###########################

i <- 'pDC_CIWG'
j <- 'pDC signature'
xlabs <- paste(levels(factor(pData(eset)$CMS)),"\n(N=",table(pData(eset)$CMS),")",sep="")
p <- ggplot(data=pData(eset), aes_string(x = 'CMS', y= as.character(i))) +
  geom_boxplot(aes(fill=factor(CMS))) + ylim(-1.5,5) +
  theme_bw() + 
  labs(title=paste0("GSE39582: ",j), x = '', y = "Expression") +
  theme(plot.title = element_text(lineheight=1.5, face="bold", hjust=0.5, size=20), 
        axis.title = element_text(lineheight=1.5, face="bold", hjust=0.5,size=16), 
        axis.text = element_text(lineheight=1.5, face="bold", hjust=0.5, colour = "black", size=14),
        legend.position="none") + 
  scale_x_discrete(labels=xlabs) 

m <- as.data.frame(with(pData(eset), pairwise.t.test(get(i), CMS))[["p.value"]]) %>% signif(.,2)
m[is.na(m)] <- "-"

pdf(paste0("../figures/Figure4E_",j,"_byCMS_GSE39582.pdf"), width=6, height=6, onefile=FALSE)
print(p + annotation_custom(tableGrob(m), xmin=3, xmax=3.5, ymin=3.5, ymax=4.5))
dev.off()


####### FIGURE S3C: PREVALENCE OF KRAS MUTATION BY CMS SUBTYPE IN GSE39582 ########

table(eset$CMS, eset$kras.mutation)
chisq.test(eset$CMS, eset$kras.mutation) # 2.74e-10
KRAS <- factor(eset$kras.mutation, levels=c("M", "WT"))

m <- 100*prop.table(xtabs(~ KRAS + eset$CMS), 2)

pdf("../figures/FigureS3C_GSE39582_KRASvsCMS.pdf", width=6, height=6, onefile=FALSE)
par(mar=c(5.1, 5.1, 4.1, 7.1), xpd=TRUE)
barplot(m, col=c('red','grey'), main='KRAS mutation: GSE39582', cex.names = 1.5, cex.axis=1.5, ylab='Percentage', width=2, cex.lab=1.5, las=2)
legend('topright', inset=c(-0.3,0), rownames(m), fill=c('red','grey'))
dev.off()


####### FIGURE S3F: PREVALENCE OF BRAF MUTATION BY CMS SUBTYPE IN GSE39582 ########

table(eset$CMS, eset$braf.mutation)
chisq.test(eset$CMS, eset$braf.mutation) # p <2.2e-16

BRAF <- factor(eset$braf.mutation, levels=c("M", "WT"))

m <- 100*prop.table(xtabs(~ BRAF + eset$CMS), 2)

pdf("../figures/FigureS3F_GSE39582_BRAFvsCMS.pdf", width=6, height=6, onefile=FALSE)
par(mar=c(5.1, 5.1, 4.1, 7.1), xpd=TRUE)
barplot(m, col=c('red','grey'), main='BRAF mutation: GSE39582', cex.names = 1.5, cex.axis=1.5, ylab='Percentage', width=2, cex.lab=1.5, las=2)
legend('topright', inset=c(-0.3,0), rownames(m), fill=c('red','grey'))
dev.off()


############## FIGURE S7A: RECURRENCE-FREE SURVIVAL BY SIGNATURE QUANTILE (WEIGHTED) #############

x <- 'RFS'
m <- rfs
fit<- survfit(rfs ~ glmnetW.quartiles, data = pData(eset))
# Drawing survival curves
p <- ggsurvplot(fit, 
   title = paste("Validation dataset : GSE39582"),
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
   legend.title = "Prediction score Quartiles",
   legend.labs = c('lowest', 'low','high','highest'),
   ggtheme = theme_bw()) # Change ggplot2 theme
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=12))

pdf("../figures/FigureS7A_AVANTsignature_GSE39582_WeightedScores.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()


############## FIGURE S7B: FOREST PLOT FOR RFS WITH SIGNATURE QUANTILE (WEIGHTED) #############

x <- 'RFS'
m <- rfs
    myforest <- data.frame(matrix(ncol=0,nrow=0))
  cols <- c("glmnet.weighted",signatures.names[-c(7,8)]) %>% paste0(., "_std")
  for(i in cols){
  # Individual gene signatures prognostic value 
  temp <- with(pData(eset), summary(coxph(m ~ get(i)))) 
  myforest[i,"HR"] <- temp$coefficients[[2]]
  myforest[i,"CI_lower"] <- temp$conf.int[[3]]
  myforest[i,"CI_upper"] <- temp$conf.int[[4]]
  myforest[i,"Cox_Pvalue"] <- signif(temp$coefficients[[5]],2)
  }
  cols <- signatures.names[-c(7,8)] %>% paste0(., "_std")
  for(i in cols){
  # Added value of individual signatures to elastic net
   myanova <- anova(
    with(pData(eset), coxph(m ~ glmnet.weighted_std)), 
    with(pData(eset), coxph(m ~ glmnet.weighted_std + get(i)))
    )
  myforest[i,"Published signature added benefit \n to AVANT signature \n (P-value)"] <- signif(myanova$`P(>|Chi|)`[2],2)
  # Added value of elastic net to individual signatures
    myanova <- anova(
    with(pData(eset), coxph(m ~ get(i))), 
    with(pData(eset), coxph(m ~ get(i) + glmnet.weighted_std))
    )
  myforest[i,"AVANT signature added benefit \n to published signatures \n (P-value)"] <- signif(myanova$`P(>|Chi|)`[2],2)
  }
  rownames(myforest) <- gsub("_std", "", rownames(myforest))
  myforest$Group <- rownames(myforest)
  myforest <- myforest[-c(4,5),]
  myforest$Group <- c("AVANT signature", "OncotypeDx", "Teffector signature", "CAF Isella et al", "F-TBRS Calon et al")
  myforest2 <- data.frame(NA, NA, NA, NA, 0.016, 4.7e-4, "CMS subtypes")
  colnames(myforest2) <- colnames(myforest)
  myforest2 <- rbind(myforest, myforest2)
  rownames(myforest2)[nrow(myforest2)] <- "CMS_subtypes"
# Forest plot with p-values
p <- ggplot(myforest, aes(x = Group, y = HR, ymin = CI_lower, ymax = CI_upper), colour = "black") +
            geom_hline(aes(yintercept = 1), colour = 'grey', linetype = 2) +
            geom_pointrange(size = 0.7, shape = 15) +
            scale_y_log10(limits = c(0.3, 3), breaks = c(0.3, 1, 3)) +
            scale_x_discrete(limits = rev(myforest$Group)) + theme_bw() + 
            ggtitle(paste(x, "\n Validation Dataset: GSE39582")) + 
            coord_flip() + 
            theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5), axis.text = element_text(size=12,face="bold"))  + annotate("text", x = 5:1, y=2.5, label=paste("P=", myforest$Cox_Pvalue))
tt <- gridExtra::ttheme_minimal(colhead=list(fg_params = list(parse=TRUE)))
tbl <- gridExtra::tableGrob(myforest2[-1,c(7,5,6)], rows=NULL, theme=tt)
# Plot chart and table into one object

pdf("../figures/FigureS7B_ForestPlot_GSE39582_WeightedScores_OncotypeDx7genes_WEIGHTED.pdf", width=7.2, height=6.5)
gridExtra::grid.arrange(p, tbl,
             ncol=1,nrow=2,
             as.table=TRUE,
             heights=c(5,3))
dev.off()


##################### FIGURE S10: RFS BY AVANT SUB-SIGNATURES #####################

for(i in names(mygenes)[c(1,4,5)]){
j <- mygenes[i]
genes <- mygenes[[i]]
k <- paste0(names(j),".levels")
x <- 'RFS'
m <- rfs
fit <- survfit(m ~ pData(eset)[,k], data = pData(eset))
# Drawing survival curves
p <- ggsurvplot(fit, 
   title = paste(names(j), "\n Validation dataset : GSE39582"),
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

if (i=="Proliferative_signature") pdf(paste0("../figures/FigureS10A_",i,"_GSE39582.pdf"), width=6, height=6, onefile=FALSE) else if (i=="Immune_signature") pdf(paste0("../figures/FigureS10B_",i,"_GSE39582.pdf"), width=6, height=6, onefile=FALSE) else if (i=="Stromal_signature") pdf(paste0("../figures/FigureS10C_",i,"_GSE39582.pdf"), width=6, height=6, onefile=FALSE)
print(p)
dev.off()
}


######## FIGURE S12B RESULTS ########

x <- 'RFS'
m <- rfs
print(with(pData(eset), coxph(m ~ glmnet.unweighted + age.at.diagnosis + tnm.stage + tnm.n + Sex)))


######################## FIGURE S13: RFS IN CMS2, T-EFFECTOR LOW TUMORS ######################

### Subset of T effector signature low CMS2 tumors
x <- 'RFS'
i <- "CMS2"
dsub <- subset(pData(eset), CMS==i)
dsub <- subset(dsub, Teffector_signature_noGZMB.MeanLevels=="low")
rfs <- with( 
  dsub,
  Surv(
    time  = rfs.delay,
    event = rfs.event
  )
)

# Make survival curves 
m <- rfs
fit<- survfit(m ~ GZMB.Levels, data = dsub)
p <- ggsurvplot(fit, 
	title = paste0("GZMB\nValidation set: GSE39582 (",i,", n=",nrow(dsub),")"),
   font.main = c(14, "bold", "black"),
   font.x = c(14, "bold", "black"),
   font.y = c(14, "bold", "black"),
   font.tickslab = 14, 
   palette = c(1:2),
   xlab = "Time (months)", 
   ylab = x,
   conf.int = FALSE, 
   pval=TRUE, 
   # linetype = "strata", change line type based on groups - solid/dashed
   #linetype = c(1,2),
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
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5, size=20),legend.text=element_text(size=14), legend.title=element_text(size=14)) #+ scale_color_manual(values=c(CMS.cols[[i]],CMS.cols[[i]])) #+ scale_alpha_discrete(range=c(0.4, 1))

pdf("../figures/FigureS13_GZMBexpression_CMS2andTeffLow_GSE39582.pdf", width=6, height=6, onefile=FALSE)
print(p)
dev.off()

