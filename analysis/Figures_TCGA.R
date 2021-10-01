source("Load_Libraries.R")
source("Functions.R")

####### CODE FOR FIGURES 4F, 4G
####### CODE FOR SUPPLEMENTARY FIGURE S3B, S3E, S3H, S14D


### Load TCGA data: expression, annotation
load("../data/TCGA_GRCh38.RData")

### TCGA CMS subtypes
tcga.pData <- read.csv("../data/TCGA_CMS_Labels.csv", sep=",", header=TRUE)
table(tcga.pData$RFscore)
tcga.pData$RFscore <- as.character(tcga.pData$RFscore)
tcga.pData$RFscore[tcga.pData$RFscore=='CMS1,CMS2'| tcga.pData$RFscore=='CMS1,CMS3'| tcga.pData$RFscore=='CMS1,CMS4'] <- 'unclassified'   # Rename 

col1 <- colorRampPalette(c("red","white","blue"))


##### FIGURE 4F: T-EFFECTOR SIGNATURE CORRELATION IN TCGA #####

goi <- c("GZMA", "CXCL9", "CD8A", "PRF1", "GZMB", "IFNG", "TBX21", "CXCL10")
data.goi <- log2(exprs(eset)[match(goi, fData(eset)$symbol), ]+1)
rownames(data.goi) <- fData(eset)[rownames(data.goi), "symbol"]

### Tissues of interest
toi <- c("Colon", "Breast", "Lung", "Head_Neck", "Bladder", "Ovary", "Kidney", "Skin", "Thyroid", "Cervix", "Liver", "Prostate", "Stomach", "Uterus")
data.goi <- data.goi[, which(eset$Tissue %in% toi)]
tissues <- as.character(eset$Tissue[which(eset$Tissue %in% toi)])

### Calculate correlation of each gene with the remaining genes, per indication
data.scaled <- t(scale(t(data.goi), center=TRUE, scale=TRUE))

avg.no <- vector("list", length(goi))
names(avg.no) <- goi
for (g in goi) {
	avg.no[[g]] <- colMeans(data.scaled[-match(g, rownames(data.scaled)), ])
}

cors.breast <- cor.values("Breast")
cors.bladder <- cor.values("Bladder")
cors.colon <- cor.values("Colon")
cors.head_neck <- cor.values("Head_Neck")
cors.lung <- cor.values("Lung")
cors.ovary <- cor.values("Ovary")
cors.kidney <- cor.values("Kidney")
cors.skin <- cor.values("Skin")
cors.thyroid <- cor.values("Thyroid")
cors.cervix <- cor.values("Cervix")
cors.liver <- cor.values("Liver")
cors.prostate <- cor.values("Prostate")
cors.stomach <- cor.values("Stomach")
cors.uterus <- cor.values("Uterus")

### Calculate correlation of each gene with the remaining genes, per CMS subtype in COAD
common <- intersect(gsub("_","-",as.character(tcga.pData$X)),
									substr(colnames(data.goi),1,16))
indices.coad <- match(common, substr(colnames(data.goi),1,16))
CMS <- rep(NA, length(tissues))
CMS[indices.coad] <- tcga.pData[match(common, gsub("_","-",as.character(tcga.pData$X))), "RFscore"]
CMS[CMS %in% "unclassified"] <- NA

### Calculate correlation of each gene with the remaining genes, per CMS subtype
cors.colon_CMS1 <- cor.values.bysubtype("Colon", "CMS1")
cors.colon_CMS2 <- cor.values.bysubtype("Colon", "CMS2")

cors.values <- data.frame(Breast=cors.breast, Bladder=cors.bladder, CRC=cors.colon,
						`Head & Neck`=cors.head_neck, Lung=cors.lung, Ovary=cors.ovary,
						`CRC CMS1`=cors.colon_CMS1, `CRC CMS2`=cors.colon_CMS2, 
						Kidney=cors.kidney, Skin=cors.skin, Thyroid=cors.thyroid, Cervix=cors.cervix,
						Liver=cors.liver, Prostate=cors.prostate, Stomach=cors.stomach, Uterus=cors.uterus)
colnames(cors.values) <- c("Breast" ,"Bladder", "CRC", "Head & Neck", "Lung", "Ovary", "CRC CMS1", "CRC CMS2",
											"Kidney", "Skin", "Thyroid", "Cervix", "Liver", "Prostate", "Stomach", "Uterus")
cors.values <- cors.values[,c(3,7,8,5,6,1,4,2,9:16)]
cors.values <- cors.values[c(5,1:4,6:8), ]

pdf("../figures/Figure4F_TCGAfullpancancer_GZMBvsTeffector_CorrelationHeatmap_CMSsubtype1and2.pdf")#, width=6, height=5)
corrplot(as.matrix(cors.values), is.corr=FALSE, cl.lim=c(0,1), col=col1(20), title = NA,mar=c(0,0,1,0), tl.col="black")
dev.off()


##### FIGURE 4G: NK SIGNATURE CORRELATION IN TCGA #####

goi <- c("KLRC3", "KLRK1", "NKG7", "KLRC2", "KLRD1", "GZMB", "FCGR3A") 
data.goi <- log2(exprs(eset)[match(goi, fData(eset)$symbol), ]+1)
rownames(data.goi) <- fData(eset)[rownames(data.goi), "symbol"]
goi <- rownames(data.goi)

### Tissues of interest
toi <- c("Colon", "Breast", "Lung", "Head_Neck", "Bladder", "Ovary", "Kidney", "Skin", "Thyroid", "Cervix", "Liver", "Prostate", "Stomach", "Uterus")
data.goi <- data.goi[, which(eset$Tissue %in% toi)]
tissues <- as.character(eset$Tissue[which(eset$Tissue %in% toi)])

### Calculate correlation of each gene with the remaining genes, per indication
data.scaled <- t(scale(t(data.goi), center=TRUE, scale=TRUE))

avg.no <- vector("list", length(goi))
names(avg.no) <- goi
for (g in goi) {
	avg.no[[g]] <- colMeans(data.scaled[-match(g, rownames(data.scaled)), ])
}

cors.breast <- cor.values("Breast")
cors.bladder <- cor.values("Bladder")
cors.colon <- cor.values("Colon")
cors.head_neck <- cor.values("Head_Neck")
cors.lung <- cor.values("Lung")
cors.ovary <- cor.values("Ovary")
cors.kidney <- cor.values("Kidney")
cors.skin <- cor.values("Skin")
cors.thyroid <- cor.values("Thyroid")
cors.cervix <- cor.values("Cervix")
cors.liver <- cor.values("Liver")
cors.prostate <- cor.values("Prostate")
cors.stomach <- cor.values("Stomach")
cors.uterus <- cor.values("Uterus")

### Calculate correlation of each gene with the remaining genes, per CMS subtype in COAD
common <- intersect(gsub("_","-",as.character(tcga.pData$X)),
									substr(colnames(data.goi),1,16))
indices.coad <- match(common, substr(colnames(data.goi),1,16))
CMS <- rep(NA, length(tissues))
CMS[indices.coad] <- tcga.pData[match(common, gsub("_","-",as.character(tcga.pData$X))), "RFscore"]
CMS[CMS %in% "unclassified"] <- NA

### Calculate correlation of each gene with the remaining genes, per CMS subtype
cors.colon_CMS1 <- cor.values.bysubtype("Colon", "CMS1")
cors.colon_CMS2 <- cor.values.bysubtype("Colon", "CMS2")

cors.values <- data.frame(Breast=cors.breast, Bladder=cors.bladder, CRC=cors.colon,
						`Head & Neck`=cors.head_neck, Lung=cors.lung, Ovary=cors.ovary,
						`CRC CMS1`=cors.colon_CMS1, `CRC CMS2`=cors.colon_CMS2, 
						Kidney=cors.kidney, Skin=cors.skin, Thyroid=cors.thyroid, Cervix=cors.cervix,
						Liver=cors.liver, Prostate=cors.prostate, Stomach=cors.stomach, Uterus=cors.uterus)
colnames(cors.values) <- c("Breast" ,"Bladder", "CRC", "Head & Neck", "Lung", "Ovary", "CRC CMS1", "CRC CMS2",
											"Kidney", "Skin", "Thyroid", "Cervix", "Liver", "Prostate", "Stomach", "Uterus")
cors.values <- cors.values[,c(3,7,8,5,6,1,4,2,9:16)]
cors.values <- cors.values[c(6,2,4,1,5,3,7), ]
cors.values[cors.values<0] <- 0

pdf("../figures/Figure4G_TCGAfullpancancer_GZMBvsNKsignature_CorrelationHeatmap_CMSsubtype1and2.pdf")#, width=6, height=5)
corrplot(as.matrix(cors.values), is.corr=FALSE, cl.lim=c(0,1), col=col1(20), title = NA,mar=c(0,0,1,0), tl.col="black")
dev.off()


##### FIGURE S3B: PREVALENCE OF KRAS MUTATION BY CMS SUBTYPE IN TCGA COLON #####

tcga <- get(load("../data/TCGA.expression.rdata"))
tumor.samples <- which(pData(tcga)$Cancer_Status == "tumor")
tcga <- tcga[,tumor.samples]
# Downloaded KRAS and BRAF mutations from cBioPortal, July 18, 2018
mut.KRAS <- read.table("../data/TCGA CRC cbioportal KRAS mutations.txt")
mut.BRAF <- read.table("../data/TCGA CRC cbioportal BRAF mutations.txt")
mutations.KRAS <- mut.KRAS[,2]
names(mutations.KRAS) <- gsub("-", "_", sapply(strsplit(as.character(mut.KRAS[,1]), ":"), "[[", 2))
mutations.BRAF <- mut.BRAF[,2]
names(mutations.BRAF) <- gsub("-", "_", sapply(strsplit(as.character(mut.BRAF[,1]), ":"), "[[", 2))
table(mutations.KRAS) # 96 MUT, 127 WT
table(mutations.BRAF) # 22 MUT, 201 WT
sampleNames(tcga) <- substr(sampleNames(tcga), 1, 15)

### Add downloaded information to pData(tcga)
pData(tcga)$KRAS.cBioPortal <- NA
pData(tcga)$BRAF.cBioPortal <- NA
mutations.KRAS[mutations.KRAS==1] <- "MUT"
mutations.KRAS[mutations.KRAS=="0"] <- "WT"
mutations.BRAF[mutations.BRAF==1] <- "MUT"
mutations.BRAF[mutations.BRAF=="0"] <- "WT"

common <- intersect(sampleNames(tcga), names(mutations.KRAS))
pData(tcga)[common, "KRAS.cBioPortal"] <- mutations.KRAS[common] # 54 MUT, 93 WT, 266 unknown
pData(tcga)[common, "BRAF.cBioPortal"] <- mutations.BRAF[common] # 20 MUT, 127 WT, 266 unknown

exprs.tcga <-  t(scale(t(exprs(tcga)), center = TRUE, scale = TRUE))
rownames(exprs.tcga) <- fData(tcga)$EntrezID

# Runs the classifier on the dataset and provides CMS labels on the TCGA samples
Rfcms <- CMSclassifier::classifyCMS(exprs.tcga,method="RF")[[3]] # Genes are rows, samples are columns. Rownames need to be entrezIDs. 
pData(tcga)$RFscore <- Rfcms$RF.nearestCMS
pData(tcga)$RFscore[pData(tcga)$RFscore=='CMS1,CMS2'|pData(tcga)$RFscore=='CMS1,CMS3'|pData(tcga)$RFscore=='CMS1,CMS4'] <- 'unclassified'   # Rename unclassified samples or having more than 1 class as unclassified

pData(tcga)$microsatellite_instability[pData(tcga)$microsatellite_instability=='MSI-L'] <- 'MSS'
# Remove factors from MSI
pData(tcga)$microsatellite_instability <- sapply(pData(tcga)$microsatellite_instability, as.character, stringsAsFactors=FALSE)

# Percentage of KRAS mutation in CMS subtypes
m <- 100*prop.table(xtabs(~ KRAS.cBioPortal + RFscore, data = pData(tcga)[pData(tcga)$RFscore !='unclassified',]), 2)

pdf("../figures/FigureS3B_TCGA_KRASvsCMS.pdf", width=6, height=6, onefile=FALSE)
par(mar=c(5.1, 5.1, 4.1, 7.1), xpd=TRUE)
barplot(m, col=c('red','grey'), main='KRAS mutation: TCGA', cex.names = 1.5, cex.axis=1.5, ylab='Percentage', width=2, cex.lab=1.5, las=2)
legend('topright', inset=c(-0.3,0), rownames(m), fill=c('red','grey'))
dev.off()

table(pData(tcga)$RFscore[pData(tcga)$RFscore !='unclassified'], pData(tcga)$KRAS.cBioPortal[pData(tcga)$RFscore !='unclassified'])
chisq.test(pData(tcga)$RFscore[pData(tcga)$RFscore !='unclassified'], pData(tcga)$KRAS.cBioPortal[pData(tcga)$RFscore !='unclassified'])
# p 0.01182


##### FIGURE S3E: PREVALENCE OF BRAF MUTATION BY CMS SUBTYPE IN TCGA COLON #####

# Percentage of BRAF mutation in CMS subtypes
m <- 100*prop.table(xtabs(~ BRAF.cBioPortal + RFscore, data = pData(tcga)[pData(tcga)$RFscore !='unclassified',]), 2)

pdf("../figures/FigureS3E_TCGA_BRAFvsCMS.pdf", width=6, height=6, onefile=FALSE)
par(mar=c(5.1, 5.1, 4.1, 7.1), xpd=TRUE)
barplot(m, col=c('red','grey'), main='BRAF mutation: TCGA', cex.names = 1.5, cex.axis=1.5, ylab='Percentage', width=2, cex.lab=1.5, las=2)
legend('topright', inset=c(-0.3,0), rownames(m), fill=c('red','grey'))
dev.off()

table(pData(tcga)$RFscore[pData(tcga)$RFscore !='unclassified'], pData(tcga)$BRAF.cBioPortal[pData(tcga)$RFscore !='unclassified'])
chisq.test(pData(tcga)$RFscore[pData(tcga)$RFscore !='unclassified'], pData(tcga)$BRAF.cBioPortal[pData(tcga)$RFscore !='unclassified'])
# p 8.212e-10


##### FIGURE S3H: PREVALENCE OF MSI STATUS BY CMS SUBTYPE IN TCGA COLON #####

# Percentage of MSI/MSS distribution in CMS subtypes
m <- 100*prop.table(xtabs(~ microsatellite_instability + RFscore, data = pData(tcga)[pData(tcga)$RFscore !='unclassified',]), 2)

pdf("../figures/FigureS3H_TCGA_MSIvsCMS.pdf", width=6, height=6, onefile=FALSE)
par(mar=c(5.1, 5.1, 4.1, 7.1), xpd=TRUE)
barplot(m, col=c('red','grey'), main='MSI status: TCGA', cex.names = 1.5, cex.axis=1.5, ylab='Percentage', width=2, cex.lab=1.5, las=2)
legend('topright', inset=c(-0.3,0), rownames(m), fill=c('red','grey'))
dev.off()

table(pData(tcga)$RFscore[pData(tcga)$RFscore !='unclassified'], pData(tcga)$microsatellite_instability[pData(tcga)$RFscore !='unclassified'])
chisq.test(pData(tcga)$RFscore[pData(tcga)$RFscore !='unclassified'], pData(tcga)$microsatellite_instability[pData(tcga)$RFscore !='unclassified'])
# p <2.2e-16


##### FIGURE S14D: pDC SIGNATURE CORRELATION IN TCGA #####

goi <- c("GZMB","CLEC4C", "IL3RA", "TCF4", "NRP1") 
data.goi <- log2(exprs(eset)[match(goi, fData(eset)$symbol), ]+1)
rownames(data.goi) <- fData(eset)[rownames(data.goi), "symbol"]
goi <- rownames(data.goi)

### Tissues of interest
toi <- c("Colon", "Breast", "Lung", "Head_Neck", "Bladder", "Ovary", "Kidney", "Skin", "Thyroid", "Cervix", "Liver", "Prostate", "Stomach", "Uterus")
data.goi <- data.goi[, which(eset$Tissue %in% toi)]
tissues <- as.character(eset$Tissue[which(eset$Tissue %in% toi)])

### Calculate correlation of each gene with the remaining genes, per indication
data.scaled <- t(scale(t(data.goi), center=TRUE, scale=TRUE))

avg.no <- vector("list", length(goi))
names(avg.no) <- goi
for (g in goi) {
	avg.no[[g]] <- colMeans(data.scaled[-match(g, rownames(data.scaled)), ])
}

cors.breast <- cor.values("Breast")
cors.bladder <- cor.values("Bladder")
cors.colon <- cor.values("Colon")
cors.head_neck <- cor.values("Head_Neck")
cors.lung <- cor.values("Lung")
cors.ovary <- cor.values("Ovary")
cors.kidney <- cor.values("Kidney")
cors.skin <- cor.values("Skin")
cors.thyroid <- cor.values("Thyroid")
cors.cervix <- cor.values("Cervix")
cors.liver <- cor.values("Liver")
cors.prostate <- cor.values("Prostate")
cors.stomach <- cor.values("Stomach")
cors.uterus <- cor.values("Uterus")

### Calculate correlation of each gene with the remaining genes, per CMS subtype in COAD
common <- intersect(gsub("_","-",as.character(tcga.pData$X)),
									substr(colnames(data.goi),1,16))
indices.coad <- match(common, substr(colnames(data.goi),1,16))
CMS <- rep(NA, length(tissues))
CMS[indices.coad] <- tcga.pData[match(common, gsub("_","-",as.character(tcga.pData$X))), "RFscore"]
CMS[CMS %in% "unclassified"] <- NA

### Calculate correlation of each gene with the remaining genes, per CMS subtype
cors.colon_CMS1 <- cor.values.bysubtype("Colon", "CMS1")
cors.colon_CMS2 <- cor.values.bysubtype("Colon", "CMS2")

cors.values <- data.frame(Breast=cors.breast, Bladder=cors.bladder, CRC=cors.colon,
						`Head & Neck`=cors.head_neck, Lung=cors.lung, Ovary=cors.ovary,
						`CRC CMS1`=cors.colon_CMS1, `CRC CMS2`=cors.colon_CMS2, 
						Kidney=cors.kidney, Skin=cors.skin, Thyroid=cors.thyroid, Cervix=cors.cervix,
						Liver=cors.liver, Prostate=cors.prostate, Stomach=cors.stomach, Uterus=cors.uterus)
colnames(cors.values) <- c("Breast" ,"Bladder", "CRC", "Head & Neck", "Lung", "Ovary", "CRC CMS1", "CRC CMS2",
											"Kidney", "Skin", "Thyroid", "Cervix", "Liver", "Prostate", "Stomach", "Uterus")
cors.values <- cors.values[,c(3,7,8,5,6,1,4,2,9:16)]

pdf("../figures/FigureS14D_TCGAfullpancancer_GZMBvspDCsignature_CorrelationHeatmap_CMSsubtype1and2.pdf")#, width=6, height=5)
corrplot(as.matrix(cors.values), is.corr=FALSE, cl.lim=c(0,1), col=col1(20), title = NA,mar=c(0,0,1,0), tl.col="black")
dev.off()

