source("Load_Libraries.R")

####### CODE FOR SUPPLEMENTARY FIGURES S4A-F

### Setup 
get(load("../data/avant_analysis_dataset.rdata"))
names(pData(avant))[names(pData(avant))=="CMSnew"] <- "CMS_Subtypes"

CMS.cols <- setNames(
  c('orange','blue','magenta','darkgreen', 'grey'),
  c('CMS1', 'CMS2', 'CMS3', 'CMS4', 'unclassified')
)
CMS.genes = setNames(
    c('orange','blue','magenta','darkgreen','grey'),
  c('CMS1', 'CMS2', 'CMS3', 'CMS4','NA')
)
Marisa.cols <- setNames(
  c("orange","blue","green","grey","purple","red",  "black"),
  c("C4", "C2", "C3","C6", "C5", "C1", "Unclassified")
)


######### FIG S4B: PAM HEATMAP IN TCGA #########

## load preprocessed TCGA RNAseq data from Thomas
tcga <- get(load("../data/TCGA.expression.rdata"))

# Subset only tumor samples from TCGA
tumor.samples <- which(pData(tcga)$Cancer_Status == "tumor")
tcga <- tcga[,tumor.samples]

exprs.tcga <-  t(scale(t(exprs(tcga)), center = TRUE, scale = TRUE))
rownames(exprs.tcga) <- fData(tcga)$EntrezID

# Runs the classifier on the dataset and provides CMS labels on the TCGA samples
Rfcms <- CMSclassifier::classifyCMS(exprs.tcga,method="RF")[[3]] # Genes are rows, samples are columns. Rownames need to be entrezIDs. 
pData(tcga)$RFscore <- Rfcms$RF.nearestCMS
pData(tcga)$RFscore[pData(tcga)$RFscore=='CMS1,CMS2'|pData(tcga)$RFscore=='CMS1,CMS3'|pData(tcga)$RFscore=='CMS1,CMS4'] <- 'unclassified'   # Rename unclassified samples or having more than 1 class as unclassified

# Subset tcga dataset to only contain the 829 genes from AVANT nanostring
test <- tcga[match(fData(avant)$EntrezID, fData(tcga)$EntrezID),]
test <- test[rowSums(is.na(fData(test)))==0,] # removes rows with NA values
test <- test[,pData(test)$RFscore !="unclassified"]

### Obtain PAM classifier genes from TCGA
a <- read.csv('../data/PAMR_Classifier_genes_TCGA_v2.csv',sep=",", header=TRUE)
b <- a[a$cms != 'NA',]
# Additional filtering based on prediction scores for each gene
c <- b %>% mutate_each(funs(as.numeric(as.character(.))),2:5)
c$Symbol <- b$Symbol
c$cms <- b$cms
c <- c[c$CMS1.score > abs(0.009)|
       c$CMS2.score > abs(0.009)|
       c$CMS3.score > abs(0.009)|
       c$CMS4.score > abs(0.009),]
c <- c[c$cms !="NA",]
test <- test[,order(pData(test)$RFscore)]

m <- exprs(test)[fData(test)$EntrezID %in% c$id,] # Subset to classifier genes 
# Z score matrix
m <- t(scale(t(m), center = TRUE, scale = TRUE))
m[m < -3] <- -3
m[m >  3] <- 3
rownames(m) <- gsub("GeneID:","",rownames(m))
pData(test)$RFscore <- factor(pData(test)$RFscore)

#Column annotation
ha_column = HeatmapAnnotation(
    df = data.frame(
    CMS_Subtypes=pData(test)[colnames(m), "RFscore"]),
    col = list(
    CMS_Subtypes = CMS.cols), 
    show_legend = TRUE,which='column')

# Row order
row_dend <- as.dendrogram(hcluster(m,method='euclidean',link='complete'))

### Optimal order for TCGA, AVANT and GSE39582 heatmaps
subtypes <- b$cms[match(labels(row_dend), b$id)]
order.rows <- c(labels(row_dend)[subtypes %in% "CMS1"], labels(row_dend)[subtypes %in% "CMS2"],
						labels(row_dend)[subtypes %in% "CMS3"], labels(row_dend)[subtypes %in% "CMS4"]) 

ha_row = HeatmapAnnotation(
  df = data.frame(CMS_genes=b$cms[match(order.rows, b$id)]),
  col = list(CMS_genes = CMS.genes),
  show_legend = TRUE,which='row')

ht1 = Heatmap(m[order.rows, ], cluster_rows =FALSE, 
    cluster_columns = FALSE, 
    show_row_names=FALSE, 
    show_column_names = FALSE,
    column_title = "CMS classes: TCGA patients", 
    name='Z-score', 
    column_title_side = "bottom", 
    row_title = 'PAM Classifier genes', 
    show_row_dend = FALSE, 
    top_annotation = ha_column)
    
ht_list = ht1 + ha_row 

pdf("../figures/FigureS4B_TCGA_PAMRheatmap_FixedGeneOrdering.pdf", width=6, height=6, onefile=FALSE)
draw(ht_list,  heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()


######### FIG S4D: CMS SUBTYPE PREVALENCE IN TCGA #########

pdf("../figures/FigureS4D_TCGA_CMSprevalence.pdf", width=6, height=6, onefile=FALSE)
par(mar=c(8, 5.1, 4.1, 7.1), xpd=TRUE)
barplot(table(tcga$RFscore)*100/sum(table(tcga$RFscore)), col=CMS.cols, main='TCGA', cex.names = 1.5, cex.axis=1.5, ylab='Percentage', width=2, cex.lab=1.2, las=2, ylim=c(0,50))
legend('topright', inset=c(-0.3,0), c("CMS1", "CMS2", "CMS3", "CMS4", "unclassified"), fill=CMS.cols)
dev.off()


######### FIG S4E: PAM HEATMAP IN GSE39582 #########

## load preprocessed microarray data for GSE59382
eset <- get(load("../data/GSE39582.rdata"))
pData(eset)$cit.molecularsubtype <- factor(pData(eset)$cit.molecularsubtype)
#Z-score transformation
exprs(eset) <- t(scale(t(exprs(eset)), center=TRUE, scale=TRUE))

### CMS subtypes
SScms <- classifyCMS(exprs(eset),method="SSP")[[3]] # Genes are rows, samples are columns
pData(eset)$CMS <- SScms$SSP.nearestCMS

### Heatmap with PAM classifier genes
m <- exprs(eset)[match(order.rows,rownames(fData(eset))),] # rownames are entrez ids in fdata
m <- na.omit(m)
m[m < -3] <- -3
m[m >  3] <- 3
m <- m[,order(pData(eset)[colnames(m),'CMS'])]

#Column annotation
ha_column = HeatmapAnnotation(
    df = data.frame(
    Marisa_subtypes = pData(eset)[colnames(m),"cit.molecularsubtype"], 
    CMS_subtypes = pData(eset)[colnames(m),'CMS']),
    col = list( 
    Marisa_subtypes = Marisa.cols, 
    CMS_subtypes = CMS.cols), 
    show_legend = TRUE,which='column') # , show_annotation_name = TRUE, annotation_name_side = "right"

# Row annotation
ha_row = HeatmapAnnotation(
  df = data.frame(CMS_genes=b$cms[match(rownames(m), b$id)]),
  col = list(CMS_genes = CMS.genes),
  show_legend = TRUE,which='row')

# Heatmap using predefined clustering of patients
ht1 = Heatmap(m, cluster_rows =FALSE, 
    cluster_columns = FALSE, 
    show_row_names=FALSE, 
    show_column_names = FALSE,
    column_title = "CMS classes: GSE39582 patients", 
    name='Z-score', 
    column_title_side = "bottom", 
    row_title = 'PAM Classifier genes', 
    show_row_dend = FALSE, 
    top_annotation = ha_column)
    
ht_list = ht1 + ha_row 

pdf("../figures/FigureS4E_GSE39582_PAMRheatmap_FixedGeneOrdering.pdf", width=6, height=6, onefile=FALSE)
draw(ht_list,  heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()


######### FIG S4F: CMS SUBTYPE PREVALENCE IN GSE39582 #########

pdf("../figures/FigureS4F_GSE39582_CMSprevalence.pdf", width=6, height=6, onefile=FALSE)
par(mar=c(8, 5.1, 4.1, 7.1), xpd=TRUE)
barplot(table(eset$CMS)*100/(sum(table(eset$CMS))), col=CMS.cols, main='GSE39582', cex.names = 1.5, cex.axis=1.5, ylab='Percentage', width=2, cex.lab=1.2, las=2, ylim=c(0,50))
legend('topright', inset=c(-0.3,0), c("CMS1", "CMS2", "CMS3", "CMS4"), fill=CMS.cols)
dev.off()


######### FIG S4A: PAM HEATMAP IN AVANT #########

### Heatmap with PAM classifier genes
m <- exprs(avant)[match(order.rows,fData(avant)$EntrezID),] # rownames are entrez ids in fdata
dim(m)
m <- na.omit(m)
m[m < -3] <- -3
m[m >  3] <- 3
m <- m[,order(pData(avant)[colnames(m),'CMS_Subtypes'])]

#Column annotation
ha_column = HeatmapAnnotation(
    df = data.frame(
    CMS_subtypes = pData(avant)[colnames(m),'CMS_Subtypes']),
    col = list( 
    CMS_subtypes = CMS.cols), 
    show_legend = TRUE,which='column') # , show_annotation_name = TRUE, annotation_name_side = "right"

# Row annotation
ha_row = HeatmapAnnotation(
  df = data.frame(CMS_genes=b$cms[match(rownames(m), b$id)]),
  col = list(CMS_genes = CMS.genes),
  show_legend = TRUE,which='row')

# Heatmap using predefined clustering of patients
ht1 = Heatmap(m, cluster_rows =FALSE, 
    cluster_columns = FALSE, 
    show_row_names=FALSE, 
    show_column_names = FALSE,
    column_title = "CMS classes: AVANT patients", 
    name='Z-score', 
    column_title_side = "bottom", 
    row_title = 'PAM Classifier genes', 
    show_row_dend = FALSE, 
    top_annotation = ha_column)
    
ht_list = ht1 + ha_row 

pdf("../figures/FigureS4A_AVANT_PAMRheatmap_FixedGeneOrdering.pdf", width=6, height=6, onefile=FALSE)
draw(ht_list,  heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()


######### FIG S4C: CMS SUBTYPE PREVALENCE IN AVANT #########

pdf("../figures/FigureS4C_AVANT_CMSprevalence.pdf", width=6, height=6, onefile=FALSE)
par(mar=c(8, 5.1, 4.1, 7.1), xpd=TRUE)
barplot(table(avant$CMS_Subtypes)*100/(sum(table(avant$CMS_Subtypes))), col=CMS.cols, main='AVANT', cex.names = 1.5, cex.axis=1.5, ylab='Percentage', width=2, cex.lab=1.2, las=2, ylim=c(0,50))
legend('topright', inset=c(-0.3,0), c("CMS1", "CMS2", "CMS3", "CMS4", "unclassified"), fill=CMS.cols)
dev.off()
