source("Load_Libraries.R")
source("Load_Signatures.R")


####### DATA PREPARATION FOR AVANT ########

### load AVANT nanostring data
get(load("../data/avant_analysis_dataset.rdata"))
names(pData(avant))[names(pData(avant))=="CMSnew"] <- "CMS_Subtypes"
rownames(exprs(avant)) <- fData(avant)$Symbol
CMS.subtypes <- avant$CMS_Subtypes
CMS.subtypes[avant$CMS_Subtypes %in% "unclassified"] <- NA

pData(avant)$MSI[pData(avant)$MSI=="MSI-L"] <- "MSS"
pData(avant)$MSI <- droplevels(pData(avant)$MSI)

### Add expression values of genes of interest to pData
for(i in mylist){
  if(is.na(table(fData(avant)$Symbol==i)[2])==TRUE) {next}
  # Add gene expression to pdata
  pData(avant)[,i] <- exprs(avant)[i,sampleNames(avant)]
  # Cut gene expression into 2 levels at median and add to pdata
  pData(avant)[paste0(i,".Levels")] <- factor(cut(pData(avant)[,i], include.lowest = TRUE,
    breaks=c(-Inf,median(pData(avant)[,i]), +Inf), labels = paste(c("Low", "High"))))
}

### Add values of pathway z-score to pData for AVANT sub-signatures
for(i in names(mygenes)){
  j <- mygenes[i]
  genes <- mygenes[[i]]
  myd <- exprs(avant)[rownames(exprs(avant)) %in% genes,]# extracts expression matrix from the ExpressionSet
  if(is.null(dim(myd))) {next
  } else{
  pData(avant)[[names(j)]] <- colMeans(myd) # adds a column with column means of the signature genes to pData. Note that AVANT data is already median centered and scaled
  # Median cut levels
  pData(avant)[paste0(names(j),".Levels")] <- factor(cut(pData(avant)[,names(j)], include.lowest = TRUE,
    breaks=c(-Inf,median(pData(avant)[,names(j)]), +Inf), labels = paste(c("Low", "High"))))
  # Quartile levels
    pData(avant)[paste0(names(j),".quartiles")] <- cut(
          pData(avant)[,names(j)],
          breaks = quantile(
         pData(avant)[,names(j)],
         probs = seq(0, 1, 0.25)),
        labels = paste( c("lowest", "low", "high", "highest"), "quartile"),
        include.lowest = TRUE
      )
  }
}

### Pathway mean expression values for other signatures of interest
for(i in names(signatures)[-c(7,8)]){
  d <- exprs(avant)[rownames(exprs(avant)) %in% signatures[[i]],]# extracts expression matrix from the ExpressionSet
  pData(avant)[,i] <- colMeans(d) # adds a column with column means of the signature genes to pData
  # Get a standardized signature score
  j <- paste0(i,"_std")
  pData(avant)[,j] <- scale(pData(avant)[,i], center = TRUE, scale = TRUE)
  pData(avant)[paste0(i,".Levels")] <- factor(cut(pData(avant)[,i], include.lowest = TRUE,
    breaks=c(-Inf,median(pData(avant)[,i]), +Inf), labels = paste(c("Low", "High"))))
  pData(avant)[paste0(i,".MeanLevels")] <- factor(cut(pData(avant)[,i], include.lowest = TRUE,
    breaks=c(-Inf,mean(pData(avant)[,i]), +Inf), labels = paste(c("Low", "High"))))
}

### Calculate oncotypeDx score according to weighted approach described in Clark-Langone et al, BMC Cancer 2010. In that reference using PCR-based, housekeeping-normalized gene expression, the authors first calculate an average expression value for the three cell cycle genes (cell cycle group score), an average value for the three stromal genes (stromal group score), and gene GADD45B, reflective of an early response or genotoxic stress pathway, is considered on its own. An unscaled recurrence score RSU is subsequently calculated as 0.1263 x stromal group score - 0.3158 x cell cycle group score + 0.3406 x GADD45B. Lastly, the unscaled recurrence score is rescaled to be between 0 and 100, as follows: RS = 0 if 44.16 x (RSU + 0.30) < 0; RS = 44.16 x (RSU + 0.30) if 0 <= 44.16 x (RSU + 0.30) <=100; and RS = 100 if 44.16 x (RSU + 0.30) >=100. With four out of seven cancer-related oncotypeDx genes available on the Nanostring platform, we adopted the approach as follows: cell cycle group score = average expression of MYC and MYBL2; and stromal group score = BGN expression. Furthermore, we used unscaled recurrence score values as a continuum, for the comparison with our signature. The recurrence score scaling rules were specifically optimized for the dynamic range of PCR, and do not apply to other expression platforms such as Nanostring. 

d <- exprs(avant)[rownames(exprs(avant)) %in% c("MYC", "MYBL2"),]# extracts expression matrix from the ExpressionSet
pData(avant)$OncotypeDx_CellCycleScore <- colMeans(d)
pData(avant)$OncotypeDx_StromalScore <- exprs(avant)[rownames(exprs(avant)) %in% "BGN",]
pData(avant)$OncotypeDx_GADD45B <- exprs(avant)[rownames(exprs(avant)) %in% "GADD45B",]
pData(avant)$OncotypeDx_weighted <- (0.1263*pData(avant)$OncotypeDx_StromalScore) - (0.3158*pData(avant)$OncotypeDx_CellCycleScore) + (0.3406*pData(avant)$OncotypeDx_GADD45B)
pData(avant)$OncotypeDx_weighted_std <- scale(pData(avant)[,"OncotypeDx_weighted"], center = TRUE, scale = TRUE)
pData(avant)[paste0("OncotypeDx_weighted",".Levels")] <- factor(cut(pData(avant)[,"OncotypeDx_weighted"], include.lowest = TRUE,
    breaks=c(-Inf,median(pData(avant)[,"OncotypeDx_weighted"]), +Inf), labels = paste(c("Low", "High"))))
pData(avant)[paste0("OncotypeDx_weighted",".MeanLevels")] <- factor(cut(pData(avant)[,"OncotypeDx_weighted"], include.lowest = TRUE,
    breaks=c(-Inf,mean(pData(avant)[,"OncotypeDx_weighted"]), +Inf), labels = paste(c("Low", "High"))))

### Elastic net scores
mya <- c(0.2)
m <- paste0('alpha_',mya)
mypredict <- data.frame(read.table('../data/AVANT_Elasticnet_GZMBsearch_prediction_scores 0.2.txt')) # Read in prediction scores
pData(avant)$Elasticnet.score <- mypredict$x

# Standardize the elastic net score
  pData(avant)[,"Elasticnet.score_std"] <- scale(
    pData(avant)[,"Elasticnet.score"],
    center = TRUE,
    scale = TRUE
  )
  pData(avant)[,"Elasticnet.score.quartiles"] <- cut(
  	pData(avant)[,"Elasticnet.score_std"],
    breaks = quantile(
    pData(avant)[,"Elasticnet.score_std"],
    probs = seq(0, 1, 0.25)),
    labels = paste( c("lowest", "low", "high", "highest"), "quartile"),
    include.lowest = TRUE
  )

### 3-yr event rate
dfs <- with(
  pData(avant),
  Surv(
	time  = TTMDFS,
    event = CSDFS
  )
)

os <- with(
  pData(avant),
  Surv(
	time  = TTMDIED,
    event = CSDIED
  )
)
