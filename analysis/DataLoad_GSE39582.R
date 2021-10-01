source("Load_Libraries.R")
source("Load_Signatures.R")


####### DATA PREPARATION FOR GSE39582, MARISA ET AL ########

### Load preprocessed microarray data
eset <- get(load("../data/GSE39582.rdata"))
pData(eset)$cit.molecularsubtype <- factor(pData(eset)$cit.molecularsubtype)
#Z-score transformation
exprs(eset) <- t(scale(t(exprs(eset)), center=TRUE, scale=TRUE))

### CMS subtypes in GSE39582
SScms <- classifyCMS(exprs(eset),method="SSP")[[3]]# Genes are rows, samples are columns
pData(eset)$CMS <- SScms$SSP.nearestCMS

### Add expression values of genes of interest to pData
for(i in mylist){
  if(is.na(table(fData(eset)$SYMBOL==i)[2])==TRUE) {next}
  # Get the gene ID for gene symbol
  j <- rownames(fData(eset))[fData(eset)$SYMBOL==i]
  # Add gene expression to pdata
  pData(eset)[,i] <- exprs(eset)[j,sampleNames(eset)]
  # Cut gene expression into 2 levels at median and add to pdata
  pData(eset)[paste0(i,".Levels")] <- factor(cut(pData(eset)[,i], 
    breaks = c(-Inf, 
    median(pData(eset)[,i]) ,
    +Inf),
  labels =  c("low", "high"),
  include.lowest = TRUE))
}

### Pathway mean expression values for signatures of interest
for(i in 1:length(signatures)){
  j <- signatures[i]
  genes <- signatures[[i]]
  myd <- exprs(eset)[fData(eset)$SYMBOL %in% genes,]# extracts expression matrix from the ExpressionSet
  dim(myd)
  # Z-score
  myd <- t(scale(t(myd), center=TRUE, scale=TRUE))
  pData(eset)[[names(j)]] <- colMeans(myd) # adds a column with column means of the signature genes to pData
  # Median cut signature
  pData(eset)[[paste0(names(j),".levels")]] <- cut(
  pData(eset)[,names(j)],
  breaks = c(-Inf, 
    median(pData(eset)[,names(j)]) ,
    +Inf),
  labels =  c("low", "high"),
  include.lowest = TRUE)
  # Mean cut signature
  pData(eset)[[paste0(names(j),".MeanLevels")]] <- cut(
  pData(eset)[,names(j)],
  breaks = c(-Inf, 
    mean(pData(eset)[,names(j)]) ,
    +Inf),
  labels =  c("low", "high"),
  include.lowest = TRUE)
  # Quartile cut signature
  pData(eset)[[paste0(names(j),".quartiles")]] <- cut(
  pData(eset)[,names(j)],
  breaks = quantile(
    pData(eset)[,names(j)],
    probs = seq(0, 1, 0.25)
  ),
  labels = paste( c("lowest", "low", "high", "highest"), "quartile"),
  include.lowest = TRUE
)
  i <- paste0(names(j),"_std")
  pData(eset)[,i] <- scale(
    pData(eset)[,names(j)],
    center = TRUE,
    scale = TRUE
  )

}

### Add values of pathway z-score to pData for AVANT sub-signatures
for(i in names(mygenes)){
  j <- mygenes[i]
  genes <- mygenes[[i]]
  myd <- exprs(eset)[fData(eset)$SYMBOL %in% genes,]# extracts expression matrix from the ExpressionSet
  dim(myd)
  #myd[1:5,1:5]
  pData(eset)[[names(j)]] <- colMeans(myd) # adds a column with column means of the signature genes to pData
  pData(eset)[[paste0(names(j),".levels")]] <- cut(
  pData(eset)[,names(j)],
  breaks = c(-Inf, 
    median(pData(eset)[,names(j)]) ,
    +Inf),
  labels =  c("Low", "High"),
  include.lowest = TRUE)
}

### Recurrence free survival
rfs <- with(
  pData(eset),
  Surv(
    rfs.delay,
    rfs.event
  )
)

### Read in gene signature for selected alpha = 0.2
m <- paste0('alpha_',0.2)
predictor.genes <- data.frame(read.table('../data/AVANT_Elasticnet_GZMBsearch_prognosticgenes 0.2.txt'))
predictor.genes <- predictor.genes[predictor.genes$gene %in% fData(eset)$SYMBOL,]
rownames(predictor.genes) <- predictor.genes$gene
predictor.genes <- predictor.genes[c(2)]

# Subset eset to match the predictor genes
test <- na.omit(exprs(eset)[match(rownames(predictor.genes),fData(eset)$SYMBOL),])
rownames(test) <- fData(eset)$SYMBOL[match(rownames(test),rownames(fData(eset)))]

# weighted_directional_score
glmnet <- t(as.matrix(predictor.genes)) %*% test
glmnet <- glmnet[1,]
pData(eset)$glmnet.weighted <- glmnet

pData(eset)$glmnetW.quartiles <- cut(
  pData(eset)$glmnet.weighted,
  breaks = quantile(
    pData(eset)$glmnet.weighted,
    probs = seq(0, 1, 0.25)
  ),
  labels = paste( c("lowest", "low", "high", "highest"), "quartile"),
  include.lowest = TRUE
)

# unweighted_directional_score
predictor.sign <- sign(predictor.genes)
glmnet <- t(as.matrix(predictor.sign)) %*% test  ## use only the sign of the coefficients
glmnet <- glmnet[1,]/length(predictor.sign)
pData(eset)$glmnet.unweighted <- glmnet

pData(eset)$glmnetU.quartiles <- cut(
  pData(eset)$glmnet.unweighted,
  breaks = quantile(
    pData(eset)$glmnet.unweighted,
    probs = seq(0, 1, 0.25)
  ),
  labels = paste( c("lowest", "low", "high", "highest"), "quartile"),
  include.lowest = TRUE
)

# Standardize the elastic net scores
pData(eset)[,"glmnet.unweighted_std"] <- scale(
    pData(eset)[,"glmnet.unweighted"],
    center = TRUE,
    scale = TRUE
  )
pData(eset)[,"glmnet.weighted_std"] <- scale(
    pData(eset)[,"glmnet.weighted"],
    center = TRUE,
    scale = TRUE
  )

### Calculate oncotypeDx score according to weighted approach described in Clark-Langone et al, BMC Cancer 2010. In that reference using PCR-based, housekeeping-normalized gene expression, the authors first calculate an average expression value for the three cell cycle genes (cell cycle group score), an average value for the three stromal genes (stromal group score), and gene GADD45B, reflective of an early response or genotoxic stress pathway, is considered on its own. An unscaled recurrence score RSU is subsequently calculated as 0.1263 x stromal group score - 0.3158 x cell cycle group score + 0.3406 x GADD45B. Lastly, the unscaled recurrence score is rescaled to be between 0 and 100, as follows: RS = 0 if 44.16 x (RSU + 0.30) < 0; RS = 44.16 x (RSU + 0.30) if 0 <= 44.16 x (RSU + 0.30) <=100; and RS = 100 if 44.16 x (RSU + 0.30) >=100. We used unscaled recurrence score values as a continuum, for the comparison with our signature. The recurrence score scaling rules were specifically optimized for the dynamic range of PCR, and do not apply to other expression platforms such as Nanostring. 

d <- exprs(eset)[fData(eset)$SYMBOL %in% c("MYC", "MYBL2", "MKI67"),]
pData(eset)$OncotypeDx_CellCycleScore <- colMeans(d)
d <- exprs(eset)[fData(eset)$SYMBOL %in% c("BGN", "FAP", "INHBA"),]
pData(eset)$OncotypeDx_StromalScore <- colMeans(d)
pData(eset)$OncotypeDx_GADD45B <- exprs(eset)[fData(eset)$SYMBOL %in% "GADD45B",]
pData(eset)$OncotypeDx_weighted <- (0.1263*pData(eset)$OncotypeDx_StromalScore) - (0.3158*pData(eset)$OncotypeDx_CellCycleScore) + (0.3406*pData(eset)$OncotypeDx_GADD45B)
pData(eset)$OncotypeDx_weighted_std <- scale(pData(eset)[,"OncotypeDx_weighted"], center = TRUE, scale = TRUE)
pData(eset)[paste0("OncotypeDx_weighted",".MeanLevels")] <- factor(cut(pData(eset)[,"OncotypeDx_weighted"], include.lowest = TRUE,
    breaks=c(-Inf,mean(pData(eset)[,"OncotypeDx_weighted"]), +Inf), labels = paste(c("Low", "High"))))
