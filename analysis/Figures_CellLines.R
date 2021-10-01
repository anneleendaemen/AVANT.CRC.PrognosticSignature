source("Load_Libraries.R")

####### CODE FOR FIGURES 4H, 4I
####### CODE FOR SUPPLEMENTARY FIGURES S15A-H


######## FIGURE 4H, S15A, S15C, S15E, S15G: Gene expression by indication in a pan-cancer cell line cohort ########

pData.eset <- read.csv("../data/gCell_Colorectal_CMS_Labels.csv", sep=",", header=TRUE, check.names=FALSE)

goi <- c("GZMB", "GZMA", "PRF1", "FCGR3A", "CLEC4C")
pData.eset$SSscore.predicted <- as.character(pData.eset$SSscore.predicted)
pData.eset$SSscore.predicted[is.na(pData.eset$SSscore.predicted)] <- "unclassified"

for(j in goi){
xlabs <- paste(levels(factor(pData.eset$SSscore.predicted)),"\n(N=",table(pData.eset$SSscore.predicted),")",sep="")
p <- ggplot(data=pData.eset, aes_string(x = 'SSscore.predicted', y= as.character(j))) +
  geom_boxplot(aes(fill=factor(SSscore.predicted))) + ylim(0,9) +
  theme_bw() + 
  labs(title=paste0("CRC cell lines: ",j), x = '', y = "Expression") +
  theme(plot.title = element_text(lineheight=1.5, face="bold", hjust=0.5, size=20), 
        axis.title = element_text(lineheight=1.5, face="bold", hjust=0.5,size=16), 
        axis.text = element_text(lineheight=1.5, face="bold", hjust=0.5, colour = "black", size=14),
        legend.position="none") + 
  scale_x_discrete(labels=xlabs) 
m <- as.data.frame(with(pData.eset, pairwise.t.test(get(j), SSscore.predicted))[["p.value"]]) %>% signif(.,2)
m[is.na(m)] <- "-"

if (j=="GZMB") pdf(paste0("../figures/Figure4H_",j,"_AbsoluteLog2Expression_byCMS_CRCcellLines.pdf"), width=6, height=6, onefile=FALSE) else if (j=="GZMA") pdf(paste0("../figures/FigureS15A_",j,"_AbsoluteLog2Expression_byCMS_CRCcellLines.pdf"), width=6, height=6, onefile=FALSE) else if (j=="PRF1") pdf(paste0("../figures/FigureS15C_",j,"_AbsoluteLog2Expression_byCMS_CRCcellLines.pdf"), width=6, height=6, onefile=FALSE) else if (j=="FCGR3A") pdf(paste0("../figures/FigureS15E_",j,"_AbsoluteLog2Expression_byCMS_CRCcellLines.pdf"), width=6, height=6, onefile=FALSE) else if (j=="CLEC4C") pdf(paste0("../figures/FigureS15G_",j,"_AAbsoluteLog2Expression_byCMS_CRCcellLines.pdf"), width=6, height=6, onefile=FALSE)
print(p + annotation_custom(tableGrob(m), xmin=3.5, xmax=4, ymin=7, ymax=8))
dev.off()
}


######## FIGURE 4I, S15B, S15D, S15F, S15H: Gene expression by indication in a pan-cancer cell line cohort ########

### gCell
load("../data/gCell_AllTissues.RData")
annot <- fData(eset)
data <- t(exprs(eset))
colnames(data) <- annot$symbol
pData(eset) <- cbind(pData(eset), data)

goi <- c("GZMB", "GZMA", "PRF1", "FCGR3A", "CLEC4C")
for (i in goi) {
xlabs <- paste(levels(factor(pData(eset)$Tissue)),"\n(N=",table(pData(eset)$Tissue),")",sep="")
p <- ggplot(data=pData(eset), aes_string(x = "Tissue", y= as.character(i))) +
  geom_boxplot(aes(fill=factor(Tissue))) + ylim(0,max(pData(eset)[,goi])) +
  theme_bw() + 
  labs(title=paste0("Cell lines: ",i), x = '', y = "Expression") +
  theme(plot.title = element_text(lineheight=1.5, face="bold", hjust=0.5, size=20), 
        axis.title = element_text(lineheight=1.5, face="bold", hjust=0.5,size=16), 
        axis.text = element_text(lineheight=1.5, face="bold", hjust=0.5, colour = "black", size=14),
        legend.position="none") + 
  scale_x_discrete(labels=xlabs)

if (i=="GZMB") pdf(paste0("../figures/Figure4I_",i,"_AbsoluteLog2Expression_CellLines_LargeCohort.pdf"), width=11, height=6, onefile=FALSE) else if (i=="GZMA") pdf(paste0("../figures/FigureS15B_",i,"_AbsoluteLog2Expression_CellLines_LargeCohort.pdf"), width=11, height=6, onefile=FALSE) else if (i=="PRF1") pdf(paste0("../figures/FigureS15D_",i,"_AbsoluteLog2Expression_CellLines_LargeCohort.pdf"), width=11, height=6, onefile=FALSE) else if (i=="FCGR3A") pdf(paste0("../figures/FigureS15F_",i,"_AbsoluteLog2Expression_CellLines_LargeCohort.pdf"), width=11, height=6, onefile=FALSE) else if (i=="CLEC4C") pdf(paste0("../figures/FigureS15H_",i,"_AbsoluteLog2Expression_CellLines_LargeCohort.pdf"), width=11, height=6, onefile=FALSE)
print(p)
dev.off()
}
