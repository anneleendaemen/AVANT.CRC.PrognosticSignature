### Immune genes of interest
mylist <- c('CD28','CXCL1','TNF', 'FCRL5','TAP1','CXCR4','GZMB','TP73','PDCD1', # AVANT immune signature genes
'CD8A', 'GZMA','IFNG','CXCL9','CXCL10','PRF1','TBX21','KLRK1', # T-effector genes
   'TCF4', 'NRP1') # pDC genes

### AVANT signature
myd <- read.csv('../data/Predictive_signature_genes_pathways_alpha0.2.csv')
Proliferative_signature <- myd$gene[myd$Signatures=='Proliferative_signature']
Stromal_signature_TGFb_enriched <- myd$gene[myd$Signatures=='Stromal_signature_TGFb_enriched']
Stromal_signature_ECM_endocrine <- myd$gene[myd$Signatures=='Stromal_signature_ECM_endocrine']
Stromal_signature <- myd$gene[myd$Signatures=='Stromal_signature_ECM_endocrine'|myd$Signatures=='Stromal_signature_TGFb_enriched']
Immune_signature <- myd$gene[myd$Signatures=='Immune_signature']

mygenes <- list(Proliferative_signature=Proliferative_signature,
  Stromal_signature_TGFb_enriched=Stromal_signature_TGFb_enriched,
  Stromal_signature_ECM_endocrine = Stromal_signature_ECM_endocrine,
  Stromal_signature=Stromal_signature,
  Immune_signature = Immune_signature)

### Other signatures of interest
CAF_Isella_et_al <- c(
  "INHBA","PTN","LOC100132116","MEG3","CCL2","LOXL2","GPC6","BNC2",
  "PDPN","RAB31","GUCY1B3","EFEMP1","TGFB3","ADAM12","C1R","SGCE",
  "GLT8D2","COL12A1","HGF","DPT","CPXM1","SFRP2","THY1","C1S","FBN1",
  "CDR1","FHL1","LMOD1","NCAM1","DDR2","ANTXR1","CYGB","TRO","CCDC80",
  "HTRA3","CCL11","EDIL3","GALNTL2","LHFP","GGT5","ACAN","UNC5C",
  "WISP1","ABCA9","PDE1A","ADAMTSL3","EVC","MOXD1","SGCD","MAB21L2",
  "C9orf47","CCL13","MFAP5","SCARA5","PLXDC1","PLAT","FST","SFRP1",
  "PDGFRA","BOC","DZIP1","PLEKHH2","COL14A1","DIO2","SLC26A10",
  "COL8A1","ISLR","COL11A1","TNFAIP6","F2RL2","PCDH18","OLFML2B",
  "KCNJ8","KIAA1755","4-Sep","FOXF1","MRGPRF","ASPN","CH25H",
  "FOXL1","FNDC1","WNT2B","PTGIS","WBSCR17","FBN2","ITGA11",
  "SFRP4","CPZ","LOX","MGC24103","WNT2","FIBIN","RASSF8","SMO",
  "MMP3","MEIS1","WNT5A","CHRDL1","RGMA","TNFSF11","EPHA3",
  "KCNE4","PAMR1","HHIP","GLI2","PPAPDC1A","SMARCD3","GRP",
  "LRRC17","LRRN4CL","COLEC12","TFPI2","SPOCK1","TSHZ3","PTGER3",
  "ADAMDEC1","MMP19","PCDH7","ABCA6","TMEM119","MASP1","GAS1"
  ,"CYP7B1","CCDC102B","NEXN","FAM65C","OLFML1","DSEL","CD302",
  "ADH1B","OGN"
)
OncotypeDx <- c("MKI67","MYC","MYBL2","FAP","BGN","INHBA", "GADD45B")

Teffector_signature <- c('CD8A', 'GZMA','GZMB','IFNG','CXCL9','CXCL10','PRF1','TBX21') 

Teffector_signature_modified <- c('CD8A', 'IFNG','CXCL9','CXCL10','PRF1','TBX21') # T-effector signature excluding GZMB and GZMA

Teffector_signature_noGZMB <- c('CD8A', 'GZMA','IFNG','CXCL9','CXCL10','PRF1','TBX21') # T-effector signature excluding GZMB

F_TBRs_Calon <- c('AADACP1','ANGPTL4','ANKRD44','ANOS1','APBB2','ARHGEF3','ARHGEF40','AUTS2','BHLHE40','BPGM','C3orf52','C4orf26','CACHD1','CALB2','CALD1','CDH6','CDKN2B','CHST11','CILP','CLDN14','CLDN4','CNNM4','CNTN1','COL10A1','COL27A1','CTGF','DAAM1','DACT1','DDX10','DHRS2','DLX2','DNAJB5','DNAJC18','DNM3OS','DOCK10','E2F7','EDN1','EFNB2','EGR2','ELMOD1','EPHA4','ERN1','ESM1','ETV6','F2RL1','FBXO32','FGF1','FGF18','FLT1','FN1','FNIP2','FOXP1','FRMD4A','FUT4','GADD45B','GAS1','GOPC','GPR161','GPR183','GRB14','GZMK','HAS2','HBEGF','HEY1','HIC1','HIVEP2','HS3ST3B1','IFIH1','IGFBP3','IL11','IL6','INHBA','ITGB6','JARID2','JUNB','KANK4','KDM6B','KDM7A','KIAA1755','KIF26B','KLF7','LIF','LINC00312','LMCD1','LMO4','LOC101927482','LOC101928955','LOC339260','LRRC8C','MBOAT2','MEDAG','MEGF9','MEOX1','MEX3B','MIR143HG','MIR181A2HG','MIR503HG','MSC','MSX2','MTSS1','MURC','NEDD9','NET1','NGF','NKX3-1','NOX4','NREP','NUAK1','OSGIN2','PALLD','PCDH9','PDGFC','PDLIM4','PDPN','PELI1','PGBD5','PGM2L1','PIK3CD','PKNOX2','PLAUR','PLCE1','PLEK2','PMEPA1','PODXL','PPP1R14C','PRDM1','PRR5L','PRR9','PTGS2','PTHLH','RAB30','RASD1','RASGRP3','RASL12','RHOU','RNF150','RYBP','S1PR5','SDC1','SEMA7A','SERPINE1','SETBP1','SKIL','SLC19A2','SLC35F2','SLC35F3','SLC46A3','SMAD7','SNAI1','SNCAIP','SNORD114-3','SNX30','SORBS2','SOX4','SOX6','SPSB1','STEAP2','STK17B','STK38L','SYNE1','TAGLN3','TBX3','TGFB2','TIMP3','TMEM2','TNC','TNFAIP6','TRIB1','TSHZ3','TSPAN2','TUFT1','VEGFA','VEPH1','VMP1','WNT2','WNT9A','YIPF5','ZEB1','ZNF365')

NK_signature <- c("KLRC3", "KLRK1", "NKG7", "KLRC2", "KLRD1")

pDC_CIWG <- c("BDCA2", "BDCA4", "CD123", "CLEC4C", "HLA-DR", "IL3RA", "NRP1", "TCF4")

signatures <- list(OncotypeDx = OncotypeDx, 
                   Teffector_signature = Teffector_signature,
                   Teffector_signature_modified = Teffector_signature_modified,
  				   Teffector_signature_noGZMB = Teffector_signature_noGZMB,
                   CAF_Isella_et_al=CAF_Isella_et_al,
                   F_TBRs_Calon = F_TBRs_Calon,
                   NK_signature = NK_signature,
                   pDC_CIWG = pDC_CIWG
                   )
