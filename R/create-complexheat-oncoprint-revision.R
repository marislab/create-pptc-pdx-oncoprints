###call new matrix
#comp.heat.mat <- New.New.onco.matrix
#comp.heat.mat[1:3, 1:3]


##subset clinical file per histology
samples <- subset(clin.maf, Histology.Onco.New == broad.hist)

###subset samples in oncoprint - removing samples without WES:
sample.subset<- subset(samples, Have.maf == "yes")

#subset colors for histologies present
var.in <- names(histcol)[names(histcol) %in% sample.subset$Histology.Detailed.Onco]
histcol2 <-histcol[var.in]

##read in gene list
goi.list <- read.delim(paste0(dataDir, broad.hist, "-goi-list.txt"), sep = "\t",
                       header = F, as.is = T)

##subset maf for histology-specific samples and goi list
  sub.maf.goi = maftools::subsetMaf(maf = maf2, tsb = sample.subset$Tumor_Sample_Barcode, genes = goi.list$V1, mafObj = TRUE) #all genes for hist-specific
  maf.tsb = maftools::subsetMaf(maf = maf2, tsb = sample.subset$Tumor_Sample_Barcode, mafObj = TRUE) #all genes for hist-specific
  ###get top mutated per this subset and goi list
  gene.sum <- mafSummary(sub.maf.goi)$gene.summary  
  ##get top altered genes rather than mutated only genes - first choose genes in goi, then sort by altered, then top N
  #summ.goi <- as.data.frame(subset(gene.sum, Hugo_Symbol %in% goi.list$V1))
  #summ.goi$AlteredSamples <- as.numeric(as.character(gene.sum$AlteredSamples))
  goi.ordered <- gene.sum[order(gene.sum$AlteredSamples, decreasing = T),]
  
  ###select N top genes
  N <- ifelse(nrow(gene.sum) < 50, nrow(gene.sum), 50)
  
  goi.ordered.N <- goi.ordered[1:N,]  
  
  setwd(paste0(subDirHist))
    pdf(paste0(subDirHist,"/",Sys.Date(), "-", broad.hist, "-oncoprint", N, ".pdf"), height = 15, width = 15)
    try({
      print(maftools::oncoplot(maf = maf.tsb, genes = goi.ordered.N$Hugo_Symbol, 
                               #top = N, 
                           removeNonMutated = F, bgCol = "whitesmoke", 
                           showTumorSampleBarcodes = T, drawRowBar = T, sepwd_genes = 1, sepwd_samples = 1,
                           annotationFontSize = 1, gene_mar = 7, barcode_mar = 6,
                           sortByAnnotation = T, fontSize = 1, #legendFontSize = 2,
                           colors = colores, writeMatrix = T, showTitle = F, logColBar = T,
                           clinicalFeatures = c("Histology.Detailed.Onco", "Phase", "Sex", "tp53_score_discrete", "nf1_score_discrete", "Molecular.Subtype"),
                           annotationColor = list(Histology.Detailed.Onco = histcol, 
                                                  Sex = sexcol, Phase = phasecol,
                                                  #Age = colorRamp2(c(0, 5, 10, 15, 20, 45), 
                                                 #                  c("#4F94CD", "#48D1CC", "#FFFACD", "#FF8C00", "#EE2C2C", "#171717")),
                                                  tp53_score_discrete = tp53_score_discrete,
                                                  nf1_score_discrete = nf1_score_discrete,
                                                 Molecular.Subtype = subcol)))
      
      log("a")}, ##suppress error message with oncoprint plotting so loop continues
      silent=TRUE)
    while (!is.null(dev.list()))  dev.off()
    pdf(paste0(subDirHist,"/",Sys.Date(), "-", broad.hist, "-co-occur.pdf"), height = 6, width = 6)
    somaticInteractions(sub.maf.goi)
    dev.off()
  
    
    
    
    
  
