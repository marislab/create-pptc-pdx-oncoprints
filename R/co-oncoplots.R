###co-oncoplots for dx/relapse

##subset clinical file per histology
###leukemias
if(broad.hist == "leukemia"){
  for (each in c("T-ALL", "BCP-ALL")){
  samples <- subset(clin.maf, Histology.Detailed == each)
  ###subset samples in oncoprint - removing samples without WES:
  sample.subset<- subset(samples, Have.maf == "yes")

  #subset colors for histologies present
  var.in <- names(histcol)[names(histcol) %in% sample.subset$Histology.Detailed]
  histcol2 <-histcol[var.in]

  ##read in gene list
  all.genes <- read.delim(paste0(dataDir, "all-goi-list.txt"), sep = "\t",
                          header = F, as.is = T)
  goi.list <- read.delim(paste0(dataDir, broad.hist, "-goi-list-old.txt"), sep = "\t",
                         header = F, as.is = T)
  
  dx <- subset(samples, Phase2 == "Diagnosis")
  rel <- subset(samples, Phase2 == "Relapse")
  
  ##subset maf for histology-specific samples and all goi list - need to use all or else subset does not work for samples with no gene mutations
  dx.maf = maftools::subsetMaf(maf = maf2, tsb = dx$Tumor_Sample_Barcode, genes = all.genes$V1, mafObj = TRUE) 
  rel.maf = maftools::subsetMaf(maf = maf2, tsb = rel$Tumor_Sample_Barcode, genes = all.genes$V1, mafObj = TRUE) 
  ###get top mutated in relapse per this subset and goi list
  rel.gene.sum <- mafSummary(rel.maf)$gene.summary  
  ###subset for genes in the histology-specific list
  gene.sum <- subset(rel.gene.sum, Hugo_Symbol %in% goi.list$V1) 
  ##get top altered genes rather than mutated only genes - first choose genes in goi, then sort by altered, then top N
  #summ.goi <- as.data.frame(subset(gene.sum, Hugo_Symbol %in% goi.list$V1))
  #summ.goi$AlteredSamples <- as.numeric(as.character(gene.sum$AlteredSamples))
  goi.ordered <- gene.sum[order(gene.sum$AlteredSamples, decreasing = T),]

  ###select N top genes
  N <- ifelse(nrow(goi.ordered) > 20, 20, nrow(goi.ordered))
  #N <- 10
  goi.ordered.N <- goi.ordered[1:N,]  

  setwd(paste0(subDirHist))
  pdf(paste0(subDirHist,"/",Sys.Date(), "-", unique(samples$Histology.Detailed), "-co-oncoprint.pdf"), height = 5, width = 10)
  try({
   print(maftools::coOncoplot(m1=dx.maf, m2=rel.maf, m1Name = "Diagnosis", m2Name = "Relapse",
                              genes = goi.ordered.N$Hugo_Symbol, removeNonMutated = F, 
                              colors = colores, bgCol = "whitesmoke",
                              sepwd_genes1 = 1, sepwd_genes2= 1, 
                              sepwd_samples1 = 1,sepwd_samples2 = 1,
                              legendFontSize = 0))
                        

  log("a")}, ##suppress error message with oncoprint plotting so loop continues
  silent=TRUE)
  while (!is.null(dev.list()))  dev.off()

  }
}

if(broad.hist == "solid"){
  for (each in c("Neuroblastoma", "Osteosarcoma")){
    samples <- subset(clin.maf, Histology.Detailed == each)
    ###subset samples in oncoprint - removing samples without WES:
    sample.subset<- subset(samples, Have.maf == "yes")
    
    #subset colors for histologies present
    var.in <- names(histcol)[names(histcol) %in% sample.subset$Histology.Detailed]
    histcol2 <-histcol[var.in]
    
    ##read in gene list
    all.genes <- read.delim(paste0(dataDir, "all-goi-list.txt"), sep = "\t",
                           header = F, as.is = T)
    goi.list <- read.delim(paste0(dataDir, each, "-goi-list.txt"), sep = "\t",
                           header = F, as.is = T)
    
    dx <- subset(samples, Phase2 == "Diagnosis")
    rel <- subset(samples, Phase2 == "Relapse")
    
    ##subset maf for histology-specific samples and all goi list - need to use all or else subset does not work for samples with no gene mutations
    dx.maf = maftools::subsetMaf(maf = maf2, tsb = dx$Tumor_Sample_Barcode, genes = all.genes$V1, mafObj = TRUE) 
    rel.maf = maftools::subsetMaf(maf = maf2, tsb = rel$Tumor_Sample_Barcode, genes = all.genes$V1, mafObj = TRUE) 
    ###get top mutated in relapse per this subset and goi list
    rel.gene.sum <- mafSummary(rel.maf)$gene.summary  
    ###subset for genes in the histology-specific list
    gene.sum <- subset(rel.gene.sum, Hugo_Symbol %in% goi.list$V1)
    ##get top altered genes rather than mutated only genes - first choose genes in goi, then sort by altered, then top N
    #summ.goi <- as.data.frame(subset(gene.sum, Hugo_Symbol %in% goi.list$V1))
    #summ.goi$AlteredSamples <- as.numeric(as.character(gene.sum$AlteredSamples))
    goi.ordered <- gene.sum[order(gene.sum$AlteredSamples, decreasing = T),]
    
    ###select N top genes
    N <- ifelse(nrow(goi.ordered) > 20, 20, nrow(goi.ordered))
    #N <- 10
    goi.ordered.N <- goi.ordered[1:N,]  
    
    setwd(paste0(subDirHist))
    pdf(paste0(subDirHist,"/",Sys.Date(), "-", unique(samples$Histology.Detailed), "-co-oncoprint.pdf"), height = 3, width = 16)
    try({
      print(maftools::coOncoplot(m1=dx.maf, m2=rel.maf, m1Name = "Diagnosis", m2Name = "Relapse",
                                 genes = goi.ordered.N$Hugo_Symbol, removeNonMutated = F, 
                                 colors = colores, bgCol = "whitesmoke",
                                 sepwd_genes1 = 1, sepwd_genes2= 1, 
                                 sepwd_samples1 = 1,sepwd_samples2 = 1,
                                 legendFontSize = 0))
      
      
      log("a")}, ##suppress error message with oncoprint plotting so loop continues
      silent=TRUE)
    while (!is.null(dev.list()))  dev.off()
    
  }
}





