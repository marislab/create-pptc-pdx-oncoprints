
prefix_name <- broad.hist
#use genelist
genelist <- lapply(as.matrix(goi.list), function(x)x)
genelist <- unique(genelist)

#subset maf histology
hist.clin <- subset(clin.pptc, Histology.Oncoprints == broad.hist)

##subset maf for histology-specific samples only
sub.maf = maftools::subsetMaf(maf = maf, tsb = hist.clin$Tumor_Sample_Barcode, mafObj = TRUE) #all genes for hist-specific
sub.maf.goi = maftools::subsetMaf(maf = sub.maf, genes = genelist, mafObj = TRUE) #only goi

setwd(paste0(mainDir, subDirHist))
pdf(paste0(mainDir, subDirHist, "/tmp.pdf"))
oncoplot(maf = sub.maf, genes = genelist, 
         removeNonMutated = F, drawRowBar = F,
         showTumorSampleBarcodes = F, titleFontSize = 0, legendFontSize = 8,
         annotationFontSize = 8,  sortByAnnotation = T, fontSize = 14,
         colors = col, writeMatrix = T, #bgCol = "whitesmoke",
         clinicalFeatures = c("Histology.Detailed", "Phase", "Sex"),
         annotationColor = list(Histology.Detailed = histcol, Sex = sexcol, Phase = phasecol))
dev.off()
unlink(paste0(mainDir, subDirHist, "/tmp.pdf"))
