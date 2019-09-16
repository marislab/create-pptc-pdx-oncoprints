
#use genelist
genelist <- lapply(as.matrix(goi.list), function(x)x)
genelist <- unique(genelist)

#subset maf histology
hist.clin <- subset(clin.maf, Histology.Oncoprints == broad.hist)

##subset maf for histology-specific samples only
sub.maf = maftools::subsetMaf(maf = maf, tsb = hist.clin$Tumor_Sample_Barcode, mafObj = TRUE) #all genes for hist-specific
#sub.maf = maftools::subsetMaf(maf = maf, tsb = hist.clin$Model, mafObj = TRUE)
sub.maf.goi = maftools::subsetMaf(maf = sub.maf, genes = genelist, mafObj = TRUE) #only goi

setwd(paste0(subDirHist))
pdf(paste0(subDirHist,"/", each, "-tmp.pdf"), height = 25, width = 25)
try({
  print(maftools::oncoplot(maf = sub.maf.goi, 
         removeNonMutated = F, bgCol = "whitesmoke",
         showTumorSampleBarcodes = F, drawRowBar = F, logColBar = F,
         annotationFontSize = 0, cohortSize = T,
         sortByAnnotation = T, fontSize = 2, legendFontSize = 3,
         colors = col, writeMatrix = T, showTitle = F,
         clinicalFeatures = c("Histology.Detailed", "Phase", "Sex"),
         annotationColor = list(Histology.Detailed = histcol, Sex = sexcol, Phase = phasecol)))
  log("a")}, ##suppress error message with oncoprint plotting so loop continues
  silent=TRUE)
while (!is.null(dev.list()))  dev.off()

unlink(paste0(subDirHist, "/", each, "-tmp.pdf"))
