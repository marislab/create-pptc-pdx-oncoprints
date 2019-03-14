##add clinical data
clin.maf = clin
clin.hist <- clin.maf[,c("Model", "Histology.Detailed")]

##maftools recognizes TSB as ID, so swap with model for accurate matrix printing
colnames(clin.maf)[colnames(clin.maf) == "Tumor_Sample_Barcode"] <- "TSB"
colnames(clin.maf)[colnames(clin.maf) == "Model"] <- "Tumor_Sample_Barcode"

#clin.pptc <- subset(clin.maf, Model.Part.of.PPTC == "yes")


maf = read.maf(maf = pptc.merge, clinicalData = clin.maf, vc_nonSyn = c("Frame_Shift_Del", 
                                                                         "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", 
                                                                         "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation",  
                                                                         "Stop_Codon_Ins", "Start_Codon_Del"))

load <- getSampleSummary(maf)
load$subs <- load$Missense_Mutation + load$Nonsense_Mutation
#used Roche Nimblegen VCRome v. 2.1 = 45.1 Mb capture
#for muts per mb, only use substitutions 
load$MutperMB <- round(load$subs/45.1,2)
colnames(load)[colnames(load)== "Tumor_Sample_Barcode"] <- "Model"

#merge with histology detailed
load.hist <- merge(clin.hist, load)

##write in sort order most to least Mut per MB
write.table(load.hist[order(-load.hist$MutperMB),], paste0(subDir, Sys.Date(), "-mutations-per-model.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)
