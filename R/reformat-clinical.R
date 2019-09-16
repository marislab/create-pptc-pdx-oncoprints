###create new clinical file for annotations, including classifier scores and ages
clin.sub <- subset(clin, Have.maf == "yes")
class <- read.delim(paste0(dataDir, "classifier_scores.tsv"), 
                    as.is = T, header = T, check.names = F)

##fix header to enable merge
colnames(class)[colnames(class) == "sample_id"] <- "Model"
class2 <- class[, c("Model", "tp53_score", "nf1_score")]
###recode TP53 and NF1 to discrete
class2$tp53_score_discrete <- apply(class2["tp53_score"], 2, function(x) ifelse(x>=0.55,"Inactive", 
                                                                                ifelse(x<=(0.45),"Active",   
                                                                                       ifelse(x<0.55&x>(0.45),"Modest",x))) ) 
###recode TP53 and NF1 to discrete
class2$nf1_score_discrete <- apply(class2["nf1_score"], 2, function(x) ifelse(x>=0.55,"Inactive", 
                                                                              ifelse(x<=(0.45),"Active",   
                                                                                     ifelse(x<0.55&x>(0.45),"Modest",x))) ) 
###merge new clinical file
clin.class <- merge(clin.sub, class2, all.x = T)

##refactor NA
clin.class$tp53_score_discrete <- ifelse(is.na(clin.class$tp53_score_discrete), "none", clin.class$tp53_score_discrete)
clin.class$nf1_score_discrete <- ifelse(is.na(clin.class$nf1_score_discrete), "none", clin.class$nf1_score_discrete)

##add clinical data
clin.maf = clin.class
clin.hist <- clin.maf[,c("Model", "Histology.Detailed")]

##maftools recognizes TSB as ID, so swap with model for accurate matrix printing
colnames(clin.maf)[colnames(clin.maf) == "Tumor_Sample_Barcode"] <- "TSB"
colnames(clin.maf)[colnames(clin.maf) == "Model"] <- "Tumor_Sample_Barcode"