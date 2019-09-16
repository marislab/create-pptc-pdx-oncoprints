##convert matrix to dataframe and add rownames as column
cn.df <- as.data.frame(focal.cn.mat)
cn.df$Hugo_Symbol <- rownames(cn.df)
###convert to long dataframe
#cn.df.long <- reshape::melt(cn.df, meas)
cn.df.long <- reshape::melt(cn.df, measure.vars=1:ncol(cn.df)-1, id.vars="Hugo_Symbol")
names(focal.cn.mat)
head(cn.df.long)
table(cn.df.long$value)
##subset, change colnames
colnames(cn.df.long) <- c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification")
##remove NA and blanks
cn.df.na.rm <- cn.df.long %>% tidyr::drop_na(Variant_Classification)
cn.df.complete <- subset(cn.df.na.rm, Variant_Classification != "")
table(cn.df.complete$Variant_Classification)

cntable <- cn.df.complete
cntable$Variant_Classification <- gsub("Amplification", "Amp", cntable$Variant_Classification)
cntable$Variant_Classification <- gsub("Hom_Deletion", "Del", cntable$Variant_Classification)
cntable$Variant_Classification <- gsub("Hem_Deletion", "Del", cntable$Variant_Classification)

