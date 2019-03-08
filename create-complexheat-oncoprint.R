library(ComplexHeatmap)
library(deconstructSigs)

###add arm level dels/amps to annotation
ifelse(broad.hist == "neuroblastoma", 
arm <- read.delim(paste0(dataDir, "arm-lesions.txt"), 
                  as.is = T, header = T, check.names = F), NA)

clin.arm <- merge(clin.sub, arm)


###create mutational burden matrix
sig.df <- pptc.merge[,c("Tumor_Sample_Barcode","Chromosome","Start_position","Reference_Allele","Tumor_Seq_Allele2")]
#require Sample, chr, pos, ref, alt
names(sig.df) <- c("Sample", "chr", "pos", "ref", "alt")
#### Convert to deconstructSigs input
sigs.input <- mut.to.sigs.input(mut.ref = sig.df, 
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt")

mut.sum <-apply(sigs.input,1,sum)
log.mut <- log2(mut.sum)

mat <- New.New.onco.matrix
mat[1:3, 1:3]

###create sample order for oncoprint - by histology detailed and then phase of therapy
samples <- subset(clin.arm, Histology.Oncoprints2 == narrow.hist)
###subset samples in oncoprint - removing those without WES:
sample.sub<- intersect(colnames(mat), samples$Model)
sample.sub.df <- subset(samples, Model%in% sample.sub)
mat2 <- mat[,sample.sub]


hists <- as.list(unique(sample.sub.df$Histology.Oncoprints))
hist <- "neuroblastoma"
#for (hist in hists){
#samples <- subset(clin.sub, Histology.Detailed == hist)
samples <- subset(sample.sub.df, Histology.Oncoprints == hist)
#sample.sub<- intersect(colnames(mat), samples$Model)
mat2 <- mat[,samples$Model %in% colnames(mat)]
#mat2 <- mat[,sample.sub]

#remove empty rows or rows with NA
mat.clean<- mat2[apply(mat2, 1, function(y) !all(is.na(y) | y=="" )),]
tail(mat.clean)
#mat.clean<- mat2[apply(mat2, 1, function(y) !all(is.na(y))),]
#tail(mat.clean)


rownames(samples) <- samples$Model
#sort by histology if > 1 hist, else sort by phase
ifelse(length(unique(samples$Histology.Detailed)) >1, 
       sample_df <- samples[order(samples$Histology.Detailed, samples$Phase),],
       sample_df <- samples[order(samples$Phase),])
sample_order <- rownames(sample_df)

#check for discrepancies
print(setdiff(sample_order, colnames(mat)))
print(setdiff(colnames(mat2),sample_order))

#for annotations, must create new df of rows in the same order of the current matrix prior to sample sorting to match the annotations, else does not print correctly
df_hist <- sample_df[colnames(mat.clean), "Histology.Detailed", drop = FALSE]

#rename to Histology for legend printing
colnames(df_hist)[colnames(df_hist) == "Histology.Detailed"] <- "Histology"
df_phase <- sample_df[colnames(mat.clean), "Phase", drop = FALSE]
df_sex <- sample_df[colnames(mat.clean), "Sex", drop = FALSE]
df_1pdel <- sample_df[colnames(mat.clean), "1p36.33", drop = FALSE]
df_1qgain <- sample_df[colnames(mat.clean), "1q24.3", drop = FALSE]
df_17qgain <- sample_df[colnames(mat.clean), "17q24.1", drop = FALSE]
df_11qdel <- sample_df[colnames(mat.clean), "11q24.3", drop = FALSE]

df_anno <- cbind(df_hist, df_phase, df_sex, df_11qdel, df_17qgain, df_1pdel, df_1qgain)

#create annotation objects
heat_anno = HeatmapAnnotation(df = df_anno,
                              col = list(Histology = histcol,
                                         Phase = phasecol,
                                         Sex = sexcol,
                                         `17q24.1` = `17q24.1`,
                                         `11q24.3` = `11q24.3`,
                                         `1p36.33` = `1p36.33`,
                                         `1q24.3` = `1q24.3`), 
                              annotation_height = 1,
                              na_col = "whitesmoke")
#bar_anno = HeatmapAnnotation(mutations = anno_barplot(log.mut.order), 
#                            gp = c(show_legend = F, annotation_height = 1,
#                                  axis = TRUE, border = FALSE))

pdf(paste(pptc.folder,"Manuscript/figures/oncoprints/", each, "/", Sys.Date(), "-", hist, "-pdx-oncoprint-goi-mut-cn.pdf", sep = ""),
    height = 11, width = 9)
store.plot <- oncoPrint(mat.clean, get_type = function(x) strsplit(x, ";")[[1]],
                        #column_order = sample_order,
                        alter_fun = function(x, y, w, h, v) {
                          n = sum(v)
                          h = h*0.75
                          #background color
                          grid.rect(x, y, w, h, gp = gpar(fill = "whitesmoke", col = NA))
                          # use `names(which(v))` to correctly map between `v` and `col`
                          if(n) grid.rect(x, y - h*0.5 + 1:n/n*h, w*0.75, 1/n*h, 
                                          gp = gpar(fill = col[names(which(v))], col = "whitesmoke"), just = "top")
                        },
                        col = col, 
                        heatmap_legend_param = list(title = "Alteration", at = names(col), 
                                                    labels = mut.labels),
                        show_column_names = T,
                        #top_annotation = NULL,
                        show_row_barplot = F,
                        top_annotation = heat_anno) #or add as bottom_annotation
print(store.plot)

write.table(sample_df[,c("Model", "Histology.Detailed", "Phase", "Sex")], 
            paste0(paste(pptc.folder,"Manuscript/figures/oncoprints/", each, "/", Sys.Date(), 
                         "-", hist, "-oncoprint-sampleorder.txt", sep = "")), sep = "\t", col.names = T, 
            row.names = F, quote = F)

###get order of samples from new matrix if not pre-ordering
df.cols <- as.data.frame(column_order(store.plot)$matrix)
names(df.cols) <- "columns"
mat.order <- mat.clean[,df.cols$columns]
head(mat.clean,2)
head(mat.order,2)

#mutation order for barplot
log.mut.order <-log.mut[colnames(mat.order), drop = FALSE]
log.mut.order.df <- as.data.frame(log.mut.order)
log.mut.order.df$order <- rownames(log.mut.order.df)
log.mut.order.df$order <- factor(log.mut.order.df$order, levels = log.mut.order.df$order)

#make nicer barplots for above
print(ggplot(log.mut.order.df, aes(x = order, y = log.mut.order)) +
        geom_bar(stat="identity", fill = "black") +
        theme_minimal() + scale_y_continuous(limits = (c(0,15)))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1)))

write.table(as.data.frame(log.mut.order), 
            paste0(pptc.folder,"Manuscript/figures/oncoprints/", each, "/", Sys.Date(), 
                   "-", each, "-", hist, "-oncoprint-mutorder.txt"), sep = "\t", col.names = F, 
            row.names = T, quote = F)
dev.off()
#}


