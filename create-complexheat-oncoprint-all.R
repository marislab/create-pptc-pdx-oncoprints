require(circlize)
clin.sub <- subset(clin, Have.maf == "yes")

###add arm level dels/amps to annotation if NBL or renal
if(broad.hist == "neuroblastoma" | broad.hist == "renal") {
       arm <- read.delim(paste0(dataDir, broad.hist, "-specific-lesions.txt"), 
                         as.is = T, header = T, check.names = F)
       clin.arm <- merge(clin.sub, arm, all.x = T)
       } else { 
       clin.arm <- clin.sub
         }


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

clin.class <- merge(clin.arm, class2, all.x = T)

###call new matrix
comp.heat.mat <- New.New.onco.matrix
comp.heat.mat[1:3, 1:3]

###create sample order for oncoprint - by histology detailed and then phase of therapy
samples1 <- subset(clin.class, Histology.Oncoprints == broad.hist)
sub.hists <- as.list(unique(samples1$Histology.Oncoprints2))

for (each in sub.hists){

###create sample order for oncoprint
  samples <- subset(samples1, Histology.Oncoprints2 == each)
  ###subset samples in oncoprint - removing those without WES:
sample.sub<- intersect(colnames(comp.heat.mat), samples$Model)
sample.sub.df <- subset(samples, Model%in% sample.sub)
mat2 <- comp.heat.mat[,sample.sub]

#remove rows with all NA
mat.clean<- mat2[apply(mat2, 1, function(y) !all(is.na(y) | y=="" )),]
tail(mat.clean)

rownames(samples) <- samples$Model
#sort by histology if > 1 hist, else sort by phase
ifelse(length(unique(samples$Histology.Detailed)) >1, 
       sample_df <- samples[order(samples$Histology.Detailed, samples$Phase),],
       sample_df <- samples[order(samples$Phase),])
sample_order <- rownames(sample_df)
  
 
#for annotations, must create new df of rows in the same order of the current matrix prior to sample sorting to match the annotations, else does not print correctly
 df_hist <- sample_df[colnames(mat.clean), "Histology.Detailed", drop = FALSE]
  
#rename to Histology for legend printing
 colnames(df_hist)[colnames(df_hist) == "Histology.Detailed"] <- "Histology"  

 df_phase <- sample_df[colnames(mat.clean), "Phase", drop = FALSE]
 df_sex <- sample_df[colnames(mat.clean), "Sex", drop = FALSE]
 df_tp53 <- sample_df[colnames(mat.clean), "tp53_score_discrete", drop = FALSE]
 df_nf1 <- sample_df[colnames(mat.clean), "nf1_score_discrete", drop = FALSE]
 df_age <- sample_df[colnames(mat.clean), "Age", drop = FALSE]
 
if(broad.hist != "neuroblastoma" & broad.hist != "renal") {
  #combine all
  df_anno <- cbind(df_hist, df_phase, df_sex, df_age, df_nf1, df_tp53)
  #create annotation objects
  heat_anno = HeatmapAnnotation(df = df_anno,
                                col = list(Histology = histcol,
                                           Phase = phasecol,
                                           Sex = sexcol,
                                           Age = colorRamp2(c(0, 5, 10, 15, 20, 45), 
                                                            c("#4F94CD", "#48D1CC", "#FFFACD", "#FF8C00", "#EE2C2C", "#171717")),
                                           tp53_score_discrete = tp53_score_discrete,
                                           nf1_score_discrete = nf1_score_discrete),
                                annotation_height = 1,
                                na_col = "whitesmoke")
  bar_anno = HeatmapAnnotation(mutations = anno_barplot(log.mut.order), 
                               gp = c(show_legend = F, annotation_height = 1,
                                      axis = TRUE, border = FALSE))
}
 
if(broad.hist == "neuroblastoma") {
  df_1pdel <- sample_df[colnames(mat.clean), "1p36.33", drop = FALSE]
  df_17qgain <- sample_df[colnames(mat.clean), "17q24.1", drop = FALSE]
  df_11qdel <- sample_df[colnames(mat.clean), "11q24.3", drop = FALSE]
  df_anno <- cbind(df_hist, df_phase, df_sex, df_age, df_nf1, df_tp53, df_11qdel, df_17qgain, df_1pdel)
  #create annotation objects
  heat_anno = HeatmapAnnotation(df = df_anno,
                                col = list(Histology = histcol,
                                           Phase = phasecol,
                                           Sex = sexcol,
                                           Age = colorRamp2(c(0, 5, 10, 15, 20, 45), 
                                                            c("#4F94CD", "#48D1CC", "#FFFACD", "#FF8C00", "#EE2C2C", "#171717")), 
                                           tp53_score_discrete = tp53_score_discrete,
                                           nf1_score_discrete = nf1_score_discrete,
                                           `17q24.1` = `17q24.1`,
                                           `11q24.3` = `11q24.3`,
                                           `1p36.33` = `1p36.33`), 
                                annotation_height = 1,
                                na_col = "whitesmoke")
 } else { df_anno <- df_anno
  }

if(broad.hist == "renal"){
     df_1p <- sample_df[colnames(mat.clean), "1p", drop = FALSE]
     df_1q <- sample_df[colnames(mat.clean), "1q", drop = FALSE]
     df_11p13 <- sample_df[colnames(mat.clean), "11p13", drop = FALSE]
     df_11p15 <- sample_df[colnames(mat.clean), "11p15.5", drop = FALSE]
     df_16q <- sample_df[colnames(mat.clean), "16q", drop = FALSE]
     df_Xq11 <- sample_df[colnames(mat.clean), "Xq11.2", drop = FALSE]
     df_anno <- cbind(df_hist, df_phase, df_sex, df_age,df_nf1, df_tp53,df_1p, df_1q, df_11p13, df_11p15, df_16q, df_Xq11)
     #create annotation objects
     heat_anno = HeatmapAnnotation(df = df_anno,
                                   col = list(Histology = histcol,
                                              Phase = phasecol,
                                              Sex = sexcol,
                                              Age = colorRamp2(c(0, 5, 10, 15, 20, 45), 
                                                               c("#4F94CD", "#48D1CC", "#FFFACD", "#FF8C00", "#EE2C2C", "#171717")), 
                                              tp53_score_discrete = tp53_score_discrete,
                                              nf1_score_discrete = nf1_score_discrete,
                                              `1p` = `1p`,
                                              `1q` = `1q`,
                                              `11p13` = `11p13`,
                                              `11p15.5` = `11p15.5`,
                                              `16q` = `16q`,
                                              `Xq11.2` = `Xq11.2`), 
                                   annotation_height = 1,
                                   na_col = "whitesmoke")
} else { df_anno <- df_anno
}




 pdf(paste(subDirHist, "/", Sys.Date(), "-", each, "-oncoprint-goi-mut-cn-fusions.pdf", sep = ""),
     height = 20, width = 14)
store.plot <-  oncoPrint(mat.clean, get_type = function(x) strsplit(x, ";")[[1]],
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
                  show_row_barplot = T,
                  top_annotation = heat_anno) #or add as bottom_annotation
  
print(store.plot)
  
  
write.table(sample_df[,c("Model", "Histology.Detailed", "Phase", "Sex")], 
            paste0(paste(subDirHist, "/", Sys.Date(), 
                         "-", each, "-oncoprint-sampleorder.txt", sep = "")), sep = "\t", col.names = T, 
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
              paste0(subDirHist, "/", Sys.Date(), 
                     "-", each, "-oncoprint-mutorder.txt"), sep = "\t", col.names = F, 
              row.names = T, quote = F)
  dev.off()
  
}

