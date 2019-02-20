require(circlize)
###subset samples with DNA assays
clin.sub <- subset(clin, DNA.Part.of.PPTC == "yes")
#call matrix
mat <- New.New.onco.matrix

###create sample order for oncoprint - by histology detailed and then phase of therapy
samples1 <- subset(clin.sub, Histology.Oncoprints == broad.hist)
sub.hists <- as.list(unique(samples1$Histology.Oncoprints2))

for (each in sub.hists){
  #samples <- subset(clin.sub, Histology.Detailed == hist)
  samples <- subset(samples1, Histology.Oncoprints2 == each)
  mat2 <- mat[,colnames(mat) %in% samples$Model]
  
  #remove empty rows or rows with NA
  mat.clean<- mat2[apply(mat2, 1, function(y) !all(is.na(y) | y=="" )),]
  
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
  df_age <- sample_df[colnames(mat.clean), "Age", drop = FALSE]

  df_anno <- cbind(df_hist, df_phase, df_sex, df_age)
  
  #create annotation objects
  heat_anno = HeatmapAnnotation(df = df_anno,
                                col = list(Histology = histcol,
                                           Phase = phasecol,
                                           Sex = sexcol,
                                           Age = colorRamp2(c(0, 5, 10, 15, 20, 45), 
                                                            c("#4F94CD", "#48D1CC", "#FFFACD", "#FF8C00", "#EE2C2C", "#171717"))), 
                                annotation_height = 1, na_col = "whitesmoke")
 
  pdf(paste(mainDir,subDirHist, "/", Sys.Date(), "-", each, "-oncoprint-goi-mut-cn-fusions.pdf", sep = ""),
      height = 17, width = 11)
  store.plot <- oncoPrint(mat.clean, get_type = function(x) strsplit(x, ";")[[1]],
                          ##keep rhabdos separated by histology
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
                          top_annotation = heat_anno) #or add as bottom_annotation,
  
  print(store.plot)
  
  write.table(sample_df[,c("Model", "Histology.Detailed", "Phase", "Sex")], 
              paste0(paste(mainDir,subDirHist, "/", Sys.Date(), 
                           "-", each, "-oncoprint-sampleorder.txt", sep = "")), sep = "\t", col.names = T, 
              row.names = F, quote = F)
  
  ###get order of samples from new matrix if not pre-ordering
  df.cols <- as.data.frame(column_order(store.plot)$matrix)
  names(df.cols) <- "columns"
  mat.order <- mat.clean[,df.cols$columns]
  
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
              paste0(mainDir, subDirHist, "/", Sys.Date(), 
                     "-", each, "-oncoprint-mutorder.txt"), sep = "\t", col.names = F, 
              row.names = T, quote = F)
  dev.off()
  
}


