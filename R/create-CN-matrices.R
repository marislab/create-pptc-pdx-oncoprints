  genelist <- lapply(as.matrix(goi.list), function(x)x)
  genelist <- unique(genelist)

#make matrices for broad.hist histology
  sub.clin <- subset(clin, Histology.Oncoprints == broad.hist) 
  genecn2 <- focal.cn.mat[,which(colnames(focal.cn.mat) %in% sub.clin$Model)]
  ##subset for goi
  cn.matrix <- genecn2[rownames(genecn2) %in% genelist,] 
  ##subset for only models of interest
  final.cn.mat <- cn.matrix[,colnames(cn.matrix) %in% sub.clin$Model]
  
  ##Write new CN matrix
  write.table(final.cn.mat, paste0(subDirHist, "/", broad.hist, "-CN-matrix.txt"),
              sep = "\t", col.names = T, row.names = T, quote = F)
  
  
