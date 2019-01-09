### specify to not use X, Y chromosomes
chrlengths <- seqlengths(Hsapiens)[paste("chr",c(1:22) ,sep="")]
chr_genome_pos<-list()
for(i in 1:length(chrlengths)) chr_genome_pos[[i]] <- sum(as.numeric(c(0,chrlengths)[1:(i) ]))
chr_genome_pos<-unlist(chr_genome_pos)
names(chr_genome_pos) <-names(chrlengths)

###select only pertinent columns, make gene names rownames
genecn.gistic.out <- gistic.out[,c(4:ncol(gistic.out))]
rownames(genecn.gistic.out)<-gistic.out[,1]
genecn_annot <- genecn.gistic.out[,1:3]

#remove rows containing chr
genecn <- genecn.gistic.out[grep("chr",rownames(genecn.gistic.out),invert=T),]

####harmonize gene IDs - keep approved, fix old
genecn$Previous.symbols <- rownames(genecn)

##subset gene IDs file
gene.ids.subset <- gene.ids[,c("Approved.symbol", "Previous.symbols")]

###split previous symbols into separate rows
expanded <- gene.ids.subset %>% 
  mutate(Previous.symbols = strsplit(as.character(Previous.symbols), ",")) %>% 
  unnest(Previous.symbols)


#matrix of old IDs, update rownames to new IDs
unmatched <- genecn[rownames(genecn) %in% expanded[,"Previous.symbols"],]
old.mer <- merge(unmatched, expanded, all.x = T)
##n rows match - no
identical(nrow(unmatched), nrow(old.mer))
#nrow(unmatched)
#nrow(old.mer)
##which do not match - issue at previous symbols - fix by doing unique
#subset(data.frame(table(old.mer$Approved.symbol)), Freq >1)
#subset(data.frame(table(old.mer$Previous.symbols)), Freq >1)

rownames(old.mer) <- old.mer$Approved.symbol
old.mer$Approved.symbol <- NULL

#do unique before removing symbol column
old.mer <- unique(old.mer)
old.mer$Previous.symbols <- NULL
identical(nrow(unmatched), nrow(old.mer)) #now identical


##subset rows with gene IDs that match current
matched.newIDs <- genecn[genecn[,"Previous.symbols"] %in% gene.ids[,"Approved.symbol"],]
###do all rownames match approved symbols? YES
identical(rownames(matched.newIDs), matched.newIDs$Previous.symbols)
#remove column
matched.newIDs$Previous.symbols <- NULL

##merge matrices back together
updated.cn.mat <- rbind(matched.newIDs, old.mer)



###remove TERT from CN because probes are noisy and falsely calling
genecn.notert <- updated.cn.mat[grep("TERT",rownames(updated.cn.mat),invert=T),]


#broad.hist <- as.list(unique(clin$Histology.Oncoprints))
#oncoprint of genes of interest 
#broad.hist <- "neuroblastoma"

#for (broad.hist in hists){
  ##select correct gene of interest list
  #goi.list <- read.delim(paste0(pptc.folder, "Data/genelists/", broad.hist, "-goi-list.txt"), sep = "\t",
   #               header = F, as.is = T)
  genelist <- lapply(as.matrix(goi.list), function(x)x)
  genelist <- unique(genelist)

#make matrices for broad.hist histology
  sub.clin <- subset(clin, Histology.Oncoprints == broad.hist) 
  genecn2 <- genecn.notert[,which(colnames(genecn.notert) %in% sub.clin$Model)]
  genecn2[genecn2 == 2] <- "Amplification"  
  genecn2[genecn2 == -2] <- "Hom_Deletion"
  genecn2[genecn2 == 0] <- ""
  #genecn2[genecn2 == 1] <- "Gain"
  genecn2[genecn2 == 1] <- ""
  #genecn2[genecn2 == -1] <- "Shallow_Del"

    ###Add het del for TP53 only for all
 for(i in 1:(ncol(genecn2))){ 
    genecn2["TP53",i][genecn2["TP53",i] == -1] = "Hem_Deletion"
 }
  
  ###Add het del for TP53 only for all
  for(i in 1:(ncol(genecn2))){
    genecn2["SMARCB1",i][genecn2["SMARCB1",i] == -1] = "Hem_Deletion"
  }
  
  ###add homozygous deletions for only leukemia for purposes of CDKN2A/B and for TP53 - else, too busy
  ifelse(broad.hist == "leukemia", 
         genecn2[genecn2 == -1] <- "Hem_Deletion", 
         genecn2[genecn2 == -1] <- "")


  mat.goi <- genecn2[rownames(genecn2) %in% genelist,] 
  #remove empty rows
  cn.matrix <- mat.goi[!apply(mat.goi == "", 1, all),]  ##throws an error if only one or zero alterations (carcinoma/rhabdoid)
  #cn.matrix[1:3, 1:3]
  ###add SMARCB1 deletion in 6753ATRT due to manual inspection
  if("ICb-6753ATRT" %in% colnames(cn.matrix)) {
    cn.matrix["SMARCB1","ICb-6753ATRT"] <- "Hom_Deletion"
  }
    ###fix ALL-50 - should be hom del for both CDKN2A/B
  if("ALL-50" %in% colnames(cn.matrix)) {
    cn.matrix["CDKN2B","ALL-50"] <- "Hom_Deletion"
  }
  ###fix PALKTY - should be het del for both CDKN2B
  if("PALKTY" %in% colnames(cn.matrix)) {
    cn.matrix["CDKN2B","PALKTY"] <- "Hem_Deletion"
  }
  
  ##Also, fix CHLA-79 sample to be non-amp - issue with CN seg file
  if("CHLA-79" %in% colnames(cn.matrix)) {
    cn.matrix["MYCN","CHLA-79"] <- ""
  }
  ###NOTE: COG-N-471x is MYCN-amp by pathology, but the amp segment starts after MYCN and was not called
  if("COG-N-471x" %in% colnames(cn.matrix)) {
    cn.matrix["MYCN","COG-N-471x"] <- "Amplification"
  }
  ###add WT1 deletion due to manual inspection
  if("KT-13" %in% colnames(cn.matrix)) {
    cn.matrix["WT1","KT-13"] <- "Hem_Deletion"
    }
  if("NCH-WT-4" %in% colnames(cn.matrix)) {
    cn.matrix["WT1","NCH-WT-4"] <- "Hem_Deletion"
  }
  if("NCH-WT-7" %in% colnames(cn.matrix)) {
    cn.matrix["WT1","NCH-WT-7"] <- "Hem_Deletion"
  }
  if("KT-6" %in% colnames(cn.matrix)) {
    cn.matrix["WT1","KT-6"] <- "Hem_Deletion"
  }  
  if("KT-18" %in% colnames(cn.matrix)) {
    cn.matrix["WT1","KT-18"] <- "Hem_Deletion"
  }  
  if("KT-11" %in% colnames(cn.matrix)) {
    cn.matrix["WT1","KT-11"] <- "Hem_Deletion"
  }
  if("NCH-WT-5" %in% colnames(cn.matrix)) {
    cn.matrix["WT1","NCH-WT-5"] <- "Hem_Deletion"
  }
  if("NCH-WT-6-S13-1506" %in% colnames(cn.matrix)) {
    cn.matrix["WT1","NCH-WT-6-S13-1506"] <- "Hem_Deletion"
  }
  if("NCH-HEP1" %in% colnames(cn.matrix)) {
    cn.matrix["WT1","NCH-HEP1"] <- "Hem_Deletion"
  }
  
  write.table(cn.matrix, paste0(mainDir,subDirHist, "/", broad.hist, "-CN-matrix.txt"),
              sep = "\t", col.names = T, row.names = T, quote = F)
  