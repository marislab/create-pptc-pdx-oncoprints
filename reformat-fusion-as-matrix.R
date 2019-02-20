#### Read collapsed Fusion File ####
Fusion_File = read.delim("~/Box Sync/PPTC-genomics-collaboration/Manuscript/tables/fusion-results/DriverFusions_Collapsed.txt",
                         sep = "\t", header = T, as.is = T)
#### Get unique list of models from the 3rd column (Models) ####
Model_List = strsplit(as.character(Fusion_File$Models),split = ", ")
Model_List.unique = unique(unlist(Model_List))

#### Change RP11- to RP11_ genes so they separate properly ####
Fusion_File$Fused.Genes <- gsub("RP11-", "RP11_", Fusion_File$Fused.Genes)
Fusion_File$Fused.Genes <- gsub("TRD-GTC9-1", "TRD_GTC9_1", Fusion_File$Fused.Genes)
Fusion_File$Fused.Genes <- gsub("IGHVIII-13-1", "IGHVIII_13_1", Fusion_File$Fused.Genes)

#### Get unique list of individual Genes from the 1st column (Gene fusions) ####
Gene_List = strsplit(as.character(Fusion_File$Fused.Genes),split = "-")
Gene_List.unique = unique(unlist(Gene_List))

#### Replace RP11_ to RP11- genes ####
Gene_List <- lapply(Gene_List, function(X) gsub("RP11_", "RP11-", X))
Gene_List.unique <- gsub("RP11_", "RP11-", Gene_List.unique)
Gene_List <- lapply(Gene_List, function(X) gsub("IGHVIII_13_1", "IGHVIII-13-1", X))
Gene_List.unique <- gsub("IGHVIII_13_1", "IGHVIII-13-1", Gene_List.unique)
Gene_List <- lapply(Gene_List, function(X) gsub("TRD_GTC9_1", "TRD-GTC9-1", X))
Gene_List.unique <- gsub("TRD_GTC9_1", "TRD-GTC9-1", Gene_List.unique)


#### Get line indices for each model ####
Model.Indx = c()
for(i in 1:length(Model_List.unique))
{
  MDL.inx = which(unlist(lapply(Model_List, function(x) any(x==Model_List.unique[i]))))
  Model.Indx = c(Model.Indx,list(MDL.inx))
}

#### Get line indices for each Individual Gene ####
Gene.Indx = c()
for(i in 1:length(Gene_List.unique))
{
  gene.inx = which(unlist(lapply(Gene_List, function(x) any(x==Gene_List.unique[i]))))
  Gene.Indx = c(Gene.Indx,list(gene.inx))
}

#### Create blank fusion matrix ####
Fusion_Matrix = matrix(nrow=length(Gene_List.unique),ncol = length(Model_List.unique),dimnames = list(Gene_List.unique,Model_List.unique))


#### Find Genes with one or mode fusions in the same model ####
#### Loops through all genes and compares the gene indices with the line indices of the models ####
#### A overlap between the indices (%in%) indicates a hit ####
for(i in 1:length(Gene_List.unique))
{
  Fusion_Matrix_Line = lapply(Model.Indx, function(X) 
    {Line.Intersect=intersect(X,Gene.Indx[[i]]);
    if (length(Line.Intersect>0))
    {
      #### Finds intersecting Lines between Genes and Models ####
      #### reorders the genes in the list and finds the unique gene pairs to remove recipricol fusions ####
      L = length(unique(lapply(Gene_List[Line.Intersect],function(X) X[order(X)])))
      if(L==1)return("Fusion");
      if(L>1)return("Multi_Hit_Fusion");
    }
    if (length(Line.Intersect)==0)return("");
    })
  
  #### Add Results for gene iteration to Fusion Matrix ####
  Fusion_Matrix[i,] = unlist(Fusion_Matrix_Line)
  
}#for(i in 1:length(length(Gene_List.unique)))


write.table(Fusion_Matrix, paste0(mainDir, subDir, "/fusion-matrix.txt"), quote = F,col.names = T,row.names = T, sep = "\t")




