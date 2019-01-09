#### Read fusion file ####
#Fusion_Matrix <- as.matrix(read.delim(paste0(pptc.folder, "Manuscript/figures/oncoprints/Fusions_Matrix.txt"),
#                           sep = "\t", header = T, as.is = T, check.names = F))

colnames(New.onco.matrix) %in% colnames(Fusion_Matrix)
dim(New.onco.matrix)
dim(Fusion_Matrix)
#### Rearranges CN file into same order as the onco file ####
#### Removes Genes that aren't present in Onco file ####
Fusion_Matrix.rearranged=Fusion_Matrix[match(rownames(New.onco.matrix),rownames(Fusion_Matrix)),]

#### Determine Genes in cn but not in onco ####
Unmatched_Genes = setdiff(rownames(Fusion_Matrix),rownames(New.onco.matrix))

#### Create matrix with unmatched brain cn genes to be appended to rearranged matrix ####
Unmatched_cn_Gene_matrix = Fusion_Matrix[rownames(Fusion_Matrix) %in% Unmatched_Genes,]

#### Create blank matrix of unmatched genes to be appeneded to onco matrix ####
Onco_Unmatched_Gene_append.matrix = matrix(nrow=nrow(Unmatched_cn_Gene_matrix),ncol=ncol(New.onco.matrix))
rownames(Onco_Unmatched_Gene_append.matrix) = rownames(Unmatched_cn_Gene_matrix)
colnames(Onco_Unmatched_Gene_append.matrix) = colnames(New.onco.matrix)

#### Append unmatched genes matrix to Brain cn matrix ####
Fusion_Matrix.rearranged = rbind(Fusion_Matrix.rearranged,Unmatched_cn_Gene_matrix)

#### Make new matrix coverting Factors to Characters and add rows for unmatched Brain cn genes ####
New.New.onco.matrix = rbind(as.matrix(New.onco.matrix),Onco_Unmatched_Gene_append.matrix)


#### Replace NAs with "" ####
New.New.onco.matrix[is.na(New.New.onco.matrix)] = ""

#### Check Each onco column name for a corresponding column name in the Brain CN file ####
for(i in 1:ncol(New.onco.matrix))
{
  
  COL.inx = which(colnames(Fusion_Matrix.rearranged) == colnames(New.onco.matrix)[i])
  
  #### If the column index == 0 then there is no matching column in brain cn matrix - proceed to next column ####  
  #### If the column index == 1 then a matched column was found #### 
  if(length(COL.inx) > 0)
  {
    #### Replace NAs with "" ####
    Fusion_Matrix.rearranged[is.na(Fusion_Matrix.rearranged[,COL.inx]),COL.inx] = "" 
    
    #### Find rows with text in matrices ####
    ONCO.TXT.CHECK.inx = nchar(as.character(New.New.onco.matrix[,i])) >= 1
    cn.TXT.CHECK.inx = nchar(as.character(Fusion_Matrix.rearranged[,COL.inx])) >= 1
    
    New.onco.Column = New.New.onco.matrix[,i]
    
    #### Append records with ";" for matching rows ####
    
    if(any(ONCO.TXT.CHECK.inx & cn.TXT.CHECK.inx))
    {
      Common.inx = which(ONCO.TXT.CHECK.inx & cn.TXT.CHECK.inx)
      New.onco.Column[Common.inx] = paste(as.character(New.New.onco.matrix[Common.inx,i]),as.character(Fusion_Matrix.rearranged[Common.inx,COL.inx]),sep=";")
      New.onco.Column[Common.inx] = gsub(" ","",New.onco.Column[Common.inx])
    }#if(any(ONCO.TXT.CHECK.inx & cn.TXT.CHECK.inx))
    
    #### Just add Brain CN info for umatching rows ####
    if(any(!(ONCO.TXT.CHECK.inx & cn.TXT.CHECK.inx)))
    {
      #### Determines index of the rows with no text in onco but text in Brain cn ####
      Single.Insert.inx = which(!ONCO.TXT.CHECK.inx & cn.TXT.CHECK.inx) 
      New.onco.Column[Single.Insert.inx] = as.character(Fusion_Matrix.rearranged[Single.Insert.inx,COL.inx])
    }#if(any(!(ONCO.TXT.CHECK.inx & cn.TXT.CHECK.inx)))
    
    #### Insert new column into Matrix ####
    New.New.onco.matrix[,i] = New.onco.Column
  }#if(length(COL.inx) > 0)
  
}#for(i in 1:ncol(New.onco.matrix))


##shift ID columns over and write table
write.table(as.data.frame(New.New.onco.matrix), 
            paste0(mainDir, subDirHist, "/", broad.hist, "-mut-cn-fus-matrix.txt"), 
            quote = F,col.names = T,row.names = T, sep = "\t")
New.New.onco.matrix[1:10,1:4]





