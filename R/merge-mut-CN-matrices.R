#### Read Onco File ####
onco.matrix = read.delim(paste0(subDirHist, "/", broad.hist, "_onco_matrix.txt"), 
                         sep = "\t", header = T, as.is = T, check.names = F)

colnames(onco.matrix) %in% colnames(final.cn.mat)
dim(onco.matrix)
dim(final.cn.mat)
#### Rearranges CN file into same order as the onco file ####
#### Removes Genes that aren't present in Onco file ####
cn.matrix.rearranged=final.cn.mat[match(rownames(onco.matrix),rownames(final.cn.mat)),]

#### Determine Genes in cn but not in onco ####
Unmatched_Genes = setdiff(rownames(final.cn.mat),rownames(onco.matrix))

#### Create matrix with unmatched brain cn genes to be appended to rearranged matrix ####
Unmatched_cn_Gene_matrix = final.cn.mat[rownames(final.cn.mat) %in% Unmatched_Genes,]

#### Create blank matrix of unmatched genes to be appeneded to onco matrix ####
Onco_Unmatched_Gene_append.matrix = matrix(nrow=nrow(Unmatched_cn_Gene_matrix),ncol=ncol(onco.matrix))
rownames(Onco_Unmatched_Gene_append.matrix) = rownames(Unmatched_cn_Gene_matrix)
colnames(Onco_Unmatched_Gene_append.matrix) = colnames(onco.matrix)

#### Append unmatched genes matrix to Brain cn matrix ####
cn.matrix.rearranged = rbind(cn.matrix.rearranged,Unmatched_cn_Gene_matrix)

#### Make new matrix coverting Factors to Characters and add rows for unmatched Brain cn genes ####
New.onco.matrix = rbind(as.matrix(onco.matrix),Onco_Unmatched_Gene_append.matrix)


#### Replace NAs with "" ####
New.onco.matrix[is.na(New.onco.matrix)] = ""

#### Check Each onco column name for a corresponding column name in the Brain CN file ####
for(i in 1:ncol(onco.matrix))
{
  
  COL.inx = which(colnames(cn.matrix.rearranged) == colnames(onco.matrix)[i])
  
  #### If the column index == 0 then there is no matching column in brain cn matrix - proceed to next column ####  
  #### If the column index == 1 then a matched column was found #### 
  if(length(COL.inx) > 0)
  {
    #### Replace NAs with "" ####
    cn.matrix.rearranged[is.na(cn.matrix.rearranged[,COL.inx]),COL.inx] = "" 
    
    #### Find rows with text in matrices ####
    ONCO.TXT.CHECK.inx = nchar(as.character(New.onco.matrix[,i])) >= 1
    cn.TXT.CHECK.inx = nchar(as.character(cn.matrix.rearranged[,COL.inx])) >= 1
    
    New.onco.Column = New.onco.matrix[,i]
    
    #### Append records with ";" for matching rows ####
    
    if(any(ONCO.TXT.CHECK.inx & cn.TXT.CHECK.inx))
    {
      Common.inx = which(ONCO.TXT.CHECK.inx & cn.TXT.CHECK.inx)
      New.onco.Column[Common.inx] = paste(as.character(New.onco.matrix[Common.inx,i]),as.character(cn.matrix.rearranged[Common.inx,COL.inx]),sep=";")
      New.onco.Column[Common.inx] = gsub(" ","",New.onco.Column[Common.inx])
    }#if(any(ONCO.TXT.CHECK.inx & cn.TXT.CHECK.inx))
    
    #### Just add Brain CN info for umatching rows ####
    if(any(!(ONCO.TXT.CHECK.inx & cn.TXT.CHECK.inx)))
    {
      #### Determines index of the rows with no text in onco but text in Brain cn ####
      Single.Insert.inx = which(!ONCO.TXT.CHECK.inx & cn.TXT.CHECK.inx) 
      New.onco.Column[Single.Insert.inx] = as.character(cn.matrix.rearranged[Single.Insert.inx,COL.inx])
    }#if(any(!(ONCO.TXT.CHECK.inx & cn.TXT.CHECK.inx)))
    
    #### Insert new column into Matrix ####
    New.onco.matrix[,i] = New.onco.Column
  }#if(length(COL.inx) > 0)
  
}#for(i in 1:ncol(onco.matrix))

###remove empty rows
New.onco.matrix <- New.onco.matrix[!apply(New.onco.matrix == "", 1, all),]  

##shift ID columns over and write table
#write.table(as.data.frame(New.onco.matrix), paste0(mainDir, subDirHist, "/", broad.hist,
 #                                   "_full_matrix.txt"), quote = F,col.names = T,row.names = T, sep = "\t")


