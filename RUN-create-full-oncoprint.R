###### create cohesive oncoprints to regenerate Figure 2########
#
#     Authors: Jo Lynne Rokita, Alvin Farrel, Khushbu Patel     
#
################################################################

# Set working directory
mainDir <- "~/pptc-pdx-oncoprints/"
script.folder <- "~/create-pptc-pdx-oncoprints/" # path to your git cloned repo
setwd(mainDir)
dataDir <- "~/pptc-pdx-oncoprints/data/"

####Dependencies
source(paste0(script.folder, "install-packages.R"))

library(devtools)
library(maftools)
library(rmatio)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(deconstructSigs)
library(data.table)
library(circlize)
library(gplots)




# create new directories in mainDir
subDir <- paste0(mainDir,"onco-out/")
ifelse(!dir.exists(file.path(subDir)), dir.create(file.path(subDir)), "Directory exists!")


#load file for harmonization of gene IDs
gene.ids <- read.delim(paste0(dataDir,"2019-02-14-Hugo-Symbols-approved.txt"),
                       sep = "\t", as.is = T, header = T)
##load clinical file
clin <- read.delim(paste0(dataDir, "pptc-pdx-clinical-web.txt"), as.is = T, header = T)

##specify histology categorizations
broad.hists <- as.list(unique(clin$Histology.Oncoprints)) ## use for generation of mutation and CN matrices

###load color functions
source(paste0(script.folder, "mutation-color-function.R"))
source(paste0(script.folder, "demog-color-function.R"))

###load MAF file into WD and into maftools
load(paste0(dataDir,"2019-02-14-allpdx-clean-maf-240.rda"), verbose = T)

#pptc.merge <- rna.maf
source(paste0(script.folder, "load-maf-maftools.R"))

###load RNA expression matrix
load(paste0(dataDir,"2019-02-14-PPTC_FPKM_matrix_withModelID-244.rda"), verbose = T) 

###create mutational signatures burden matrix
source(paste0(script.folder, "create-mut-sigs-matrix.R"))

###focal CN matrix
focal.cn.mat <- read.delim(paste0(dataDir, "short_cn_matrix_fpkm1.txt"),as.is=TRUE,check.names=FALSE)

#### Read fusion file ####
source(paste0(script.folder, "reformat-fusion-as-matrix.R"))


###load gene of interest list
for (broad.hist in broad.hists){
    print(paste0("creating ", broad.hist, " matrices and oncoprints"))
    ##create directory for results
    subDirHist <- paste0(subDir, broad.hist)
    #dir.create(file.path(subDir, broad.hist))
    ifelse(!dir.exists(file.path(subDir, broad.hist)), dir.create(file.path(subDir, broad.hist)), 
           "Directory exists!")
    ##read in gene list
    goi.list <- read.delim(paste0(dataDir, broad.hist, "-goi-list.txt"), sep = "\t",
                       header = F, as.is = T)

    ###use maftools for appropriate oncoprint matrix
    source(paste0(script.folder, "create-mut-matrices.R"))
  
    ###create CN matrices
    source(paste0(script.folder, "create-CN-matrices.R"))
  
    ###merge mut and CN matrices
    source(paste0(script.folder, "merge-mut-CN-matrices.R"))
  
    ###merge mut/CN and fusion matrices
    source(paste0(script.folder, "merge-mut-CN-fusion-matrices.R"))
    
    ##plot oncoprints
    source(paste0(script.folder, "create-complexheat-oncoprint-all.R"))
    
}

##write session info
sink(paste0(subDir,Sys.Date(), "sessionInfo.txt"))
sessionInfo()
sink()
