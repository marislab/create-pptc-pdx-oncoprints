###create cohesive oncoprints 

####Dependencies
#devtools::install_github(repo = "jharenza/maftools")
#install.extras('NMF')
require(maftools)
require(NMF)
require(rmatio)
require(BSgenome.Hsapiens.UCSC.hg19)
require(dplyr)
require(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(deconstructSigs)
library(data.table)

# Setting working directory
setwd("~")
mainDir <- "~/pptc-pdx-oncoprints/"
dataDir <- "~/pptc-pdx-oncoprints/data/"
script.folder <- "~/create-pptc-pdx-oncoprints/"

# create new directories in mainDir
dir.create(file.path(mainDir,"onco-out"))
subDir <- paste0(mainDir,"onco-out/")

###create directories for saving files
#mainDir <- "~/Box Sync/PPTC-genomics-collaboration/Manuscript/scripts/"
#subDir <- "onco-out/"
#subDir <- "test/"

ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), "Directory exists!")

#set directories for saving files, specify histology of interest
#pptc.folder <- "~/Box Sync/PPTC-genomics-collaboration/"
#script.folder <- "~/Box Sync/PPTC-genomics-collaboration/Manuscript/scripts/oncoprint-r-scripts/"

#load file for harmonization of gene IDs
gene.ids <- read.delim(paste0(dataDir,"2019-02-14-Hugo-Symbols-approved.txt"),
                       sep = "\t", as.is = T, header = T)
##load clinical file
clin <- read.delim(paste0(dataDir, "pptc-pdx-clinical-web.txt"), as.is = T, header = T)
clin.pptc <- subset(clin, Model.Part.of.PPTC == "yes")
#write.table(clin.pptc, paste0(pptc.folder, "Data/clinical/2019-02-09-pdx-clinical-final-for-paper.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

##specify histology categorizations
broad.hists <- as.list(unique(clin.pptc$Histology.Oncoprints)) ## use for generation of mutation and CN matrices

###load color functions
source(paste0(script.folder, "mutation-color-function.R"))
source(paste0(script.folder, "demog-color-function.R"))

###load MAF file into WD and into maftools
load(paste0(dataDir,"2019-02-14-allpdx-clean-maf-240.rda"), verbose = T)

sort(unique(pptc.merge$Tumor_Sample_Barcode))

#pptc.merge <- rna.maf
source(paste0(script.folder, "load-maf-maftools.R"))

###load RNA expression matrix
load(paste0(dataDir,"2019-02-14-PPTC_FPKM_matrix_withModelID-244.rda"), verbose = T) 

###create mutational signatures burden matrix
source(paste0(script.folder, "create-mut-sigs-matrix.R"))

###use GISTIC output for hemizygous deletions
#gistic.out <- read.delim(paste0(pptc.folder, "Data/GISTIC-results/all-pdx/2018-08-09-gistic-results-256pdx-noXY-snpfast2-nomirna/all_thresholded.by_genes.txt"),as.is=TRUE,check.names=FALSE)
#colnames(gistic.out)[colnames(gistic.out) == "IC-2264PNET"] <- "IC-2664PNET"
###fix IC-2664

###focal CN matrix
focal.cn.mat <- read.delim(paste0(dataDir, "short_cn_matrix_fpkm1.txt"),as.is=TRUE,check.names=FALSE)

#### Read fusion file ####
source(paste0(script.folder, "reformat-fusion-as-matrix.R"))

###load gene of interest list
for (broad.hist in broad.hists){
    print(paste0("creating ", broad.hist, " matrices and oncoprints"))
    ##create directory for results
    subDirHist <- paste0(subDir, broad.hist)
    ifelse(!dir.exists(file.path(mainDir, subDirHist)), dir.create(file.path(mainDir, subDirHist)), "Directory exists!")
    
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
  
    ###plot oncoprint
    ifelse(broad.hist == "neuroblastoma" | broad.hist == "osteosarcoma" | broad.hist == "renal", 
    source(paste0(script.folder, "create-complexheat-oncoprint-", broad.hist, ".R")), 
           source(paste0(script.folder, "create-complexheat-oncoprint-other.R"))
  )
}

##write session info
sink(paste0(subDir,Sys.Date(), "sessionInfo.txt"))
sessionInfo()
sink()


