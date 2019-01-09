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

###create directories for saving files
mainDir <- "~/Box Sync/PPTC-genomics-collaboration/Manuscript/scripts/"
subDir <- "onco-out/"
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), "Directory exists!")

#set directories for saving files, specify histology of interest
pptc.folder <- "~/Box Sync/PPTC-genomics-collaboration/"
script.folder <- "~/Box Sync/PPTC-genomics-collaboration/Manuscript/scripts/oncoprint-r-scripts/"

#load file for harmonization of gene IDs
gene.ids <- read.delim(paste0(pptc.folder, "Data/Hugo_Symbols/2019-01-03-Hugo-Symbols-edited.txt"),
                       sep = "\t", as.is = T, header = T)
##load clinical file
clin <- read.delim(paste0(pptc.folder, "Data/clinical/2018-12-28-pdx-clinical-final-for-paper.txt"), as.is = T, header = T)
clin.pptc <- subset(clin, Model.Part.of.PPTC == "yes")


##specify histology categorizations
broad.hists <- as.list(unique(clin$Histology.Oncoprints)) ## use for generation of mutation and CN matrices

###load color functions
source(paste0(pptc.folder, "Manuscript/figures/oncoprints/mutation-color-function.R"))
source(paste0(pptc.folder, "Manuscript/figures/oncoprints/demog-color-function.R"))

###load MAF file into WD and into maftools
load("~/Box Sync/PPTC-genomics-collaboration/Pedcbio-upload/2019-01-03-allpdx-clean-maf-241.rda", verbose = T)
source(paste0(script.folder, "load-maf-maftools.R"))

###create mutational signatures burden matrix
source(paste0(script.folder, "create-mut-sigs-matrix.R"))

###use GISTIC output to generate copy number matrix
gistic.out <- read.delim(paste0(pptc.folder, "Data/GISTIC-results/all-pdx/2018-08-09-gistic-results-256pdx-noXY-snpfast2-nomirna/all_thresholded.by_genes.txt"),as.is=TRUE,check.names=FALSE)

#### Read fusion file ####
source(paste0(script.folder, "reformat-fusion-as-matrix.R"))

###load gene of interest list
for (broad.hist in broad.hists){
    print(broad.hist)
    ##create directory for results
    subDirHist <- paste0(subDir, broad.hist)
    ifelse(!dir.exists(file.path(mainDir, subDirHist)), dir.create(file.path(mainDir, subDirHist)), "Directory exists!")
    
    ##read in gene list
    goi.list <- read.delim(paste0(pptc.folder, "Data/genelists/", broad.hist, "-goi-list.txt"), sep = "\t",
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


