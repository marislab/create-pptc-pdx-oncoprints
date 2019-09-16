###### create cohesive oncoprints to regenerate Figures 2 and 3########
#
#     Authors: Jo Lynne Rokita, Alvin Farrel, Khushbu Patel     
#     Updated 2019-09-16
################################################################

# working directory (created with git clone)
mainDir <- "~/create-pptc-pdx-oncoprints/"
mainDir <- "~/Documents/GitHub/create-pptc-pdx-oncoprints/"
dataDir <- paste0(mainDir,"data/")
# set path to your git cloned repo
script.folder <- paste0(mainDir, "R/") 
##create directory for output files
ifelse(!dir.exists(file.path(paste0(mainDir, "onco-out/"))), dir.create(file.path(paste0(mainDir, "onco-out/"))), 
       "Directory exists!")

##set wd
setwd(mainDir)

####Dependencies
source(paste0(script.folder, "install-packages.R"))

library(devtools)
library(maftools)
library(rmatio)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(tidyr)
library(ggplot2)
library(deconstructSigs)
library(data.table)
library(circlize)
library(gplots)
library(reshape)


# create new directories in mainDir
subDir <- paste0(mainDir,"onco-out/")
ifelse(!dir.exists(file.path(subDir)), dir.create(file.path(subDir)), "Directory exists!")


#load file for harmonization of gene IDs
gene.ids <- read.delim(paste0(dataDir,"2019-02-14-Hugo-Symbols-approved.txt"),
                       sep = "\t", as.is = T, header = T)
##load clinical file
clin <- read.delim(paste0(dataDir, "pptc-pdx-clinical-web.txt"), as.is = T, header = T)

##specify histology categorizations
broad.hists <- as.list(unique(clin$Histology.Onco.New)) ## use for generation of mutation and CN matrices

###load color functions
source(paste0(script.folder, "mutation-color-function.R"))
source(paste0(script.folder, "demog-color-function.R"))

###load MAF file into WD and into maftools
load(paste0(dataDir,"2019-02-14-allpdx-clean-maf-240.rda"), verbose = T)

###reformat clinical file
source(paste0(script.folder, "reformat-clinical.R"))

###load RNA expression matrix
load(paste0(dataDir,"2019-02-14-PPTC_FPKM_matrix_withModelID-244.rda"), verbose = T) 

###create mutational signatures burden matrix
source(paste0(script.folder, "create-mut-sigs-matrix.R"))

###focal CN matrix
focal.cn.mat <- read.delim(paste0(dataDir, "short_cn_matrix_fpkm1.txt"),as.is=TRUE,check.names=FALSE)
source(paste0(script.folder, "reformat-CN-matrix-to-df.R"))

###add fusions to MAF
source(paste0(script.folder, "reformat-fusion-for-maf-revision.R"))
###read old and new, combined maf
source(paste0(script.folder, "load-maf-maftools.R"))


###load gene of interest list
for (broad.hist in broad.hists){
  try({
    prefix_name <- broad.hist
    print(paste0("creating ", broad.hist, " matrices and oncoprints"))
    ##create directory for results
    subDirHist <- paste0(subDir, broad.hist)
    #dir.create(file.path(subDir, broad.hist))
    ifelse(!dir.exists(file.path(subDir, broad.hist)), dir.create(file.path(subDir, broad.hist)), 
           "Directory exists!")
    
    ##plot oncoprints
    source(paste0(script.folder, "create-complexheat-oncoprint-revision.R"))
    
    ###plot dx/relapse co-oncoplots
    source(paste0(script.folder, "co-oncoplots.R"))
    
    log("a")}, ##suppress error message with oncoprint plotting so loop continues
    silent=TRUE)
  }

##write session info
sink(paste0(subDir,Sys.Date(), "sessionInfo.txt"))
sessionInfo()
sink()
