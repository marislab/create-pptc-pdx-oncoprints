if (!require("data.table")){
  install.packages("https://cran.r-project.org/src/contrib/Archive/data.table/data.table_1.2.tar.gz", repos=NULL, dependencies = TRUE, type = "source")
  library(data.table)
}
if (!require("devtools")){
  install.packages("https://cran.r-project.org/src/contrib/Archive/devtools/devtools_2.0.2.tar.gz", repos=NULL, dependencies = TRUE, type = "source")
  library(devtools)
}
devtools::install_github(repo = "marislab/maftools")
library(maftools)

if (!require("rmatio")){
  install.packages("rmatio", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(rmatio)
}
if (!require("dplyr")){
  install.packages("dplyr", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(dplyr)
}
if (!require("tidyr")){
  install.packages("tidyr", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(tidyr)
}
if (!require("ggplot2")){
  install.packages("ggplot2", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(ggplot2)
}

if (!require("BSgenome.Hsapiens.UCSC.hg19")){
  install.packages("https://bioconductor.org/packages/release/data/annotation/src/contrib/BSgenome.Hsapiens.UCSC.hg19_1.4.0.tar.gz", repo=NULL, type="source", dependencies = TRUE)
library(BSgenome.Hsapiens.UCSC.hg19)
  }

if (!require("deconstructSigs")){
  install.packages("https://cran.r-project.org/src/contrib/deconstructSigs_1.8.0.tar.gz", repo=NULL, type="source", dependencies = TRUE)
library(deconstructSigs)
  }
if (!require("circlize")){
  install.packages("circlize", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(circlize)
}
if (!require("gplots")){
  install.packages("gplots", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(gplots)
}
if (!require("reshape")){
  install.packages("reshape", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(gplots)
}
