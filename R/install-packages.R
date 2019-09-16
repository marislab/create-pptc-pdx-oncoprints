if (!require("devtools")){
  install.packages("devtools", repos='http://cran.us.r-project.org', dependencies = TRUE)
}
devtools::install_github(repo = "marislab/maftools")

if (!require("rmatio")){
  install.packages("rmatio", repos='http://cran.us.r-project.org', dependencies = TRUE)
}
if (!require("dplyr")){
  install.packages("dplyr", repos='http://cran.us.r-project.org', dependencies = TRUE)
}
if (!require("tidyr")){
  install.packages("tidyr", repos='http://cran.us.r-project.org', dependencies = TRUE)
}
if (!require("ggplot2")){
  install.packages("ggplot2", repos='http://cran.us.r-project.org', dependencies = TRUE)
}
if (!require("data.table")){
  install.packages("data.table", repos='http://cran.us.r-project.org', dependencies = TRUE)
}
if (!require("BSgenome.Hsapiens.UCSC.hg19")){
  install.packages("https://bioconductor.org/packages/release/data/annotation/src/contrib/BSgenome.Hsapiens.UCSC.hg19_1.4.0.tar.gz", repo=NULL, type="source", dependencies = TRUE)
}
if (!require("ComplexHeatmap")){
  install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/ComplexHeatmap_1.20.0.tar.gz", repo=NULL, type="source",dependencies = TRUE)
}
if (!require("deconstructSigs")){
  install.packages("https://cran.r-project.org/src/contrib/deconstructSigs_1.8.0.tar.gz", repo=NULL, type="source", dependencies = TRUE)
}
if (!require("circlize")){
  install.packages("circlize", repos='http://cran.us.r-project.org', dependencies = TRUE)
}
if (!require("gplots")){
  install.packages("gplots", repos='http://cran.us.r-project.org', dependencies = TRUE)
}
