library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(gridExtra)

pptc.folder <- "~/Box Sync/PPTC-genomics-collaboration/"
git.folder <- "~/Documents/GitHub/pdx-classification/"
git.folder2 <- "~/Documents/GitHub/create-pptc-pdx-oncoprints/"

source(paste0(git.folder2, "mutation-color-function.R"))
source(paste0(pptc.folder, "Manuscript/scripts/theme.R"))

##load and merge with clinical file
clin = read.delim(paste0(pptc.folder, "Data/clinical/2019-07-05-clin.txt"), sep = "\t", header = T, as.is = T)

##read in scores
scores <- read.delim(paste0(git.folder, "results/classifier_scores.tsv"),
                     sep = "\t", header = T, as.is = T)
colnames(scores)[colnames(scores) == "sample_id"] <- "Model"
scores.clin <- merge(scores, clin[,c("Model", "Histology.Detailed2")], all.x = T)
write.table(scores.clin, "~/Downloads/scores-clin.txt", sep = "\t", col.names = T, row.names = F, quote = F)
###breakpoint correlation:
bkpt <- read.delim(paste0(pptc.folder, "Manuscript/figures/SV-figure/Breakpoints_rawdata.txt"), sep = "\t",
                   header = T, as.is = T)

bp.scores <- merge(bkpt[,c("Model","n.breakpoints")], scores, all = T)
t.scores = setNames(data.frame(t(bp.scores[,-1])), bp.scores[,1])
t.scores$features <- rownames(t.scores)
t.scores <- t.scores[,c(ncol(t.scores), 1:ncol(t.scores)-1)]
keep <- t.scores[c("tp53_score", "tp53_shuffle", 'n.breakpoints'),]
#write.table(keep, "~/Downloads/scores.txt", col.names = T, row.names = F, sep = "\t", quote = F)
att <- scores.clin[,c("Model", "Histology.Detailed")]
#att[1] <- "node"
#write.table(att, "~/Downloads/att.txt", col.names = T, row.names = F, sep = "\t", quote = F)

names(clin.class)
clust.df <- bp.scores[,c("Model", "tp53_score", "n.breakpoints")]
rownames(clust.df) <- clust.df$Model
clust.df$Model <- NULL
clust.df2 <- na.omit(clust.df)
clust.df2 <- scale(clust.df2)
head(clust.df2)
distance <- get_dist(clust.df2)
k2 <- kmeans(clust.df2, centers = 2, nstart = 25)
str(k2)
fviz_cluster(k2, data = clust.df2)
k3 <- kmeans(clust.df2, centers = 3, nstart = 25)
k4 <- kmeans(clust.df2, centers = 4, nstart = 25)
k5 <- kmeans(clust.df2, centers = 5, nstart = 25)

# plots to compare
p1 <- fviz_cluster(k2, geom = "point", data = clust.df2) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point",  data = clust.df2) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = clust.df2) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point",  data = clust.df2) + ggtitle("k = 5")

grid.arrange(p1, p2, p3, p4, nrow = 2)

library(gplots) 
y <- clust.df2 
library(pheatmap); library("RColorBrewer")
pdf("~/Downloads/test.pdf", height = 20, width = 5)
pheatmap(y, color=brewer.pal(9,"Blues"))
hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete")
mycol <- colorpanel(40, "darkblue", "yellow", "white")
heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol,
          scale="row", density.info="none", trace="none", 
          RowSideColors=as.character(mycl))
dev.off()
library(ComplexHeatmap)
Heatmap(clust.df2, name = "mat", row_km = 3, row_gap = unit(5, "mm"))

annot <- as.data.frame(clust.df2)
annot$Model <- rownames(annot)
annot <- merge(annot, clin[,c("Model", "Histology.Detailed")], all.x = T)
annot$tp53_score <- NULL
annot$n.breakpoints <- NULL
color.df <- as.data.frame(histcol)
color.df$Histology.Detailed <- rownames(color.df)
annot.col <- merge(annot, color.df, all.x = T, by = "Histology.Detailed")


ht1 = Heatmap(clust.df2)
row_order(ht1)
ht2 = draw(ht1)
row_order(ht2)
model.order <- annot[row_order(ht2),]
col.order <- annot.col[match(model.order$Model, annot.col$Model), ]
col.order$histcol <- as.character(col.order$histcol)

pdf("~/Downloads/colortest.pdf", width = 20, height = 3)
print(pal(col.order$histcol)) ###not printing the right colors
dev.off()

###annotation
pdf("~/Downloads/test.pdf", height = 30, width = 5)
col_fun = colorRamp2(c(-1, 0, 4), c("deeppink3", "white", "cornflowerblue"))

ht <- Heatmap(clust.df2, name = "mat", row_km = 3, col = col_fun)
draw(ht)
dev.off()

