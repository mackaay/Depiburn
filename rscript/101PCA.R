setwd("./DEpiBurn_EMseq/")
library(GenomicRanges)
library(RColorBrewer)
library(data.table)
library(edgeR)
library(limma)
library(stringr)
library(ggplot2)
library(ggpubr)


pal.graft <- brewer.pal(8, "Set1")[4:5]
pal.depth <- c("#999999", "#E69F00", "#56B4E9")
pal.coll <- brewer.pal(5, "Set3")
pal.epi <- brewer.pal(4, "Reds")
pal.cell <- brewer.pal(8,"Set3")
pal.group <- c( "#666666" , "#FFFFB3",  "darkred", "#BEBADA","#8DD3C7"  ,"#80B1D3")
names(pal.group) <- c("BlisterFluid", "SA", "Plasma", "Nerve", "Skin", "WholeBlood")

#blood plasma samples
dir <- "/datasets/work/hb-meth-atlas/work/Data/level_2/public/cfDNA_aging_PRJNA1080627/"
bedgraph.files <- list.files(dir, pattern = "*_sd_CpG.bedGraph.gz", full.names = T, recursive = T)
plasma.df <- NULL
tmp <- read.delim(bedgraph.files[1], stringsAsFactors = F, skip = 1, header = F)
print(bedgraph.files[1])
rownames(tmp) <- paste(tmp$V1, tmp$V2+1, sep = "_")
plasma.df <- cbind(plasma.df, tmp[,"V4"])
rownames(plasma.df) <- rownames(tmp)
plasma.df <- as.data.frame(plasma.df)
plasma.df$V1 <- plasma.df$V1/100
for (i in 2:length(bedgraph.files)) {
  print(i)
  tmp <- read.delim(bedgraph.files[i], stringsAsFactors = F, skip = 1, header = F)
  print(bedgraph.files[i])
  rownames(tmp) <- paste(tmp$V1, tmp$V2+1, sep = "_")
  keep <- intersect(rownames(plasma.df), rownames(tmp))
  tmp <- tmp[keep,]
  plasma.df <- plasma.df[keep,]
  tmp$V4 <- tmp$V4/100
  plasma.df <- cbind(plasma.df, tmp$V4)
  rownames(plasma.df) <- keep
  print(nrow(plasma.df))
  colnames(plasma.df)[i] <- bedgraph.files[i]
  print("Merge Done")
}
colnames(plasma.df) <- basename(bedgraph.files)

pdata <- data.frame(id = colnames(plasma.df) , group = c(rep("MiddleAge",10), rep("Old", 14), rep("Young", 11)))
save(plasma.df, file = "./data/plasma_WGBS_CpGl.rds")



#Blister fluid samples
dir <- "./data/"
bedgraph.files <- list.files(dir, pattern = "*_sd_CpG.bedGraph.gz", full.names = T)
id.name <- gsub("_", "",substr(bedgraph.files, start = 32L, stop = 38L) )
phenoData <- read.csv("./SampleMeta_20220222.csv", stringsAsFactors = F)
rownames(phenoData) <- paste0("QUT", phenoData$SAMPLE_ID)
phenoData <- phenoData[id.name,]
phenoData$DAYS.TO.RE.EP <- gsub("\\+", "", phenoData$DAYS.TO.RE.EP)
phenoData$DAYS.TO.RE.EP <- as.numeric(phenoData$DAYS.TO.RE.EP)
names(pal.graft) <- levels(as.factor(phenoData$graft))
names(pal.depth) <- levels(factor(phenoData$DEPTH, levels = c("superficial", "deep", "full_thickness")))
names(pal.epi) <- levels(factor(phenoData$reepi))

blister.df <- NULL
tmp <- read.delim(bedgraph.files[1], stringsAsFactors = F, skip = 1, header = F)
print(bedgraph.files[1])
rownames(tmp) <- paste(tmp$V1, tmp$V2+1, sep = "_")
blister.df <- cbind(blister.df, tmp[,"V4"])
rownames(blister.df) <- rownames(tmp)
blister.df <- as.data.frame(blister.df)
blister.df$V1 <- blister.df$V1/100
for (i in 2:length(bedgraph.files)) {
  print(i)
  tmp <- read.delim(bedgraph.files[i], stringsAsFactors = F, skip = 1, header = F)
  print(bedgraph.files[i])
  rownames(tmp) <- paste(tmp$V1, tmp$V2+1, sep = "_")
  keep <- intersect(rownames(blister.df), rownames(tmp))
  tmp <- tmp[keep,]
  blister.df <- blister.df[keep,]
  tmp$V4 <- tmp$V4/100
  blister.df <- cbind(blister.df, tmp$V4)
  rownames(blister.df) <- keep
  print(nrow(blister.df))
  colnames(blister.df)[i] <- bedgraph.files[i]
  print("Merge Done")
}
colnames(blister.df) <- gsub("_", "",substr(bedgraph.files, start = 32L, stop = 38L) )

save(blister.df , file ="./data/blister_EMseq_CpG_all.rds")



#Skin, Nerve, SA
dir <- "/datasets/work/hb-meth-atlas/work/Data/level_2/"
bedgraph.files <- list.files(dir, pattern = "*_sd_CpG.bedGraph.gz", full.names = T, recursive = T)
bedgraph.files <- bedgraph.files[grep("Skin|Nerve|Subcutaneous|WholeBlood", bedgraph.files)]
bedgraph.files <- bedgraph.files[grep("Blister", bedgraph.files, invert = T)]

tissue.df <- NULL
tmp <- read.delim(bedgraph.files[1], stringsAsFactors = F, skip = 1, header = F)
print(bedgraph.files[1])
rownames(tmp) <- paste(tmp$V1, tmp$V2+1, sep = "_")
tissue.df <- cbind(tissue.df, tmp[,"V4"])
rownames(tissue.df) <- rownames(tmp)
tissue.df <- as.data.frame(tissue.df)
tissue.df$V1 <- tissue.df$V1/100
for (i in 2:length(bedgraph.files)) {
  print(i)
  tmp <- read.delim(bedgraph.files[i], stringsAsFactors = F, skip = 1, header = F)
  print(bedgraph.files[i])
  rownames(tmp) <- paste(tmp$V1, tmp$V2+1, sep = "_")
  keep <- intersect(rownames(tissue.df), rownames(tmp))
  tmp <- tmp[keep,]
  tissue.df <- tissue.df[keep,]
  tmp$V4 <- tmp$V4/100
  tissue.df <- cbind(tissue.df, tmp$V4)
  rownames(tissue.df) <- keep
  print(nrow(tissue.df))
  colnames(tissue.df)[i] <- bedgraph.files[i]
  print("Merge Done")
}
colnames(tissue.df) <- basename(bedgraph.files)
colnames(tissue.df) <- gsub("_sd_CpG.bedGraph.gz", "", colnames(tissue.df))

save(tissue.df, file = "./data/tissue_WGBS_CpG.rds")






##global, remove ChrM , Sexual Chr and Chr Unknown  ###########
blister.df <- blister.df[grep("chrM|chrX|chrY|Un", rownames(blister.df), invert= T),]
blister.df.Mvalue <- log2(blister.df+0.0000001/(1-blister.df))
blister.df.Mvalue[is.infinite(blister.df.Mvalue)] <- 23
plotMDS(blister.df.Mvalue, top=2000, gene.selection="common", 
        col=pal.depth[phenoData$DEPTH], dim=c(1,2), pch = 19)
plotMDS(blister.df.Mvalue, top=2000, gene.selection="common", 
        col=pal.depth[phenoData$DEPTH], dim=c(1,2))
legend("topright", legend=levels(factor(phenoData$DEPTH, levels = c("superficial", "deep", "full_thickness"))), 
       text.col=pal.depth, cex=0.7, bg="white")


plasma.df <- plasma.df[grep("chrM|chrX|chrY|Un", rownames(plasma.df), invert= T),]
plasma.df.Mvalue <- log2(plasma.df+0.0000001/(1-plasma.df))
plasma.df.Mvalue[is.infinite(plasma.df.Mvalue)] <- 23

keep <- intersect(rownames(blister.df.Mvalue), rownames(plasma.df.Mvalue))
mVal <- cbind(plasma.df.Mvalue[keep,] , blister.df.Mvalue[keep,])
pdata <- data.frame(id = colnames(plasma.df) , group = "Plasma")
pdata <- rbind(pdata, data.frame(id = colnames(blister.df), group = "BlisterFluid"))

tissue.df <- tissue.df[grep("chrM|chrX|chrY|Un", rownames(tissue.df), invert= T),]
tissue.df.Mvalue <- log2(tissue.df+0.0000001/(1-tissue.df))
tissue.df.Mvalue[is.infinite(tissue.df.Mvalue)] <- 23

keep <- intersect(rownames(tissue.df.Mvalue), rownames(mVal))
mVal <- cbind(mVal[keep,] , tissue.df.Mvalue[keep,])
pdata <- rbind(pdata, data.frame(id = colnames(tissue.df.Mvalue), group = c(rep("Nerve", 3), rep("SA", 3), rep("WholeBlood", 24), 
                                                                            rep("Nerve",2), rep("Skin",2), rep("WholeBlood", 5), 
                                                                            rep("Skin", 22))))
plotMDS(mVal, top=2000, gene.selection="common", 
        col=pal.group[factor(pdata$group)], dim=c(1,2), pch = 19)
legend("topright", legend=levels(factor(pdata$group)), text.col=pal.group,
       bg="white", cex=0.7)


library(umap)
set.seed(12345) # Ensure reproducibility
umap_results <- umap(t(mVal))
plot_df <- data.frame(umap_results$layout, condition = pdata$group)
ggplot(plot_df, aes(X1, X2, color = condition)) + geom_point() + 
  scale_color_manual(values = pal.group)+
  labs(x = "UMAP1",                                              # labelling x axis
       y = "UMAP2",                                        # labeling y axis
       title = "",        # title
       fill = "Cell Type") +                               # legend
  theme(
    axis.text.x = element_text(angle = 90,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )

##CpG island ############
cgi <- read.delim("/datasets/work/hb-diab-cfdna/work/Data/annotation/cpgisland.hg38.bed", 
                  stringsAsFactors = F)
#cgi$chromStart <- cgi$chromStart -2000
#cgi$chromEnd <- cgi$chromEnd +2000
cgi <- cgi[grep("_", cgi$chrom, invert = T),]
cgi <- cgi[grep("chrX|Y", cgi$chrom, invert = T),]
cgi.gr <- makeGRangesFromDataFrame(cgi, keep.extra.columns = F)

blister.gr <- matrix(unlist(strsplit(rownames(blister.df), split = "_") ), ncol = 2, byrow =T)
blister.gr <- as.data.frame(blister.gr)
blister.gr$V2 <- as.integer(blister.gr$V2)
blister.gr$start <- blister.gr$V2 - 1
blister.gr$end <- blister.gr$V2 + 1
colnames(blister.gr)[1] <- "chr"
blister.gr <- makeGRangesFromDataFrame(blister.gr, keep.extra.columns = F)
tmp.ov <- subsetByOverlaps(blister.gr, cgi.gr) 
tmp.ov <- as.data.frame(tmp.ov)
rownames(tmp.ov) <- paste(tmp.ov$seqnames, (tmp.ov$start + 1), sep = "_")

keep <- intersect(rownames(tmp.ov), rownames(blister.df))
blister.cgi <- blister.df[keep,]

blister.cgi.Mvalue <- log2(blister.cgi+0.0000001/(1-blister.cgi))
blister.cgi.Mvalue[is.infinite(blister.cgi.Mvalue)] <- 23

plotMDS(blister.cgi.Mvalue, top=2000, gene.selection="common", 
        col=pal.depth[phenoData$DEPTH], dim=c(1,2), pch = 19)
plotMDS(blister.cgi.Mvalue, top=2000, gene.selection="common", 
        col=pal.depth[phenoData$DEPTH], dim=c(1,2))
plotMDS(blister.cgi.Mvalue, top=2000, gene.selection="common", 
        col=pal.epi[phenoData$reepi], dim=c(1,2), pch = 19)




mVal.gr <- matrix(unlist(strsplit(rownames(mVal), split = "_") ), ncol = 2, byrow =T)
mVal.gr <- as.data.frame(mVal.gr)
mVal.gr$V2 <- as.integer(mVal.gr$V2)
mVal.gr$start <- mVal.gr$V2 - 1
mVal.gr$end <- mVal.gr$V2 + 1
colnames(mVal.gr)[1] <- "chr"
mVal.gr <- makeGRangesFromDataFrame(mVal.gr, keep.extra.columns = F)
tmp.ov <- subsetByOverlaps(mVal.gr, cgi.gr) 
tmp.ov <- as.data.frame(tmp.ov)
rownames(tmp.ov) <- paste(tmp.ov$seqnames, (tmp.ov$start + 1), sep = "_")

keep <- intersect(rownames(tmp.ov), rownames(mVal))
mVal.cgi <- mVal[keep,]





##nerve########## 
nerve <- rbind(tissue.marker.hyper[grep("nerve", tissue.marker.hyper$tissue),] , tissue.marker.un[grep("nerve", tissue.marker.un$tissue),])

marker.gr <- matrix(unlist(strsplit(nerve$loci, split = "_") ), ncol = 2, byrow =T)
marker.gr <- as.data.frame(marker.gr)
marker.gr$V2 <- as.integer(marker.gr$V2)
marker.gr$start <- marker.gr$V2 - 1
marker.gr$end <- marker.gr$V2 + 1
colnames(marker.gr)[1] <- "chr"
marker.gr <- makeGRangesFromDataFrame(marker.gr, keep.extra.columns = F)
tmp.ov <- subsetByOverlaps(blister.gr, marker.gr) 
tmp.ov <- as.data.frame(tmp.ov)
rownames(tmp.ov) <- paste(tmp.ov$seqnames, (tmp.ov$start + 1), sep = "_")

keep <- intersect(rownames(tmp.ov), rownames(blister.df))
blister.marker <- blister.df[keep,]

blister.marker.Mvalue <- log2(blister.marker+0.0000001/(1-blister.marker))
blister.marker.Mvalue[is.infinite(blister.marker.Mvalue)] <- 23

plotMDS(blister.marker.Mvalue, top=2000, gene.selection="common", 
        col=pal.depth[phenoData$DEPTH], dim=c(1,2), pch = 19)
plotMDS(blister.marker.Mvalue, top=2000, gene.selection="common", 
        col=pal.depth[phenoData$DEPTH], dim=c(1,2))
plotMDS(blister.marker.Mvalue, top=2000, gene.selection="common", 
        col=pal.epi[phenoData$reepi], dim=c(1,2), pch = 19)


library(umap)
set.seed(12345) # Ensure reproducibility
umap_results <- umap(t(blister.marker))
plot_df <- data.frame(umap_results$layout, condition = phenoData$DEPTH)
ggplot(plot_df, aes(X1, X2, color = condition)) + geom_point() + theme_minimal()



