library(data.table)
library(GenomicRanges)
library(stringr)
library(nnls)
library(RColorBrewer)
setwd("./DEpiBurn_EMseq/")


load(file = "/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/CGI_blister_marker_matrix_0.1_0.5_endo.RData")
rownames(tissue.marker.hyper) <- tissue.marker.hyper$loci
rownames(tissue.marker.un) <- tissue.marker.un$loci
##all meth marker
all.meth.marker <- c(tissue.marker.hyper$loci,tissue.marker.un$loci) 
all.tissue.mat <- rbind(tissue.hyper.mat, tissue.un.mat)
all.meth.marker.df <- data.frame(loci = c(tissue.marker.hyper$loci,tissue.marker.un$loci) , 
                                 tissue = c(tissue.marker.hyper$tissue, tissue.marker.un$tissue), 
                                 marker = c(rep("hyper", nrow(tissue.marker.hyper)), rep("un", nrow(tissue.marker.un))),
                                 stringsAsFactors = F)
rownames(all.meth.marker.df) <- all.meth.marker.df$loci
#reduce.meth.marker.df <- all.meth.marker.df[grep("adipose|macrophage|nerve|neutrophil|skin", all.meth.marker.df$tissue),]
colnames(all.tissue.mat) <- c("Macrophage", "Nerve", "Neutrophil", "SA", "Fibroblast", "Keratinocyte", "EC")


{##In Silico simulation ####
dir <- "/datasets/work/hb-burns-meth/work/chenkai/skin_tissue/"
bedgraph.files <- list.files(dir, pattern = "*_sd_CpG.bedGraph.gz", recursive = T, full.names = T)

silico.df <- c()
for (i in 1:length(bedgraph.files)) {
  print(i)
  tmp <- read.delim(bedgraph.files[i], stringsAsFactors = F, skip = 1, header = F)
  #toi.tmp <- toi.tmp[,c(1:4)]
  print(bedgraph.files[i])
  rownames(tmp) <- paste(tmp$V1, tmp$V2, sep = "_")
  keep<- intersect(rownames(tmp), all.meth.marker)
  tmp.df <- as.matrix(round(tmp[keep,"V4"]/100, 4), nrow = length(keep))
  rownames(tmp.df) <- keep
  if(i < 1.5){
    silico.df <- tmp.df
    silico.df <- as.data.frame(silico.df)
    rownames(silico.df) <- rownames(tmp.df)
    print(nrow(tmp.df))
  }else {
    keep<- intersect(rownames(tmp.df), rownames(silico.df))
    tmp.df <- tmp.df[keep,]
    silico.df <- silico.df[keep,]
    silico.df <- cbind(silico.df, tmp.df)
    print(nrow(silico.df))
  }
  print("merge Done")
}
colnames(silico.df) <-  basename(bedgraph.files)

save(silico.df, file = "/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/silico_df.rds")

}


#Deconvolution for purified tissue
dir <- "/datasets/work/hb-burns-meth/work/chenkai/skin_tissue/"
bedgraph.files <- list.files(dir, pattern = "*_sd_CpG.bedGraph.gz", recursive = T, full.names = T)

decon.df <- c()
for (i in 1:length(bedgraph.files)) {
  print(i)
  tmp <- read.delim(bedgraph.files[i], stringsAsFactors = F, skip = 1, header = F)
  #toi.tmp <- toi.tmp[,c(1:4)]
  print(bedgraph.files[i])
  rownames(tmp) <- paste(tmp$V1, tmp$V2, sep = "_")
  keep<- intersect(rownames(tmp), all.meth.marker)
  tmp.df <- as.matrix(round(tmp[keep,"V4"]/100, 4), nrow = length(keep))
  
  ref.mat <- as.matrix(all.tissue.mat[keep,])
  mod1 <- nnls(ref.mat, tmp.df)
  decon <- mod1$x/sum(mod1$x)
  #colnames(tissue.hyper.mat)
  decon.df <- cbind(decon.df, decon)
  print("merge Done")
}
rownames(decon.df) <- colnames(all.tissue.mat)
decon.all.df <- decon.df
colnames(decon.all.df)  <- sub("\\..*", "", basename(bedgraph.files))
colnames(decon.all.df)  <- sub("_sd_CpG", "", colnames(decon.all.df)  )
colnames(decon.all.df)  <- sub("\\.bedGraph.gz", "", colnames(decon.all.df)  )

save(decon.all.df, file = "/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/insilico_tissue_decon.rds")


library(reshape2)
insilico.df <- decon.all.df[,grep("CellLine|Muscle", colnames(decon.all.df), invert = T)]
insilico.df <- melt(insilico.df)
as.data.frame(matrix(unlist(str_split(insilico.df$Var2, pattern = "_")) , ncol = 4, byrow = T))
insilico.df <- cbind(insilico.df, as.data.frame(matrix(unlist(str_split(insilico.df$Var2, pattern = "_")) , ncol = 4, byrow = T)))
colnames(insilico.df) <- c("CellType", "id", "perc", "tissue", "subtissue", "status", "name")

pal.cell <- brewer.pal(8,"Set3")
names(pal.cell) <- c("Keratinocyte", "Nerve","SA", "EC","Neutrophil","Fibroblast", "Macrophage")


library(ggplot2)
library(ggpubr)
insilico.df<- insilico.df[grep("Skin", insilico.df$tissue),] 
insilico.df <- insilico.df[grep("F51|M37", insilico.df$name, invert = T),]
insilico.df %>%
  ggplot(data = ., mapping = aes(x = name, y = perc, fill = CellType)) +
  geom_col() +
  #geom_text(mapping = aes(label = paste0(Percentage*100, "%") ),              # converting the values to percent
  #          size = 4,                                             # size of the font
  #          position = position_stack(vjust = 0.5)) +             # positioning in the middle
  #scale_fill_brewer(palette = "Accent") +                           # coloring the plot
  scale_fill_manual(values = pal.cell)+
  #scale_fill_manual(values = c("#0FCC39", "#0FCC98", "#0FC5CC", "#0F96CC", "#0F66CC" ,"#0F37CC" ,"#160FCC" ,"#450FCC", "#CC0FA5", "#CC0F46", "#CC360F", "darkgrey", "orange", "darkorange"))+
  #facet_grid(.~Depth) +
  labs(x = "ID",                                              # labelling x axis
       y = "Percentage",                                        # labeling y axis
       title = "",        # title
       fill = "Cell") +                               # legend
  scale_y_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 90,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  ) #+ 
#geom_hline(yintercept=0.43, linetype="dashed",  color = "darkred", size=2)



skin.df <- t(decon.all.df[,grep("CellLine|Muscle", colnames(decon.all.df), invert = T)])
skin.df <- cbind(skin.df, as.data.frame(matrix(unlist(str_split(rownames(skin.df), pattern = "_")) , ncol = 4, byrow = T)))
colnames(skin.df)[8:11] <- c("tissue", "subtissue", "status", "name")
skin.df <- skin.df[grep("Skin", skin.df$tissue),]
skin.df <- skin.df[grep("F51|M37", skin.df$name, invert = T),]
skin.df$group <- c( rep("Young", 6), rep("Old",6), "Old", "Young", rep("Old", 3), rep("Young", 3))
skin.df$Cohort <- "Study1"
skin.df$Cohort[grep("Epidermis", skin.df$subtissue)] <- "Study2"
skin.df$Cohort[grep("NOS", skin.df$subtissue)] <- "Study3"

pal.young <- c("#FDBF6F", "#FF7F00")
names(pal.young) <- c("Young", "Old")

ggplot(skin.df , aes(x=group, y=Nerve, fill = group, label = name)) +
  geom_boxplot()+
  geom_jitter(aes(colour = Cohort))+
  #geom_text(aes(colour = subtissue), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=pal.young)+
  scale_colour_manual(values=c("#A6CEE3", "#1F78B4", "#1420B4"))+
  stat_compare_means(method = "wilcox.test", #method.args = list(alternative = "less"),
                     comparisons = list(c( "Young", "Old")),
                     label = "p.signif") + # Show asterisks
  labs(x = "Group",                                              # labelling x axis
       y = "Percentage",                                        # labeling y axis
       title = "Nerve",        # title
       fill = "Group") +                               # legend
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.15)) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )





##PCR region deconvolution ####
msp <- read.csv("./ROI_MSP.csv")
msp <- msp[,c(3,4,5,1,2)]
msp.gr <- makeGRangesFromDataFrame(msp, keep.extra.columns = T)

all.meth.marker.gr <- data.frame(chr=rep("chr", 27314), start=rep("chr", 27314), end=rep("chr", 27314))
all.meth.marker.gr$chr  <-  c(str_split_fixed(all.meth.marker.df$loci, "_", 2)[,1])
all.meth.marker.gr$start  <-  as.integer(str_split_fixed(all.meth.marker.df$loci, "_", 2)[,2])
all.meth.marker.gr$end  <-  as.integer(str_split_fixed(all.meth.marker.df$loci, "_", 2)[,2]) + 2
all.meth.marker.gr <- makeGRangesFromDataFrame(all.meth.marker.gr, keep.extra.columns = T)

overlap_CGI = subsetByOverlaps(all.meth.marker.gr, msp.gr)
msp.gr = msp.gr[msp.gr %over% overlap_CGI] #remain the overlapped CpG loci only 
msp.gr <- data.frame(gname = paste(seqnames(msp.gr), ranges(msp.gr), sep = "_"), 
                         beta = mcols(msp.gr)$beta, 
                         stringsAsFactors = F)
tissue.cg.df <- toi.tmp.df




