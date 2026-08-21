setwd("./DEpiBurn_EMseq/")
library(GenomicRanges)
library(nnls)
library(pheatmap)
library(RColorBrewer)
pal.graft <- brewer.pal(8, "Set1")[4:5]
pal.depth <- brewer.pal(4, "Greys")
pal.coll <- brewer.pal(5, "Set3")
pal.epi <- brewer.pal(4, "Reds")


#load(file = "/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/CGI_tissue_marker_matrix.RData")
#load(file = "/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/CGI_blister_marker_matrix.RData")
#load(file = "/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/CGI_blister_marker_matrix_0.1_0.5.RData")
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

pdf("./plots/101Heatmap_allCGI_signature.pdf", width = 2.5, height = 4)
breaksList = seq(0, 1, by = 0.01)
pheatmap(all.tissue.mat,   scale = "none",    clustering_method = "complete",
         show_rownames = F, show_colnames = T, fontsize_row = 12, fontsize_col = 12,
         cluster_cols = F, cluster_rows = F, annotation_legend = T,
         color = colorRampPalette(c("darkblue",  "lightgreen","yellow"))((length(breaksList))),
         breaks = breaksList, border_color = NA, 
         main = "Tissue Markers"
)
dev.off()

pheatmap(tissue.hyper.mat,   scale = "none",    clustering_method = "complete",
         show_rownames = F, show_colnames = T, fontsize_row = 12, fontsize_col = 10,
         cluster_cols = F, cluster_rows = F, annotation_legend = T,
         color = colorRampPalette(c("darkblue",  "lightgreen","yellow"))((length(breaksList))),
         breaks = breaksList, border_color = NA, 
         main = "Tissue Markers (hypermethylated)"
)
write.csv(all.tissue.mat , file = "./Skin_ref_matrix.csv")



#loading blister fluid samples
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


####EMseq Methylation pattern Deconvolution, skin blister tissue --------20251219######
dir <- "./data/"
bedgraph.files <- list.files(dir, pattern = "*_sd_CpG.bedGraph.gz", full.names = T)

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


#save(decon.all.df, file = "./BlisterFluid_WGBSdeconEMseq_skinblisterTOI.rds")
#save(decon.all.df, file = "./BlisterFluid_WGBSdeconEMseq_skinblisterTOI_0.1_0.5.rds")
save(decon.all.df, file = "./BlisterFluid_WGBSdeconEMseq_skinblisterTOI_0.1_0.5_endo.rds")

load(file = "./BlisterFluid_WGBSdeconEMseq_skinblisterTOI_0.1_0.5_endo.rds")
colnames(decon.all.df) <- substr(sub("\\..*", "", basename(bedgraph.files)), start = 24L, stop = 30L)
colnames(decon.all.df) <- gsub("_", "", colnames(decon.all.df))
phenoData <- read.csv("./SampleMeta_20220222.csv", stringsAsFactors = F)
rownames(phenoData) <- paste0("QUT", phenoData$SAMPLE_ID)
phenoData <- phenoData[colnames(decon.all.df),]
phenoData$DAYS.TO.RE.EP <- gsub("\\+", "", phenoData$DAYS.TO.RE.EP)
phenoData$DAYS.TO.RE.EP <- as.numeric(phenoData$DAYS.TO.RE.EP)

save(decon.all.df, file = "./BlisterFluid_WGBSdeconEMseq_skinblisterTOI_0.1_0.5_endo_20251219.rds")



# Cell type 
load(file = "./BlisterFluid_WGBSdeconEMseq_skinblisterTOI_0.1_0.5_endo_20251219.rds")
phenoData <- read.csv("./SampleMeta_20260108.csv", stringsAsFactors = F)
rownames(phenoData) <- paste0("QUT", phenoData$SAMPLE_ID)
phenoData <- phenoData[colnames(decon.all.df),]
phenoData$DAYS.TO.RE.EP <- gsub("\\+", "", phenoData$DAYS.TO.RE.EP)
phenoData$DAYS.TO.RE.EP <- as.numeric(phenoData$DAYS.TO.RE.EP)
data.frame(rownames(phenoData), colnames(decon.all.df))
phenoData <- cbind(phenoData, t(decon.all.df))
phenoData$tissue <- phenoData$SA + phenoData$Nerve 
idx.coll <- which(phenoData$COLL.DAY.POST.INJURY <= 1)
idx.graft <- which(phenoData$graft == "No")

phenoData$DEPTH <- gsub("superficial", "Superficial", phenoData$DEPTH)
phenoData$DEPTH <- gsub("deep", "Deep", phenoData$DEPTH)
phenoData$DEPTH <- gsub("full_thickness", "Full", phenoData$DEPTH)
phenoData$DEPTH <- factor(phenoData$DEPTH, levels = c("Superficial", "Deep", "Full"))
phenoData$reepi <- gsub("early", "Early", phenoData$reepi)
phenoData$reepi <- gsub("average", "Average", phenoData$reepi)
phenoData$reepi <- gsub("late", "Late", phenoData$reepi)
phenoData$reepi[2] <- NA
  

plot(phenoData$adipose[idx.coll] , phenoData$nerve[idx.coll] , pch = 19,
     xlim= c(0,0.3), ylim=c(0,0.3), col = pal.epi, main = "Re-Epi" )
plot( phenoData$DAYS.TO.RE.EP[idx.graft] , phenoData$adipose[idx.graft]  ,pch = 19, col = pal.epi,
      main = "Re-Epi" )
plot( phenoData$DAYS.TO.RE.EP[idx.graft] , phenoData$nerve[idx.graft]  ,pch = 19, col = pal.epi,
      main = "Re-Epi" )
plot( phenoData$DAYS.TO.RE.EP[idx.graft] , phenoData$tissue[idx.graft]  ,pch = 19, col = pal.epi,
      main = "Re-Epi" )
plot( phenoData$DAYS.TO.RE.EP[idx.graft] , phenoData$keratinocyte[idx.graft]  ,pch = 19, col = pal.epi,
      main = "Re-Epi" )
phenoData$targettissue <- phenoData$adipose + phenoData$nerve + phenoData$endo
plot( phenoData$DAYS.TO.RE.EP[idx.graft] , phenoData$targettissue[idx.graft]  ,pch = 19, #col = pal.epi,
      main = "Re-Epi" )
plot( phenoData$DAYS.TO.RE.EP[] , phenoData$targettissue[]  ,pch = 19, #col = pal.epi,
      main = "Re-Epi" )

library(ggplot2)
ggplot(data = phenoData[,], aes(x = DAYS.TO.RE.EP, y = Nerve)) + 
  geom_point(aes(colour = factor(DEPTH), size = 2)) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  geom_smooth(method = "lm", se = F) + 
  ggtitle("Nerve Conc.") +ylim(c(0,0.3)) +theme_classic()
ggplot(data = phenoData[1:19,], aes(x = DAYS.TO.RE.EP, y = SA)) + 
  geom_point(aes(colour = factor(DEPTH), size = 2)) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  geom_smooth(method = "lm", se = F) + 
  ggtitle("adipose Conc.")+ylim(c(0,0.1)) +theme_classic()
ggplot(data = phenoData[1:19,], aes(x = DAYS.TO.RE.EP, y = EC)) + 
  geom_point(aes(colour = factor(DEPTH), size = 2)) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  geom_smooth(method = "lm", se = F) + 
  ggtitle("endo Conc.")+ylim(c(0,0.1)) +theme_classic()
ggplot(data = phenoData[1:19,], aes(x = DAYS.TO.RE.EP, y = targettissue)) + 
  geom_point(aes(colour = factor(DEPTH), size = 2)) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  geom_smooth(method = "lm", se = F) + 
  ggtitle("Target Tissue Conc.")+ylim(c(0,0.3)) +theme_classic()


library(pheatmap)
library(RColorBrewer)
pal.graft <- brewer.pal(8, "Set1")[4:5]
pal.depth <- brewer.pal(3, "Greys")
#pal.coll <- brewer.pal(3, "Set3")
pal.epi <- brewer.pal(3, "Reds")
names(pal.graft) <- levels(as.factor(phenoData$graft))
names(pal.depth) <- levels(factor(phenoData$DEPTH, levels = c("Superficial", "Deep", "Full")))
names(pal.epi) <- levels(factor(phenoData$reepi))

annotation_col = data.frame(
  Graft = phenoData$graft, 
  Depth = phenoData$DEPTH, #Collection = phenoData$Coll_day,
  Reepi = phenoData$reepi
)
rownames(annotation_col) = colnames(decon.all.df)
#phenoData <- phenoData[order(phenoData$DAYS.TO.RE.EP),]
#decon.all.df <- decon.all.df[,match(rownames(phenoData), colnames(decon.all.df))]
#names(pal.coll) <- levels(as.factor(phenoData$Coll_day))
breaksList = seq(0, 1, by = 0.01)

pdf("./plots/101Heatmap_skinblisterCGI_0.15_0.3_endo.png", width = 5, height = 5)
pheatmap(decon.all.df[,order(phenoData$DEPTH)],   scale = "none",    clustering_method = "complete",
         show_rownames = T, show_colnames = T, fontsize_row = 12, fontsize_col = 10,
         cluster_cols = F, cluster_rows = F, annotation_legend = T,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList, border_color = NA,  annotation_col = annotation_col,
         #color = colorRampPalette(c("darkblue",  "lightgreen","yellow"))((length(breaksList))),
         annotation_colors  = list(
           Graft = pal.graft, Depth = pal.depth,   Reepi = pal.epi
         ),
         main = "Cell Type Deconvolution \n Blister Fluid cfDNA"
)
dev.off()



#stacked bar plot 
library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)
pal.cell <- brewer.pal(8,"Set3")
names(pal.cell) <- c("Keratinocyte", "Nerve","SA", "EC","Neutrophil","Fibroblast", "Macrophage")

decon.all_endo <- melt(phenoData[,c("Macrophage", "Nerve","Neutrophil","SA","Fibroblast","Keratinocyte", "EC", "DEPTH" , "SAMPLE_ID" )])
colnames(decon.all_endo)[3:4] <- c("Cell", "Percentage") 
decon.all_endo$Cell <- factor(decon.all_endo$Cell , levels = c("Keratinocyte", "Fibroblast","Neutrophil", "Macrophage", "EC","SA","Nerve"))
decon.all_endo$SAMPLE_ID <- paste0("Q", decon.all_endo$SAMPLE_ID)

p1 <- decon.all_endo[grep("Superficial", decon.all_endo$DEPTH),] %>%
  ggplot(data = ., mapping = aes(x = SAMPLE_ID, y = Percentage, fill = Cell)) +
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
       title = "Superficial",        # title
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
p2 <- decon.all_endo[grep("Deep", decon.all_endo$DEPTH),] %>%
  ggplot(data = ., mapping = aes(x = SAMPLE_ID, y = Percentage, fill = Cell)) +
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
       title = "Deep",        # title
       fill = "Cell") +                               # legend
  scale_y_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 90,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
p3 <- decon.all_endo[grep("Full", decon.all_endo$DEPTH),] %>%
  ggplot(data = ., mapping = aes(x = SAMPLE_ID, y = Percentage, fill = Cell)) +
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
       title = "Full Thickness",        # title
       fill = "Cell") +                               # legend
  scale_y_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 90,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
ggarrange(p1, p2, p3, nrow = 1, common.legend = T)



# individual cell type
phenoData$SAMPLE_ID <- paste0("QUT", phenoData$SAMPLE_ID)
library(ggplot2)
p1 <- ggplot(phenoData , aes(x=DEPTH, y=Nerve, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  #geom_text(aes(colour = DEPTH), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #scale_colour_manual(values=c("#444444", "#A77305", "#0690A1"))+
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less"),
                     comparisons = list(c( "Superficial", "Deep"), c("Deep", "Full")),
                     label = "p.signif") + # Show asterisks
  labs(x = "Group",                                              # labelling x axis
       y = "Percentage",                                        # labeling y axis
       title = "Nerve",        # title
       fill = "Group") +                               # legend
  scale_y_continuous(labels = scales::percent_format(), limits = c(0.05, 0.25)) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
p2 <- ggplot(phenoData , aes(x=DEPTH, y=SA, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  #geom_text(aes(colour = DEPTH), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #scale_colour_manual(values=c("#444444", "#A77305", "#0690A1"))+
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less"),
                     comparisons = list(c( "Superficial", "Deep"), c("Deep", "Full")),
                     label = "p.signif") + # Show asterisks
  labs(x = "Group",                                              # labelling x axis
       y = "Percentage",                                        # labeling y axis
       title = "SA",        # title
       fill = "Group") +                               # legend
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.05)) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
p3 <- ggplot(phenoData , aes(x=DEPTH, y=EC, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  #geom_text(aes(colour = DEPTH), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #scale_colour_manual(values=c("#444444", "#A77305", "#0690A1"))+
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less"),
                     comparisons = list(c( "Superficial", "Deep"), c("Deep", "Full")),
                     label = "p.signif") + # Show asterisks
  labs(x = "Group",                                              # labelling x axis
       y = "Percentage",                                        # labeling y axis
       title = "EC",        # title
       fill = "Group") +                               # legend
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.1)) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
p4 <- ggplot(phenoData , aes(x=DEPTH, y=Keratinocyte, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  #geom_text(aes(colour = DEPTH), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #scale_colour_manual(values=c("#444444", "#A77305", "#0690A1"))+
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less"),
                     comparisons = list(c( "Superficial", "Deep"), c("Deep", "Full")),
                     label = "p.signif") + # Show asterisks
  labs(x = "Group",                                              # labelling x axis
       y = "Percentage",                                        # labeling y axis
       title = "Keratinocyte",        # title
       fill = "Group") +                               # legend
  scale_y_continuous(labels = scales::percent_format(), limits = c(0.5, 1)) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
ggarrange(p1, p2, p3,  nrow = 1, common.legend = T)



#Re epitheliation
phenoData$reepi <- factor(phenoData$reepi, levels = c("Early", "Average", "Late"))
names(pal.epi) <- levels(factor(phenoData$reepi))
pdata.df <- phenoData[complete.cases(phenoData$reepi),]
p1 <- pdata.df[grep("No", pdata.df$graft),] %>% ggplot( aes(x=reepi , y=Nerve, fill = reepi , label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  #geom_text(aes(colour = DEPTH), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=pal.epi)+
  #scale_colour_manual(values=c("#444444", "#A77305", "#0690A1"))+
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less"),
                     comparisons = list(c("Early", "Average"), c("Average", "Late")), 
                     label = "p.signif") + # Show asterisks
  labs(x = "Group",                                              # labelling x axis
       y = "Percentage",                                        # labeling y axis
       title = "Nerve",        # title
       fill = "Days to Re-epitheliation") +                               # legend
  scale_y_continuous(labels = scales::percent_format(), limits = c(0.075, 0.25)) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
p2 <- pdata.df[grep("No", pdata.df$graft),] %>% ggplot( aes(x=reepi , y=SA, fill = reepi , label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  #geom_text(aes(colour = DEPTH), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=pal.epi)+
  #scale_colour_manual(values=c("#444444", "#A77305", "#0690A1"))+
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less"),
                     comparisons = list(c("Early", "Average"), c("Average", "Late")), 
                     label = "p.signif") + # Show asterisks
  labs(x = "Group",                                              # labelling x axis
       y = "Percentage",                                        # labeling y axis
       title = "SA",        # title
       fill = "Group") +                               # legend
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.03)) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
p3 <- pdata.df[grep("No", pdata.df$graft),]  %>% ggplot( aes(x=reepi , y=SA, fill = reepi , label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  #geom_text(aes(colour = DEPTH), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=pal.epi)+
  #scale_colour_manual(values=c("#444444", "#A77305", "#0690A1"))+
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less"),
                     comparisons = list(c("Early", "Average"), c("Average", "Late")), 
                     label = "p.signif") + # Show asterisks
  labs(x = "Group",                                              # labelling x axis
       y = "Percentage",                                        # labeling y axis
       title = "EC",        # title
       fill = "Group") +                               # legend
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.03)) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
ggarrange(p1, p2, p3,  nrow = 1, common.legend = T)




save(phenoData, file = "./phenoData_EMseq_deconvolution_percentage.rds")




####20240326 barplots of individual samples#####
write.csv(decon.all.df, file = "./tmp.csv")

library(ggplot2)
library(dplyr)
pal.cell <- brewer.pal(8,"Set3")
decon.all_endo <- read.csv("./decon_all_endo_df.csv")
decon.all_endo$Cell <- factor(decon.all_endo$Cell, levels = c("keratinocyte", "fibroblast", "macrophage", "neutrophil", "SA", "Endo_vein", "nerve"))
decon.all_endo[grep("superficial", decon.all_endo$Depth),] %>%
  ggplot(data = ., mapping = aes(x = ID, y = Percentage, fill = Cell)) +
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
       title = "Superficial",        # title
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
decon.all_endo[grep("deep", decon.all_endo$Depth),] %>%
  ggplot(data = ., mapping = aes(x = ID, y = Percentage, fill = Cell)) +
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
       title = "Deep",        # title
       fill = "Cell") +                               # legend
  scale_y_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 90,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
decon.all_endo[grep("full", decon.all_endo$Depth),] %>%
  ggplot(data = ., mapping = aes(x = ID, y = Percentage, fill = Cell)) +
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
       title = "Full Thickness",        # title
       fill = "Cell") +                               # legend
  scale_y_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 90,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )


##SA
all.tissue.mat[ grep("SA", all.meth.marker.df$tissue),]
tmp <- all.tissue.mat[ grep("SA", all.meth.marker.df$tissue),]
tmp[201:300,]
all.tissue.mat[ grep("fibroblast", all.meth.marker.df$tissue),]
tmp <- all.tissue.mat[ grep("fibroblast", all.meth.marker.df$tissue),]
tmp <- tmp[which(!tmp$fibroblast == 1 ),]
tmp <- tmp[which(tmp$fibroblast >0.5 ),]

##fibroblast
all.tissue.mat[ grep("fibroblast", all.meth.marker.df$tissue),]
tmp <- all.tissue.mat[ grep("fibroblast", all.meth.marker.df$tissue),]
fibroblast <- all.meth.marker.df[grep("fibroblast", all.meth.marker.df$tissue),]
fibroblast <- fibroblast[grep("hyper", fibroblast$marker),]
library(stringr)
fibroblast <- as.data.frame(str_split_fixed(fibroblast$loci, "_", 2))
fibroblast$V2 <- as.integer(fibroblast$V2)
colnames(fibroblast) <- c("seqname", "start")
fibroblast$end <- fibroblast$start+2

fibroblast <- makeGRangesFromDataFrame(fibroblast)
fibroblast <- reduce(fibroblast, min.gapwidth = 500L)
fibroblast <- fibroblast[width(fibroblast) > 500,]
ranges(fibroblast)
write.csv(fibroblast, file = "./fibroblast_decon_500gap.csv")

##Vascular Endo
all.tissue.mat[ grep("Endo_vein", all.meth.marker.df$tissue),]
tmp <- all.tissue.mat[ grep("Endo_vein", all.meth.marker.df$tissue),]
Endo_vein <- all.meth.marker.df[grep("Endo_vein", all.meth.marker.df$tissue),]
Endo_vein <- Endo_vein[grep("hyper", Endo_vein$marker),]
tmp <- tmp[which(tmp$Endo_vein >0.5 ),]




####EMseq Methylation pattern Deconvolution######
dir <- "./data/"
bedgraph.files <- list.files(dir, pattern = "*_sd_CpG.bedGraph.gz", full.names = T)

all.meth.marker <- c(tissue.marker.hyper$loci,tissue.marker.un$loci) 
all.tissue.mat <- rbind(tissue.hyper.mat, tissue.un.mat)

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
save(decon.all.df, file = "./BlisterFluid_WGBSdeconEMseq.rds")


load(file = "./BlisterFluid_WGBSdeconEMseq.rds")
plot(decon.all.df["adipose",] , decon.all.df["nerve",], pch = 19,
     xlim= c(0,0.1), ylim=c(0,0.1), col = pal.graft, main = "graft")
plot(decon.all.df["adipose",] , decon.all.df["nerve",], pch = 19,
     xlim= c(0,0.05), ylim=c(0,0.05), col = pal.depth, main = "depth")
plot(decon.all.df["adipose",] , decon.all.df["nerve",], pch = 19,
     xlim= c(0,0.1), ylim=c(0,0.1), col = pal.epi, main = "Re-Epi")

phenoData$adipose <- decon.all.df["adipose",]
phenoData$nerve <- decon.all.df["nerve",]
phenoData$bladder <- decon.all.df["bladder",]
phenoData$skin <- decon.all.df["skin",]
phenoData$tissue <- phenoData$adipose + phenoData$nerve
idx.coll <- which(phenoData$COLL.DAY.POST.INJURY <= 1)
plot(phenoData$adipose[idx.coll] , phenoData$nerve[idx.coll] , pch = 19,
     xlim= c(0,0.1), ylim=c(0,0.1), col = pal.epi, main = "Re-Epi" )
plot(phenoData$adipose[idx.coll]  , phenoData$bladder[idx.coll] , pch = 19,
     xlim= c(0,0.1), ylim=c(0,0.1), col = pal.epi, main = "Re-Epi" )
plot(phenoData$adipose[idx.coll]  , phenoData$DAYS.TO.RE.EP[idx.coll] , pch = 19,
     xlim= c(0,0.1), ylim=c(0,0.1), main = "Re-Epi" )

phenoData_nongraft <- phenoData[phenoData$graft == "No",]

plot(phenoData_nongraft$DAYS.TO.RE.EP ~ phenoData_nongraft$adipose)
plot(phenoData_nongraft$DAYS.TO.RE.EP ~ phenoData_nongraft$nerve)
plot(phenoData_nongraft$DAYS.TO.RE.EP ~ phenoData_nongraft$bladder)
plot(phenoData_nongraft$DAYS.TO.RE.EP ~ phenoData_nongraft$skin)
library("ggpubr")
ggqqplot(phenoData_nongraft$DAYS.TO.RE.EP, ylab = "DAYS.TO.RE.EP")
ggscatter(phenoData_nongraft[phenoData_nongraft$DAYS.TO.RE.EP<40,], x = "DAYS.TO.RE.EP", y = "adipose", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "DAYS.TO.RE.EP", ylab = "adipose")
ggscatter(phenoData_nongraft[phenoData_nongraft$DAYS.TO.RE.EP<40,], x = "DAYS.TO.RE.EP", y = "nerve", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "DAYS.TO.RE.EP", ylab = "nerve")
ggscatter(phenoData_nongraft[phenoData_nongraft$DAYS.TO.RE.EP<40,], x = "DAYS.TO.RE.EP", y = "bladder", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "DAYS.TO.RE.EP", ylab = "bladder")
ggscatter(phenoData_nongraft[phenoData_nongraft$DAYS.TO.RE.EP<40,], x = "DAYS.TO.RE.EP", y = "tissue", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "DAYS.TO.RE.EP", ylab = "tissue")

##reduced tissue type
dir <- "./data/"
bedgraph.files <- list.files(dir, pattern = "*_sd_CpG.bedGraph.gz", full.names = T)

all.meth.marker <- c(tissue.marker.hyper$loci,tissue.marker.un$loci) 
all.tissue.mat <- rbind(tissue.hyper.mat, tissue.un.mat)
reduce.meth.marker <- reduce.meth.marker.df$loci
reduce.meth.mat <- all.tissue.mat[reduce.meth.marker.df$loci, c("adipose","macrophage","nerve","neutrophil","skin")]

decon.df <- c()
for (i in 1:length(bedgraph.files)) {
  print(i)
  tmp <- read.delim(bedgraph.files[i], stringsAsFactors = F, skip = 1, header = F)
  #toi.tmp <- toi.tmp[,c(1:4)]
  print(bedgraph.files[i])
  rownames(tmp) <- paste(tmp$V1, tmp$V2, sep = "_")
  keep<- intersect(rownames(tmp), reduce.meth.marker)
  tmp.df <- as.matrix(round(tmp[keep,"V4"]/100, 4), nrow = length(keep))
  
  ref.mat <- as.matrix(reduce.meth.mat[keep,])
  mod1 <- nnls(ref.mat, tmp.df)
  decon <- mod1$x/sum(mod1$x)
  #colnames(tissue.hyper.mat)
  decon.df <- cbind(decon.df, decon)
  print("merge Done")
}
rownames(decon.df) <- colnames(reduce.meth.mat)
decon.all.df <- decon.df
colnames(decon.all.df)  <- sub("\\..*", "", basename(bedgraph.files))
save(decon.all.df, file = "./BlisterFluid_reduceWGBSdeconEMseq.rds")

phenoData <- cbind(phenoData, t(decon.all.df))
phenoData$DEPTH <- factor(phenoData$DEPTH , levels=c("superficial", "deep", "full_thickness"))
pdf("./plots/101Plot_nerve_adipose_allCGI_reducesample.pdf", width = 5, height = 5)
boxplot(phenoData$nerve ~ phenoData$DEPTH)
plot(nerve ~ AGE_months, data = phenoData[idx, ])
plot(nerve ~ DAYS.TO.RE.EP, data = phenoData[idx, ])
boxplot(phenoData$adipose ~ phenoData$DEPTH)
plot(adipose ~ AGE_months, data = phenoData[idx, ])
plot(adipose ~ DAYS.TO.RE.EP, data = phenoData[idx, ])
plot(adipose ~ nerve, data = phenoData[idx, ])
dev.off()




####Plot ##########
load("./BlisterFluid_WGBSdeconEMseq.rds")
load("./BlisterFluid_reduceWGBSdeconEMseq.rds")
colnames(decon.all.df) <- substr(colnames(decon.all.df), start = 24L, stop = 30L)
colnames(decon.all.df) <- gsub("_", "", colnames(decon.all.df))
phenoData <- read.csv("./SampleMeta_20220222.csv", stringsAsFactors = F)
rownames(phenoData) <- paste0("QUT", phenoData$SAMPLE_ID)
phenoData <- phenoData[colnames(decon.all.df),]
phenoData$DAYS.TO.RE.EP <- gsub("\\+", "", phenoData$DAYS.TO.RE.EP)
phenoData$DAYS.TO.RE.EP <- as.numeric(phenoData$DAYS.TO.RE.EP)

library(pheatmap)
library(RColorBrewer)
pal.graft <- brewer.pal(8, "Set1")[4:5]
pal.depth <- brewer.pal(4, "Greys")
pal.coll <- brewer.pal(4, "Set3")
pal.epi <- brewer.pal(4, "Reds")
annotation_col = data.frame(
  Graft = phenoData$graft, 
  Depth = phenoData$DEPTH, Collection = phenoData$Coll_day,
  Reepi = phenoData$reepi
)
rownames(annotation_col) = colnames(decon.all.df)

phenoData <- phenoData[order(phenoData$DAYS.TO.RE.EP),]
decon.all.df <- decon.all.df[,match(rownames(phenoData), colnames(decon.all.df))]

names(pal.graft) <- levels(as.factor(phenoData$graft))
names(pal.depth) <- levels(factor(phenoData$DEPTH, levels = c("superficial", "deep", "full_thickness", NA)))
names(pal.epi) <- levels(factor(phenoData$reepi))
names(pal.coll) <- levels(as.factor(phenoData$Coll_day))
breaksList = seq(0, 1, by = 0.01)
idx <- which(phenoData$Coll_day == "d0-3" & phenoData$graft == "No")
idx <- which(phenoData$reepi == "late" )

pdf("./plots/101Heatmap_allCGI_reducesample.pdf", width = 6, height = 4.5)
pheatmap(decon.all.df[,],   scale = "none",    clustering_method = "complete",
         show_rownames = T, show_colnames = T, fontsize_row = 10, fontsize_col = 5,
         cluster_cols = F, cluster_rows = F, annotation_legend = T,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList, border_color = NA,  annotation_col = annotation_col,
         annotation_colors  = list(
           Graft = pal.graft, Depth = pal.depth,  Collection = pal.coll, Reepi = pal.epi
         ),
         main = "All"
)
pheatmap(decon.all.df[,idx],   scale = "none",    clustering_method = "complete",
         show_rownames = T, show_colnames = T, fontsize_row = 10, fontsize_col = 5,
         cluster_cols = F, cluster_rows = F, annotation_legend = T,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList, border_color = NA,  annotation_col = annotation_col,
         annotation_colors  = list(
           Graft = pal.graft, Depth = pal.depth,  Collection = pal.coll, Reepi = pal.epi
         ),
         main = "Selective"
) #bladder marker
dev.off()







####Bladder marker in blister fluid######
bladder.marker <- all.meth.marker.df[grep("bladder", all.meth.marker.df$tissue),]
#chr2:176068538-176123675 EVX2 HOXDfamily
#chr7:27225121-27246295 EVX1
dir <- "./data/"
bedgraph.files <- list.files(dir, pattern = "*_sd_CpG.bedGraph.gz", full.names = T)

beta.bladder <- c()
cov.bladder <- c()
tmp <- read.delim(bedgraph.files[1], stringsAsFactors = F, skip = 1, header = F)
tmp$V7 <- tmp$V5+tmp$V6
print(bedgraph.files[1])
rownames(tmp) <- paste(tmp$V1, tmp$V2, sep = "_")
keep<- intersect(rownames(tmp), bladder.marker$loci)
tmp.beta <- tmp[keep,4]
tmp.cov <- tmp[keep,7]
beta.bladder <- cbind(beta.bladder,tmp.beta)
cov.bladder <- cbind(cov.bladder, tmp.cov)
rownames(beta.bladder) <- keep
rownames(cov.bladder) <- keep
for (i in 2:length(bedgraph.files)) {
  print(i)
  tmp <- read.delim(bedgraph.files[i], stringsAsFactors = F, skip = 1, header = F)
  tmp$V7 <- tmp$V5+tmp$V6
  print(bedgraph.files[i])
  rownames(tmp) <- paste(tmp$V1, tmp$V2, sep = "_")
  keep<- intersect(rownames(tmp), rownames(beta.bladder))
  tmp.beta <- tmp[keep,4]
  tmp.cov <- tmp[keep,7]
  beta.bladder <- beta.bladder[keep,]
  cov.bladder <- cov.bladder[keep,]
  beta.bladder <- cbind(beta.bladder,tmp.beta)
  cov.bladder <- cbind(cov.bladder, tmp.cov)
  print(nrow(beta.bladder))
}

phenoData <- read.csv("./SampleMeta_20220222.csv", stringsAsFactors = F)
rownames(phenoData) <- paste0("QUT", phenoData$SAMPLE_ID)
colnames(beta.bladder) <-  substr(bedgraph.files, start = 32L, stop = 38L)
colnames(beta.bladder) <-  gsub("_", "", colnames(beta.bladder))
colnames(cov.bladder) <-  colnames(beta.bladder) 
phenoData <- phenoData[grep("181", rownames(phenoData), invert = T),]
keep <- intersect(rownames(phenoData), colnames(beta.bladder))
phenoData <- phenoData[keep,]
beta.bladder <- beta.bladder[,keep]
library(pheatmap)
library(RColorBrewer)
pal.graft <- brewer.pal(8, "Set1")[4:5]
pal.depth <- brewer.pal(4, "Greys")
pal.coll <- brewer.pal(4, "Set3")
pal.epi <- brewer.pal(4, "Reds")
annotation_col = data.frame(
  Graft = phenoData$graft, 
  Depth = phenoData$DEPTH, Collection = phenoData$Coll_day,
  Reepi = phenoData$Re.epi
)
rownames(annotation_col) = colnames(beta.bladder)

phenoData <- phenoData[order(phenoData$graft),]
beta.bladder <- beta.bladder[,match(rownames(phenoData), colnames(beta.bladder))]

names(pal.graft) <- levels(as.factor(phenoData$graft))
names(pal.depth) <- levels(factor(phenoData$DEPTH, levels = c("superficial", "deep", "full_thickness")))
names(pal.coll) <- levels(as.factor(phenoData$Coll_day))
names(pal.epi) <- levels(factor(phenoData$Re.epi, levels = c("d5to7", "d7to21", "above21", "")))
breaksList = seq(0, 100, by = 0.01)
order(phenoData$graft)
pdf("./plots/101Heatmap_bladder.marker.pdf", width = 6, height = 4.5)
pheatmap(beta.bladder[,c(11:13,15:20)],   scale = "none",    clustering_method = "complete",
         show_rownames = T, show_colnames = T, fontsize_row = 10, fontsize_col = 5,
         cluster_cols = F, cluster_rows = F, annotation_legend = T,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList, border_color = NA,  annotation_col = annotation_col,
         annotation_colors  = list(
           Graft = pal.graft, Depth = pal.depth,  Collection = pal.coll, Reepi = pal.epi
         ),
         main = "All"
)
dev.off()
#



#####Bladder decon model #####
tmp <- read.delim(bedgraph.files[1], stringsAsFactors = F, skip = 1, header = F)
#toi.tmp <- toi.tmp[,c(1:4)]
print(bedgraph.files[1])
rownames(tmp) <- paste(tmp$V1, tmp$V2, sep = "_")
keep<- intersect(rownames(tmp), all.meth.marker)
tmp.df <- as.matrix(round(tmp[keep,"V4"]/100, 4), nrow = length(keep))

ref.mat <- as.matrix(all.tissue.mat[keep,])
mod1 <- nnls(ref.mat, tmp.df)





#####All CpG Tissue model ######
load(file = "./tissue_marker_mat.RData")
load(file = "./CellLine_marker_mat.RData")
#tissue.marker <- tissue.marker[grep("ADSC|ciPSC|endo_proge", tissue.marker$tissue , invert = T),]
rownames(tissue.marker) <- tissue.marker$loci
##all meth marker
#tissue.all.mat <- tissue.all.mat[, grep("ADSC|ciPSC|endo_proge", colnames(tissue.all.mat), invert = T)]


dir <- "./data/"
bedgraph.files <- list.files(dir, pattern = "*_sd_CpG.bedGraph.gz", full.names = T)
id.name <- gsub("_", "",substr(bedgraph.files, start = 32L, stop = 38L) )
phenoData <- read.csv("./SampleMeta_20220222.csv", stringsAsFactors = F)
rownames(phenoData) <- paste0("QUT", phenoData$SAMPLE_ID)
phenoData <- phenoData[id.name,]
phenoData$DAYS.TO.RE.EP <- gsub("\\+", "", phenoData$DAYS.TO.RE.EP)
phenoData$DAYS.TO.RE.EP <- as.numeric(phenoData$DAYS.TO.RE.EP)

all.meth.marker <- tissue.marker$loci
all.tissue.mat <- tissue.all.mat
all.meth.marker <- tissue.marker.hyper$loci
all.tissue.mat <- tissue.hyper.mat


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
colnames(decon.all.df)  <- id.name


annotation_col = data.frame(
  Graft = phenoData$graft, 
  Depth = phenoData$DEPTH, Collection = phenoData$Coll_day,
  Reepi = phenoData$Re.epi
)
rownames(annotation_col) = colnames(decon.all.df)

phenoData <- phenoData[order(phenoData$DAYS.TO.RE.EP),]
decon.all.df <- decon.all.df[,match(rownames(phenoData), colnames(decon.all.df))]
idx <- which(phenoData$Coll_day == "d0-3" & phenoData$graft == "No")
idx <- which(phenoData$DEPTH == "deep" )

names(pal.graft) <- levels(as.factor(phenoData$graft))
names(pal.depth) <- levels(factor(phenoData$DEPTH, levels = c("superficial", "deep", "full_thickness")))
names(pal.coll) <- levels(as.factor(phenoData$Coll_day))
names(pal.epi) <- levels(factor(phenoData$Re.epi, levels = c("d5to7", "d7to21", "above21", "")))
breaksList = seq(0, 1, by = 0.01)
pdf("./plots/101Heatmap_Blister_hyperCpG_Tissue.pdf", width = 6, height = 4.5)
pheatmap(decon.all.df[,],   scale = "none",    clustering_method = "complete",
         show_rownames = T, show_colnames = T, fontsize_row = 15, fontsize_col = 15,
         cluster_cols = F, cluster_rows = F, annotation_legend = T,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList, border_color = NA,  annotation_col = annotation_col,
         annotation_colors  = list(
           Graft = pal.graft, Depth = pal.depth,  Collection = pal.coll, Reepi = pal.epi
         ),
         main = "All"
)
pheatmap(decon.all.df[,idx],   scale = "none",    clustering_method = "complete",
         show_rownames = T, show_colnames = T, fontsize_row = 15, fontsize_col = 18,
         cluster_cols = F, cluster_rows = F, annotation_legend = T,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList, border_color = NA,  annotation_col = annotation_col,
         annotation_colors  = list(
           Graft = pal.graft, Depth = pal.depth,  Collection = pal.coll, Reepi = pal.epi
         ),
         main = "All"
)
dev.off()

save(decon.all.df, phenoData, file = "./Decon_hyperCpG_Tissue.RData")


load("./Decon_allCpG_Tissue.RData")
load("./Decon_hyperCpG_Tissue.RData")
decon.all.df[5,idx] /decon.all.df[1,idx] 
(decon.all.df[1,idx] +decon.all.df[5,idx] )/decon.all.df[7,idx]
barplot(decon.all.df["nerve",idx]) ###Nerve correlated with depth and re-epi time ???
tissue.marker.hyper[grep("nerve", tissue.marker.hyper$tissue),]
phenoData <- cbind(phenoData, t(decon.all.df))
phenoData$DEPTH <- factor(phenoData$DEPTH , levels=c("superficial", "deep", "full_thickness"))
pdf("./plots/101Plot_nerve_adipose_allCpG.pdf", width = 5, height = 5)
boxplot(phenoData$nerve ~ phenoData$DEPTH)
plot(nerve ~ AGE_months, data = phenoData[idx, ])
plot(nerve ~ DAYS.TO.RE.EP, data = phenoData[idx, ])
boxplot(phenoData$SA ~ phenoData$DEPTH)
plot(SA ~ AGE_months, data = phenoData[idx, ])
plot(SA ~ DAYS.TO.RE.EP, data = phenoData[idx, ])
plot(SA ~ nerve, data = phenoData[idx, ])
dev.off()

boxplot(dermis ~ graft, data = phenoData[idx,])
plot(dermis ~ AGE_months, data = phenoData[idx,])



hyper.loci <- tissue.marker.hyper[grep("nerve", tissue.marker.hyper$tissue),]
rownames(hyper.loci) <- hyper.loci$loci
tissue.hyper.mat
load(file = "./tissue_marker_mat.RData")
tissue.marker.hyper[grep("nerve", tissue.marker.hyper$tissue),]
keep <- intersect(rownames(hyper.loci),tissue.marker.hyper[grep("nerve", tissue.marker.hyper$tissue), "loci"] )
hyper.loci <- hyper.loci[keep,]


##Nerve ROI 







# remove ADSC, endo_progenitor
# Embryonic Fibroblast stands out
library(stringr)
ROI <- tissue.marker[tissue.marker$tissue =="EmbryonicFibroblast",] #Embryonic Fibroblast linked with depth and days to re-epitheliation 
ROI$seqname <- str_split_fixed(ROI$loci, "_", 2)[,1]
ROI$start <- as.integer( str_split_fixed(ROI$loci, "_", 2)[,2])
ROI$end <- ROI$start + 2
ROI <- makeGRangesFromDataFrame(ROI)
ROI <- reduce(ROI, min.gapwidth = 30L)
ROI.df <- data.frame(seq = seqnames(ROI), start = start(ROI), end = end(ROI), stringsAsFactors = F)
ROI.df <- ROI.df[which(ROI.df$end - ROI.df$start > 70),]
ROI.df$length <- ROI.df$end - ROI.df$start
write.csv(ROI.df, file = "tmp.csv")
