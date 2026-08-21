setwd("./DEpiBurn_EMseq/")
library(dplyr)
library(GenomicRanges)
library(nnls)
library(pheatmap)
library(RColorBrewer)

jet.colors =  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                 "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

dir <- "/datasets/work/hb-burns-meth/work/Data/level2/ONT_5mC/"
bed.files <- list.files(dir, pattern = "*collapsed.bed", full.names = T)
paired <- gsub("_cpg_collapsed.bed", "",basename(bed.files) )
paired <- gsub("_2", "",paired )

#Correlation of ONT and EMseq #####
dir <- "./data/"
bedgraph.files <- list.files(dir, pattern = "*_sd_CpG.bedGraph.gz", full.names = T)
bedgraph.files <- bedgraph.files[c(2:4,7,9,11,13,14,15,17,21,22)]

cor.list <- list()
for (i in 1:length(paired)) {
  print("Reading..")
  print(i)
  em <- read.delim(bedgraph.files[i], stringsAsFactors = F, skip = 1, header = F)
  print(bedgraph.files[i])
  rownames(em) <- paste(em$V1, em$V2, sep = "_")
  em$beta <- em$V4 /100
  
  ont <-  read.delim(bed.files[i], stringsAsFactors = F,  header = F)
  print(bed.files[i])
  #ont$V2[grep("\\+", ont$V6)] <- ont$V2[grep("\\+", ont$V6)]+1 
  #ont$minus1 <- ont$V2-1
  ont$name <- paste(ont$V1, ont$V2, sep = "_")
  ont.df <- ont %>%
    group_by(name) %>%
    summarise(
      cov = sum(V10, na.rm = TRUE),
      me = sum(V12, na.rm = TRUE), 
      unme = sum(V13, na.rm = TRUE)
    )
  ont.df$beta <- ont.df$me / ont.df$cov
  rownames(ont.df) <- ont.df$name
  keep<- intersect(rownames(em), rownames(ont.df) )
  
  cg.correlation <- cbind(em[keep, "beta"], ont.df[keep, "beta"])
  colnames(cg.correlation) <- c("EMseq", "ONT")
  
  print("Plotting...")
  ###smoothscatter---
  pdf( paste("./plots/109Cor_", i, ".pdf", sep = ""), width = 5, height = 5)
  #par(mfrow = c(1,2))
  smoothScatter(cg.correlation$EMseq, cg.correlation$ONT, 
                colramp = jet.colors, 
                nrpoints = 0, xlab = "EMseq", ylab = "ONT", 
                main = i)
  dev.off()
  
  cor.list[[i]] <- cor(cg.correlation, method = "spearman") 
}

cor.list
save(cor.list, file = "./ONTvsEMseq_correlation.rds")

#1. ONT deconvolution #####
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


dir <- "/datasets/work/hb-burns-meth/work/Data/level2/ONT_5mC/"
bed.files <- list.files(dir, pattern = "*collapsed.bed", full.names = T)

decon.df <- c()
for (i in 1:length(bed.files)) {
  print(i)
  tmp <- read.delim(bed.files[i], stringsAsFactors = F,  header = F)
  #toi.tmp <- toi.tmp[,c(1:4)]
  print(bedgraph.files[i])
  #tmp$V2[grep("\\+", tmp$V6)] <- tmp$V2[grep("\\+", tmp$V6)]+1 
  #tmp$minus1 <- tmp$V2-1
  tmp$name <- paste(tmp$V1, tmp$V2, sep = "_")
  #keep <- intersect(tmp$name , all.meth.marker)
  #tmp <- tmp[which(tmp$name %in% keep),]
  #tmp <- tmp[which(!duplicated(tmp$name)),]
  rownames(tmp) <- tmp$name
  keep<- intersect(rownames(tmp), all.meth.marker)  
  print(length(keep))
  tmp.df <- as.matrix(round(tmp[keep,"V11"]/100, 4), nrow = length(keep))
  
  ref.mat <- as.matrix(all.tissue.mat[keep,])
  mod1 <- nnls(ref.mat, tmp.df)
  decon <- mod1$x/sum(mod1$x)
  #colnames(tissue.hyper.mat)
  decon.df <- cbind(decon.df, decon)
  print("merge Done")
  #save(decon.df, file = "./deconvolution_ONT_percentage_tmp2.rds")
}
rownames(decon.df) <- colnames(all.tissue.mat)
decon.all.df <- decon.df
colnames(decon.all.df)  <- sub("\\..*", "", basename(bed.files))

decon.ont.df <- decon.all.df
save(decon.ont.df, file = "./deconvolution_ONTcollapsed_percentage.rds")
#save(decon.ont.df, file = "./deconvolution_ONT_percentage.rds")
#save(decon.ont.df, file = "./deconvolution_ONT_percentage_filter3Cov.rds")
#save(decon.ont.df, file = "./deconvolution_ONT_percentage_filter2Cov.rds")
#save(decon.ont.df, file = "./deconvolution_ONT_percentage_filter4Cov.rds")

#load(file = "./deconvolution_ONT_percentage_filter3Cov.rds")
#load(file = "./deconvolution_ONT_percentage_filter2Cov.rds")
#load(file = "./deconvolution_ONT_percentage.rds")
load(file = "./deconvolution_ONTcollapsed_percentage.rds")
load(file = "./BlisterFluid_WGBSdeconEMseq_skinblisterTOI_0.1_0.5_endo_20251219.rds")
decon.emseq.df <- decon.all.df

colnames(decon.emseq.df) <- gsub("-2","", colnames(decon.emseq.df))
colnames(decon.ont.df) <- gsub("_5mC_filter2Cov","", colnames(decon.ont.df))
colnames(decon.ont.df) <- gsub("_cpg_collapsed","", colnames(decon.ont.df))
colnames(decon.ont.df) <- gsub("_2","", colnames(decon.ont.df))

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

pdata.ont <- phenoData
colnames(decon.ont.df)[12] <- "QUT30-2"
keep <- intersect(rownames(pdata.ont), colnames(decon.ont.df) )
pdata.ont <- pdata.ont[keep,]
pdata.ont <- cbind(pdata.ont[,1:20], t(decon.ont.df))

library(ggplot2)
library(ggpubr)
ggplot(pdata.ont , aes(x=DEPTH, y=Nerve, fill = DEPTH, label = SAMPLE_ID)) +
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




keep <- intersect(colnames(decon.ont.df) , colnames(decon.emseq.df) )
decon.ont.df <- as.data.frame(t(decon.ont.df[,keep]) )
decon.emseq.df <- as.data.frame(t(decon.emseq.df[,keep]) )
#decon.emseq.df$Methods <- "EMseq"
#decon.ont.df$Methods <- "ONT"
colnames(decon.ont.df) <- paste0(colnames(decon.ont.df), "_ONT")
colnames(decon.emseq.df) <- paste0(colnames(decon.emseq.df), "_EMseq")
decon.correlation <- cbind(decon.emseq.df, decon.ont.df)

library(dplyr)
library(ggplot2)
library(ggpubr)
p1 <- decon.correlation %>%
  ggplot(  aes(x = Nerve_EMseq       , y = Nerve_ONT    )) + 
  geom_point(  size = 2.5) + 
  #scale_color_manual(values=pal.group)+
  geom_smooth(method = "lm", se = F) + 
  ggtitle("Nerve")  +stat_cor(method = "pearson", label.x = 0.05, label.y = 0.2, size = 4)  +
  xlab("EMseq Percentage") + ylab("ONT Percentage")+
  scale_y_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  scale_x_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
p2 <- decon.correlation %>%
  ggplot(  aes(x = SA_EMseq        , y = SA_ONT     )) + 
  geom_point(  size = 2.5) + 
  #scale_color_manual(values=pal.group)+
  geom_smooth(method = "lm", se = F) + 
  ggtitle("SA")  +stat_cor(method = "pearson", label.x = 0, label.y = 0.1, size = 4)  +
  xlab("EMseq Percentage") + ylab("ONT Percentage")+
  scale_y_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  scale_x_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
p3 <- decon.correlation %>%
  ggplot(  aes(x = EC_EMseq         , y = EC_ONT     )) + 
  geom_point(  size = 2.5) + 
  #scale_color_manual(values=pal.group)+
  geom_smooth(method = "lm", se = F) + 
  ggtitle("EC")  +stat_cor(method = "pearson", label.x = 0, label.y = 0.15, size = 4)  +
  xlab("EMseq Percentage") + ylab("ONT Percentage")+
  scale_y_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  scale_x_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
p4 <- decon.correlation %>%
  ggplot(  aes(x = Keratinocyte_EMseq   , y = Keratinocyte_ONT )) + 
  geom_point(  size = 2.5) + 
  #scale_color_manual(values=pal.group)+
  geom_smooth(method = "lm", se = F) + 
  ggtitle("Keratinocyte")  +stat_cor(method = "pearson", label.x = 0.7, label.y = 0.95, size = 4)  +
  xlab("EMseq Percentage") + ylab("ONT Percentage")+
  scale_y_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  scale_x_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
ggarrange(p1, p2, p3,p4, nrow = 1)  







#2. Metagenomics ####
##individual bar plot#####
dir <- "/datasets/work/hb-burns-meth/work/Data/level2/ONT_metagenomics/reports/"
braken.files <- list.files(dir, pattern = "*_G", full.names = T)

bracken.df <- NULL
for (i in 1:length(braken.files)) {
  tmp <- read.delim(braken.files[i])
  tmp$sample <- basename(braken.files[i])
  
  bracken.df <- rbind(bracken.df, tmp)
}
bracken.df$sample <- gsub("_bracken_G.txt", "", bracken.df$sample)

bracken.df.filter <- bracken.df[which(bracken.df$new_est_reads >= 50),]
table(bracken.df.filter$name, bracken.df.filter$sample)

library(ggplot2)
pal.g <- c( "black", "#8DD3C7", "#FFFFB3", "#BEBADA" ,"#FB8072", "#80B1D3" ,"#FDB462", "#B3DE69", "#FCCDE5")
names(pal.g) <- unique(bracken.df.filter$name)
bracken.df.filter$name <- factor(bracken.df.filter$name, levels = unique(bracken.df.filter$name) )
ggplot(bracken.df.filter, aes(x = sample, y = fraction_total_reads, fill = name )) +
  #theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  scale_fill_manual(values = pal.g) +
  xlab("Sample") +
  ylab("Fraction") 

##PCA plot#####
library(tidyr)
dir <- "/datasets/work/hb-burns-meth/work/Data/level2/ONT_metagenomics/reports/"
braken.files <- list.files(dir, pattern = "*_S.txt", full.names = T)
bracken.df <- NULL
for (i in 1:length(braken.files)) {
  tmp <- read.delim(braken.files[i])
  tmp$sample <- basename(braken.files[i])
  bracken.df <- rbind(bracken.df, tmp)
}
bracken.df$sample <- gsub("_bracken_S.txt", "", bracken.df$sample)
bracken.df.filter <- bracken.df[which(bracken.df$new_est_reads >= 30),]
table(bracken.df.filter$name, bracken.df.filter$sample)

bracken.mtx <- bracken.df.filter[,c(1,6,8)]
library(reshape2)
bracken.mtx <- dcast(bracken.mtx, name ~ sample, value.var = "new_est_reads", fill = 0)
rownames(bracken.mtx) <- bracken.mtx$name
bracken.mtx <- bracken.mtx[-5,-1]

sample_data <- data.frame(
  Sample_ID = colnames(bracken.mtx),
  group = c("Superficial", "Deep", "Deep", "Superficial", "Superficial", "Deep", "Superficial", "Deep", "Superficial", "Deep", "Superficial", "Full"),
  stringsAsFactors = FALSE
)
library(limma)
pal.depth <- c("Superficial"= "#999999","Deep"= "#E69F00","Full"= "#56B4E9")
plotMDS(bracken.mtx, gene.selection="common", 
        col=pal.depth[factor(sample_data$group)], dim=c(1,2), pch = 19)
legend("topright", legend=levels(factor(pdata$group)), text.col=pal.group,
       bg="white", cex=0.7)

library(RColorBrewer)
library(ggrepel)
pal.microbe <-  colorRampPalette(  brewer.pal(8,"Paired") )(12)
names(pal.microbe) <- rownames(bracken.mtx)
plotMDS(t(bracken.mtx), gene.selection="common", dim=c(1,2) , col = pal.microbe, pch = 19, cex = 3)
mds_data <- plotMDS(t(bracken.mtx))
bracken.mtx_coords <- data.frame(  Species = colnames(t(bracken.mtx)),  Dim1 = mds_data$x,  Dim2 = mds_data$y)
ggplot(bracken.mtx_coords, aes(x=Dim1, y=Dim2, color=Species)) +
  geom_point(alpha=0.6, size = 5) +         # changing the y axis nber format
  scale_color_manual(values = pal.microbe)+
  geom_text_repel(
    aes(label = Species),
    size = 5,
    fontface = "italic",
    # --- New arguments below ---
    max.overlaps = Inf,          # Forces R to label every single point regardless of overlap
    min.segment.length = 0,      # Always draw a line, even if the label is close to the point
    segment.color = "grey50",    # Color of the connecting line
    segment.size = 1,          # Thickness of the line
    box.padding = 0.8,           # Increases space around labels to help with layout
    point.padding = 0.5,         # Ensures the line doesn't touch the center of the dot
    force = 2                    # Increases the "repel" strength to push labels further out
  ) +
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )



##stringent theroldoff 
dir <- "/datasets/work/hb-burns-meth/work/Data/level2/ONT_metagenomics/reports/"
braken.files <- list.files(dir, pattern = "*_G", full.names = T)
bracken.df <- NULL
for (i in 1:length(braken.files)) {
  tmp <- read.delim(braken.files[i])
  tmp$sample <- basename(braken.files[i])
  bracken.df <- rbind(bracken.df, tmp)
}
bracken.df$sample <- gsub("_bracken_G.txt", "", bracken.df$sample)
bracken.df.filter <- bracken.df[which(bracken.df$new_est_reads >= 50),]
table(bracken.df.filter$name, bracken.df.filter$sample)
bracken.mtx <- bracken.df.filter[,c(1,6,8)]
bracken.mtx <- bracken.mtx[grep("Homo", bracken.mtx$name, invert = T),]

library(dplyr)
# 1. Calculate cohort statistics per genus
outlier_analysis <- bracken.mtx %>%
  group_by(name) %>%
  mutate(
    mean_abundance = mean(new_est_reads),
    sd_abundance = sd(new_est_reads),
    # Define threshold: 3.5 SD above the mean
    threshold = mean_abundance + (1.5 * sd_abundance)
  ) %>%
  ungroup() %>%
  # 2. Flag the outliers
  mutate(is_outlier = new_est_reads > threshold)
# 3. Filter to see only the high-confidence genera (the "11 genera")
significant_genera <- outlier_analysis %>%
  filter(is_outlier == TRUE)

# Assuming you have a column 'culture_positive' (TRUE/FALSE)
comparison_stats <- outlier_analysis %>%
  group_by(is_outlier) %>%
  summarise(
    total_count = n(),
    confirmed_by_culture = sum(culture_positive == TRUE),
    proportion = (confirmed_by_culture / total_count) * 100
  )
print(comparison_stats)

library(ggplot2)
library(dplyr)
library(scales)
# 1. Calculate stats per genus
plot_bracken <- bracken.mtx %>%
  group_by(name) %>%
  mutate(
    mean_val = mean(new_est_reads),
    sd_val = sd(new_est_reads),
    # Calculate Z-score for each observation
    z_score = (new_est_reads - mean_val) / sd_val,
    is_outlier = z_score > 3.5
  ) %>%
  ungroup()

# 2. Summary for Figure 7D (Culture confirmation)
# (Assuming you have a 'culture_match' column)
summary_stats <- plot_bracken %>%
  group_by(is_outlier) %>%
  summarise(
    prop = mean(culture_match) * 100,
    label = ifelse(is_outlier, "Outlying Genera", "Background")
  )

ggplot(plot_bracken, aes(x = new_est_reads)) +
  # Create the histogram
  geom_histogram(bins = 50, alpha = 0.7, color = "white") +
  # Add a vertical line for the threshold 
  # Note: Since different genera have different SDs, a single line 
  # is best shown on a Z-score plot. For raw reads, we use log scale:
  scale_x_log10(labels = label_number()) +
  #scale_fill_manual(values = c("gray70", "firebrick3")) +
  labs(
    title = "Distribution of Microbial Abundance",
    x = "Estimated Reads (log10)",
    y = "Frequency (Count)",
    fill = "Outlier (>3.5 SD)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )
