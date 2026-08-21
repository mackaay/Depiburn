setwd("./DEpiBurn_EMseq/")
library(ggplot2)
library(ggpubr)
library(grid)

#http://sciprim.com/html/copyNumb.v2.0.html
#human genome : 6217275501
#use 20 ul of 25 for ddPCR assay as the input 

ddpcr <- read.csv("./ddPCR_LOD.csv")
ddpcr <- read.csv("./ddPCR_LOD_20230602.csv")


###Placenta--MON2####
idx <- grep("MON2", ddpcr$target)
lod <- ddpcr[idx,]
idx <- c(1:4,7:10, 13:16,19:22,25:28,31:34)
# Add the regression line
ggplot(lod[idx, ], aes(x=log10(spike), y=log10(concentration*20), color=gene)) + 
  geom_point() + 
  geom_smooth(method=lm)+scale_color_brewer(palette="Dark2") +
  labs(title="MON2",
       x="SpikeIn (Copies) LOG", y = "ddPCR readout (Copies) LOG")+
  theme_classic() + xlim(0.5,4) + ylim(0.5,4)
dev.off()
cor.test(log10(lod[c(1:4,13:16,25:28),]$concentration*20), log10(lod[c(1:4,13:16,25:28),]$spike), method = "pearson")
cor.test(log10(lod[c(7:10,19:22,31:34),]$concentration*20), log10(lod[c(7:10,19:22,31:34),]$spike), method = "pearson")
# Create a text
grob1 <- grobTree(textGrob("R=0.992 \np = 2.385e-10", x=0.1,  y=0.75, hjust=0,
                          gp=gpar(col="darkgreen", fontsize=10, fontface="italic")))
grob2 <- grobTree(textGrob("R=0.998 \np = 7.98e-14", x=0.65,  y=0.1, hjust=0,
                           gp=gpar(col="darkred", fontsize=10, fontface="italic")))
png("./plots/105LOD_ddPCR_MON2.png", width = 10, height = 8, res = 300, units = "cm")
ggplot(lod[idx, ], aes(x=log10(spike), y=log10(concentration*20), color=gene)) + 
  geom_point() + 
  geom_smooth(method=lm)+scale_color_brewer(palette="Dark2") +
  labs(title="MON2",
       x="SpikeIn (Copies) LOG", y = "ddPCR readout (Copies) LOG")+
  theme_classic() + xlim(0.5,4) + ylim(0.5,4)+ 
  annotation_custom(grob1) + annotation_custom(grob2)
dev.off()



idx <- grep("MON2", ddpcr$target)
cv <- ddpcr[idx,]
cv <- cv[cv$gene == "MON2",]
sd(cv$concentration[c(1,7,13)])
sd(cv$concentration[c(1,7,13)])/mean(cv$concentration[c(1,7,13)])*100
sd(cv$concentration[c(2,8,14)])
sd(cv$concentration[c(2,8,14)])/mean(cv$concentration[c(2,8,14)])*100
sd(cv$concentration[c(3,9,15)])
sd(cv$concentration[c(3,9,15)])/mean(cv$concentration[c(3,9,15)])*100
sd(cv$concentration[c(4,10,16)])
sd(cv$concentration[c(4,10,16)])/mean(cv$concentration[c(4,10,16)])*100
sd(cv$concentration[c(5,11,17)])
sd(cv$concentration[c(5,11,17)])/mean(cv$concentration[c(5,11,17)])*100


###Placenta--FLT1####
idx <- grep("FLT1", ddpcr$target)
lod <- ddpcr[idx,]
idx <- c(1:4,7:10, 13:16,19:22,25:28,31:34)
lod$concentration[4]  <- 0.1
# Add the regression line
ggplot(lod[idx, ], aes(x=log10(spike), y=log10(concentration*20), color=gene)) + 
  geom_point() + 
  geom_smooth(method=lm)+scale_color_brewer(palette="Dark2") +
  labs(title="FLT1",
       x="SpikeIn (Copies) LOG", y = "ddPCR readout (Copies) LOG")+
  theme_classic() + xlim(0.5,4) + ylim(0.5,4)
dev.off()
cor.test(log10(lod[c(1:4,13:16,25:28),]$concentration*20), log10(lod[c(1:4,13:16,25:28),]$spike), method = "pearson")
cor.test(log10(lod[c(7:10,19:22,31:34),]$concentration*20), log10(lod[c(7:10,19:22,31:34),]$spike), method = "pearson")
# Create a text
grob1 <- grobTree(textGrob("R=0.984 \np = 8.269e-09", x=0.1,  y=0.75, hjust=0,
                           gp=gpar(col="darkgreen", fontsize=10, fontface="italic")))
grob2 <- grobTree(textGrob("R=0.990 \np = 6.627e-10", x=0.65,  y=0.1, hjust=0,
                           gp=gpar(col="darkred", fontsize=10, fontface="italic")))

png("./plots/105LOD_ddPCR_FLT1.png", width = 10, height = 8, res = 300, units = "cm")
ggplot(lod[idx, ], aes(x=log10(spike), y=log10(concentration*20), color=gene)) + 
  geom_point() + 
  geom_smooth(method=lm)+scale_color_brewer(palette="Dark2") +
  labs(title="FLT1",
       x="SpikeIn (Copies) LOG", y = "ddPCR readout (Copies) LOG")+
  theme_classic() + xlim(0,4) + ylim(0,4)+ 
  annotation_custom(grob1) + annotation_custom(grob2)
dev.off()



idx <- grep("FLT1", ddpcr$target)
cv <- ddpcr[idx,]
cv <- cv[cv$gene == "FLT1",]
sd(cv$concentration[c(1,7,13)])
sd(cv$concentration[c(1,7,13)])/mean(cv$concentration[c(1,7,13)])*100
sd(cv$concentration[c(2,8,14)])
sd(cv$concentration[c(2,8,14)])/mean(cv$concentration[c(2,8,14)])*100
sd(cv$concentration[c(3,9,15)])
sd(cv$concentration[c(3,9,15)])/mean(cv$concentration[c(3,9,15)])*100
sd(cv$concentration[c(4,10,16)])
sd(cv$concentration[c(4,10,16)])/mean(cv$concentration[c(4,10,16)])*100
sd(cv$concentration[c(5,11,17)])
sd(cv$concentration[c(5,11,17)])/mean(cv$concentration[c(5,11,17)])*100







ddpcr <- read.csv("./ddPCR_LOD_nerve.csv")
###Nerve--LAMP5AS1####
idx <- grep("LAMP5AS1_F1", ddpcr$target)
lod <- ddpcr[idx,]
idx <- c(1:4,7:10, 13:16,19:22,25:28,31:34)
# Add the regression line
ggplot(lod[idx, ], aes(x=log10(spike), y=log10(concentration*20), color=gene)) + 
  geom_point() + 
  geom_smooth(method=lm)+scale_color_brewer(palette="Dark2") +
  labs(title="Nerve",
       x="SpikeIn (Copies) \n Log10-scaled", y = "ddPCR readout (Copies) \n Log10-scaled")+
  theme_classic() + xlim(0.5,4) + ylim(0.5,4)
dev.off()
cor.test(log10(lod[c(1:4,13:16,25:28),]$concentration*20), log10(lod[c(1:4,13:16,25:28),]$spike), method = "pearson")
cor.test(log10(lod[c(7:10,19:22,31:34),]$concentration*20), log10(lod[c(7:10,19:22,31:34),]$spike), method = "pearson")
# Create a text
grob1 <- grobTree(textGrob("R=0.988 \np = 2.035e-09", x=0.1,  y=0.75, hjust=0,
                           gp=gpar(col="darkgreen", fontsize=10, fontface="italic")))
grob2 <- grobTree(textGrob("R=0.998 \np = 1.555e-13", x=0.55,  y=0.1, hjust=0,
                           gp=gpar(col="darkred", fontsize=10, fontface="italic")))

png("./plots/105LOD_ddPCR_LAMP5AS1_forward1.png", width = 10, height = 8, res = 300, units = "cm")
p1 <- ggplot(lod[idx, ], aes(x=log10(spike), y=log10(concentration*20), color=gene)) + 
  geom_point() + 
  geom_smooth(method=lm)+scale_color_brewer(palette="Dark2") +
  labs(title="Nerve Marker",
       x="SpikeIn (Copies) \n Log10-scaled", y = "ddPCR readout (Copies) \n Log10-scaled")+
  theme_classic() + xlim(0.5,4) + ylim(0.5,4)+ 
  annotation_custom(grob1) + annotation_custom(grob2)
dev.off()



idx <- grep("LAMP5AS1_F2", ddpcr$target)
lod <- ddpcr[idx,]
idx <- c(1:4,7:10, 13:16,19:22,25:28,31:34)
# Add the regression line
ggplot(lod[idx, ], aes(x=log10(spike), y=log10(concentration*20), color=gene)) + 
  geom_point() + 
  geom_smooth(method=lm)+scale_color_brewer(palette="Dark2") +
  labs(title="Forward2",
       x="SpikeIn (Copies) LOG", y = "ddPCR readout (Copies) LOG")+
  theme_classic() + xlim(0.5,4) + ylim(0.5,4)
dev.off()
cor.test(log10(lod[c(1:4,13:16,25:28),]$concentration*20), log10(lod[c(1:4,13:16,25:28),]$spike), method = "pearson")
cor.test(log10(lod[c(7:10,19:22,31:34),]$concentration*20), log10(lod[c(7:10,19:22,31:34),]$spike), method = "pearson")
# Create a text
grob1 <- grobTree(textGrob("R=0.987 \np = 2.284e-09", x=0.1,  y=0.75, hjust=0,
                           gp=gpar(col="darkgreen", fontsize=10, fontface="italic")))
grob2 <- grobTree(textGrob("R=0.995 \np = 1.702e-11", x=0.55,  y=0.1, hjust=0,
                           gp=gpar(col="darkred", fontsize=10, fontface="italic")))

png("./plots/105LOD_ddPCR_LAMP5AS1_forward2.png", width = 10, height = 8, res = 300, units = "cm")
ggplot(lod[idx, ], aes(x=log10(spike), y=log10(concentration*20), color=gene)) + 
  geom_point() + 
  geom_smooth(method=lm)+scale_color_brewer(palette="Dark2") +
  labs(title="LAMP5AS1_f2",
       x="SpikeIn (Copies) LOG", y = "ddPCR readout (Copies) LOG")+
  theme_classic() + xlim(0.5,4) + ylim(0.5,4)+ 
  annotation_custom(grob1) + annotation_custom(grob2)
dev.off()



idx <- grep("LAMP5AS1_F1", ddpcr$target)
cv <- ddpcr[idx,]
cv <- cv[cv$gene == "LAMP5AS1_F1",]
sd(cv$concentration[c(1,7,13)])
sd(cv$concentration[c(1,7,13)])/mean(cv$concentration[c(1,7,13)])*100
sd(cv$concentration[c(2,8,14)])
sd(cv$concentration[c(2,8,14)])/mean(cv$concentration[c(2,8,14)])*100
sd(cv$concentration[c(3,9,15)])
sd(cv$concentration[c(3,9,15)])/mean(cv$concentration[c(3,9,15)])*100
sd(cv$concentration[c(4,10,16)])
sd(cv$concentration[c(4,10,16)])/mean(cv$concentration[c(4,10,16)])*100
sd(cv$concentration[c(5,11,17)])
sd(cv$concentration[c(5,11,17)])/mean(cv$concentration[c(5,11,17)])*100

idx <- grep("LAMP5AS1_F2", ddpcr$target)
cv <- ddpcr[idx,]
cv <- cv[cv$gene == "LAMP5AS1_F2",]
sd(cv$concentration[c(1,7,13)])
sd(cv$concentration[c(1,7,13)])/mean(cv$concentration[c(1,7,13)])*100
sd(cv$concentration[c(2,8,14)])
sd(cv$concentration[c(2,8,14)])/mean(cv$concentration[c(2,8,14)])*100
sd(cv$concentration[c(3,9,15)])
sd(cv$concentration[c(3,9,15)])/mean(cv$concentration[c(3,9,15)])*100
sd(cv$concentration[c(4,10,16)])
sd(cv$concentration[c(4,10,16)])/mean(cv$concentration[c(4,10,16)])*100
sd(cv$concentration[c(5,11,17)])
sd(cv$concentration[c(5,11,17)])/mean(cv$concentration[c(5,11,17)])*100





ddpcr <- read.csv("./ddPCR_LOD_sa.csv")
###SA--GSC1b####
idx <- grep("GSC1b", ddpcr$target)
lod <- ddpcr[idx,]
idx <- c(1:4,7:10, 13:16,19:22,25:28,31:34)
# Add the regression line
ggplot(lod[idx, ], aes(x=log10(spike), y=log10(concentration*20), color=gene)) + 
  geom_point() + 
  geom_smooth(method=lm)+scale_color_brewer(palette="Dark2") +
  labs(title="GSC1b",
       x="SpikeIn (Copies) LOG", y = "ddPCR readout (Copies) LOG")+
  theme_classic() + xlim(0.5,4) + ylim(0.5,4)
dev.off()
tmp <- lod[c(1:4,13:16,25:28),]$concentration*20
tmp[4] <- 1
cor.test(log10(tmp), log10(lod[c(1:4,13:16,25:28),]$spike), method = "pearson")
cor.test(log10(lod[c(7:10,19:22,31:34),]$concentration*20), log10(lod[c(7:10,19:22,31:34),]$spike), method = "pearson")
# Create a text
grob1 <- grobTree(textGrob("R=0.993 \np = 1.019e-10", x=0.1,  y=0.75, hjust=0,
                           gp=gpar(col="darkgreen", fontsize=10, fontface="italic")))
grob2 <- grobTree(textGrob("R=0.952 \np = 1.772e-06", x=0.55,  y=0.1, hjust=0,
                           gp=gpar(col="darkred", fontsize=10, fontface="italic")))

png("./plots/105LOD_ddPCR_GSC1b.png", width = 10, height = 8, res = 300, units = "cm")
p2 <- ggplot(lod[idx, ], aes(x=log10(spike), y=log10(concentration*20), color=gene)) + 
  geom_point() + 
  geom_smooth(method=lm)+scale_color_brewer(palette="Dark2") +
  labs(title="SA Marker",
       x="SpikeIn (Copies) \n Log10-scaled", y = "ddPCR readout (Copies) \n Log10-scaled")+
  theme_classic() + xlim(0.5,4) + ylim(0.5,4)+ 
  annotation_custom(grob1) + annotation_custom(grob2)
dev.off()


idx <- grep("GSC1b", ddpcr$target)
cv <- ddpcr[idx,]
cv <- cv[cv$gene == "GSC1b",]
sd(cv$concentration[c(1,7,13)])
sd(cv$concentration[c(1,7,13)])/mean(cv$concentration[c(1,7,13)])*100
sd(cv$concentration[c(2,8,14)])
sd(cv$concentration[c(2,8,14)])/mean(cv$concentration[c(2,8,14)])*100
sd(cv$concentration[c(3,9,15)])
sd(cv$concentration[c(3,9,15)])/mean(cv$concentration[c(3,9,15)])*100
sd(cv$concentration[c(4,10,16)])
sd(cv$concentration[c(4,10,16)])/mean(cv$concentration[c(4,10,16)])*100
sd(cv$concentration[c(5,11,17)])
sd(cv$concentration[c(5,11,17)])/mean(cv$concentration[c(5,11,17)])*100





ddpcr <- read.csv("./ddPCR_LOD_kera.csv")
###Keratinocyte--IFFO1_1####
idx <- grep("IFFO1_1", ddpcr$target)
lod <- ddpcr[idx,]
idx <- c(1:4,7:10, 13:16,19:22,25:28,31:34)
# Add the regression line
ggplot(lod[idx, ], aes(x=log10(spike), y=log10(concentration*20), color=gene)) + 
  geom_point() + 
  geom_smooth(method=lm)+scale_color_brewer(palette="Dark2") +
  labs(title="GSC1b",
       x="SpikeIn (Copies) LOG", y = "ddPCR readout (Copies) LOG")+
  theme_classic() + xlim(0.5,4) + ylim(0.5,4)
dev.off()
cor.test(log10(lod[c(1:4,13:16,25:28),]$concentration*20), log10(lod[c(1:4,13:16,25:28),]$spike), method = "pearson")
cor.test(log10(lod[c(7:10,19:22,31:34),]$concentration*20), log10(lod[c(7:10,19:22,31:34),]$spike), method = "pearson")
# Create a text
grob1 <- grobTree(textGrob("R=0.996 \np = 5.545e-12", x=0.1,  y=0.75, hjust=0,
                           gp=gpar(col="darkgreen", fontsize=10, fontface="italic")))
grob2 <- grobTree(textGrob("R=0.987 \np = 3.115e-09", x=0.55,  y=0.1, hjust=0,
                           gp=gpar(col="darkred", fontsize=10, fontface="italic")))

png("./plots/105LOD_ddPCR_IFFO1_1.png", width = 10, height = 8, res = 300, units = "cm")
p3 <- ggplot(lod[idx, ], aes(x=log10(spike), y=log10(concentration*20), color=gene)) + 
  geom_point() + 
  geom_smooth(method=lm)+scale_color_brewer(palette="Dark2") +
  labs(title="Keratinocyte Marker",
       x="SpikeIn (Copies) \n Log10-scaled", y = "ddPCR readout (Copies) \n Log10-scaled")+
  theme_classic() + xlim(0.5,4) + ylim(0.5,4)+ 
  annotation_custom(grob1) + annotation_custom(grob2)
dev.off()


idx <- grep("IFFO1_1", ddpcr$target)
cv <- ddpcr[idx,]
cv <- cv[cv$gene == "IFFO1_1",]
sd(cv$concentration[c(1,7,13)])
sd(cv$concentration[c(1,7,13)])/mean(cv$concentration[c(1,7,13)])*100
sd(cv$concentration[c(2,8,14)])
sd(cv$concentration[c(2,8,14)])/mean(cv$concentration[c(2,8,14)])*100
sd(cv$concentration[c(3,9,15)])
sd(cv$concentration[c(3,9,15)])/mean(cv$concentration[c(3,9,15)])*100
sd(cv$concentration[c(4,10,16)])
sd(cv$concentration[c(4,10,16)])/mean(cv$concentration[c(4,10,16)])*100
sd(cv$concentration[c(5,11,17)])
sd(cv$concentration[c(5,11,17)])/mean(cv$concentration[c(5,11,17)])*100


ddpcr <- read.csv("./ddPCR_LOD_fibro.csv")
###Fibroblast--SCARF2_1####
idx <- grep("SCARF2_1", ddpcr$target)
lod <- ddpcr[idx,]
idx <- c(1:4,7:10, 13:16,19:22,25:28,31:34)
# Add the regression line
ggplot(lod[idx, ], aes(x=log10(spike), y=log10(concentration*20), color=gene)) + 
  geom_point() + 
  geom_smooth(method=lm)+scale_color_brewer(palette="Dark2") +
  labs(title="SCARF2_1",
       x="SpikeIn (Copies) LOG", y = "ddPCR readout (Copies) LOG")+
  theme_classic() + xlim(0.5,4) + ylim(0.5,4)
dev.off()
cor.test(log10(lod[c(1:4,13:16,25:28),]$concentration*20), log10(lod[c(1:4,13:16,25:28),]$spike), method = "pearson")
cor.test(log10(lod[c(7:10,19:22,31:34),]$concentration*20), log10(lod[c(7:10,19:22,31:34),]$spike), method = "pearson")
# Create a text
grob1 <- grobTree(textGrob("R=0.998 \np = 9.189e-14", x=0.1,  y=0.75, hjust=0,
                           gp=gpar(col="darkgreen", fontsize=10, fontface="italic")))
grob2 <- grobTree(textGrob("R=0.993 \np = 1.114e-10", x=0.55,  y=0.1, hjust=0,
                           gp=gpar(col="darkred", fontsize=10, fontface="italic")))

png("./plots/105LOD_ddPCR_SCARF2_1.png", width = 10, height = 8, res = 300, units = "cm")
p4 <- ggplot(lod[idx, ], aes(x=log10(spike), y=log10(concentration*20), color=gene)) + 
  geom_point() + 
  geom_smooth(method=lm)+scale_color_brewer(palette="Dark2") +
  labs(title="Fibroblast Marker",
       x="SpikeIn (Copies) \n Log10-scaled", y = "ddPCR readout (Copies) \n Log10-scaled")+
  theme_classic() + xlim(0.5,4) + ylim(0.5,4)+ 
  annotation_custom(grob1) + annotation_custom(grob2)
dev.off()


idx <- grep("SCARF2_1", ddpcr$target)
cv <- ddpcr[idx,]
cv <- cv[cv$gene == "SCARF2_1",]
sd(cv$concentration[c(1,7,13)])
sd(cv$concentration[c(1,7,13)])/mean(cv$concentration[c(1,7,13)])*100
sd(cv$concentration[c(2,8,14)])
sd(cv$concentration[c(2,8,14)])/mean(cv$concentration[c(2,8,14)])*100
sd(cv$concentration[c(3,9,15)])
sd(cv$concentration[c(3,9,15)])/mean(cv$concentration[c(3,9,15)])*100
sd(cv$concentration[c(4,10,16)])
sd(cv$concentration[c(4,10,16)])/mean(cv$concentration[c(4,10,16)])*100
sd(cv$concentration[c(5,11,17)])
sd(cv$concentration[c(5,11,17)])/mean(cv$concentration[c(5,11,17)])*100



ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)



###ddPCR blister fluid ####
load(file = "./phenoData_EMseq_deconvolution_percentage.rds")
phenoData
phenoData$LAMP_padj <- as.numeric(phenoData$LAMP_padj)
phenoData$SCARF2_padj <- as.numeric(phenoData$SCARF2_padj)
phenoData$COL2A1_padj <- as.numeric(phenoData$COL2A1_padj)
phenoData$ZFAM_padj <- as.numeric(phenoData$ZFAM_padj)
phenoData$FOXF1_padj <- as.numeric(phenoData$FOXF1_padj)
phenoData$IFFO1_padj <- as.numeric(phenoData$IFFO1_padj)
phenoData$Nerve_perc <- phenoData$LAMP_padj/phenoData$COL2A1_padj
phenoData$SA_perc <- phenoData$FOXF1_padj  /phenoData$COL2A1_padj
phenoData$Endo_perc <- phenoData$ZFAM_padj/phenoData$COL2A1_padj
phenoData$Kera_perc <- phenoData$IFFO1_padj/phenoData$COL2A1_padj 

library(ggplot2)
library(ggpubr)
library(ggfortify)
####PCA####
pca.df <- phenoData[,c( "SEX",    "graft", "DEPTH",
                        "DAYS.TO.RE.EP",  "AGE_months","TSBA", 
                       "CONCENTRATION.by.BCA.mg.mL.", "Nerve_perc",     "SA_perc" ,     "Endo_perc", "Kera_perc", 
                       "LAMP_padj", "FOXF1_padj", "ZFAM_padj", "IFFO1_padj")]
colnames(pca.df)[12:15] <- c("Nerve", "SA", "Endo", "Kera")
pca.df$SEX <- factor(pca.df$SEX, levels = c("M", "F"), labels =c("M", "F"))
pca.df$TSBA <- gsub("<", "", pca.df$TSBA)
pca.df$TSBA <- as.numeric(pca.df$TSBA)
pca.df$graft <- factor(pca.df$graft, levels = c("No", "Yes"), labels =c("No", "Yes"))
colnames(pca.df) <- gsub("CONCENTRATION.by.BCA.mg.mL.", "ProteinConc.", colnames(pca.df))
pca.df$DEPTH <- gsub("Superficial", 1, pca.df$DEPTH )
pca.df$DEPTH <- gsub("Deep", 2, pca.df$DEPTH )
pca.df$DEPTH <- gsub("Full", 3, pca.df$DEPTH )
pca.df$DEPTH <- as.numeric(pca.df$DEPTH)
pca.df$ProteinConc. <- as.numeric(pca.df$ProteinConc.)
pca.df <- pca.df[complete.cases(pca.df),]
prcomp_obj <- prcomp(pca.df[,3:11],  scale. = T)  
autoplot(prcomp_obj, data = pca.df, #colour = 'DEPTH',
              loadings = TRUE, loadings.colour = 'blue',
              loadings.label = TRUE, loadings.label.size = 4)
prcomp_obj <- prcomp(pca.df[,c(3:7,12:15)],  scale. = T)  
autoplot(prcomp_obj, data = pca.df, #colour = 'DEPTH',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 4)



####Correlation ####
phenoData$CONCENTRATION.by.BCA.mg.mL. <- as.numeric(phenoData$CONCENTRATION.by.BCA.mg.mL.)

conc <- read.csv("./data/blister_concentration.csv")
conc <- conc[grep("QUT", conc$id),1:2]

phenoData$cfDNA_Conc <-  conc[match(rownames(phenoData), conc$id), "Conc"]
phenoData$cfDNA_Conc[3]
phenoData$cfDNA_Conc[3] <-  55.6

ggplot(data = phenoData[,], aes(x = cfDNA_Conc, y = COL2A1_padj)) + 
  geom_point(aes(colour = factor(DEPTH), size = 2)) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  geom_smooth(method = "lm", se = F) + 
  ggtitle("Total Conc.") +theme_classic()+ stat_cor(method = "pearson", label.x = 20, label.y = 100)
ggplot(data = phenoData[,], aes(x = Endo_perc, y = EC)) + 
  geom_point(aes(colour = factor(DEPTH), size = 2)) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  geom_smooth(method = "lm", se = F) + 
  ggtitle("endo Conc.")+ylim(c(0,0.1)) +theme_classic() +theme_classic()+ stat_cor(method = "pearson", label.x = 0, label.y = 0.075)

p1 <- ggplot(data = phenoData[,], aes(x = LAMP_padj, y = Nerve )) + 
  geom_point(aes(colour = factor(DEPTH), size = 2)) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  geom_smooth(method = "lm", se = F) + 
  ggtitle("Nerve")  + stat_cor(method = "pearson", label.x = 0, label.y = 0.25) +
  xlab("Nerve Copies (ddPCR)") + ylab("Nerve Percentage (EMseq)")+
  scale_y_continuous(labels = scales::percent_format(), limits = c(0.05, 0.25)) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
p2 <- ggplot(data = phenoData[,], aes(x = FOXF1_padj , y = SA)) + 
  geom_point(aes(colour = factor(DEPTH), size = 2)) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  geom_smooth(method = "lm", se = F)  + 
  ggtitle(" ")  + stat_cor(method = "pearson", label.x = 0, label.y = 0.05) +
  xlab("SA Copies (ddPCR)") + ylab("SA Percentage (EMseq)")+
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.05)) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
p3 <- ggplot(data = phenoData[,], aes(x = Kera_perc , y = Keratinocyte)) + 
  geom_point(aes(colour = factor(DEPTH), size = 2)) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  geom_smooth(method = "lm", se = F)  + 
  ggtitle(" ")  + stat_cor(method = "pearson", label.x = 0, label.y = 0.05) +
  xlab("Keratinocyte Copies (ddPCR)") + ylab("Keratinocyte Percentage (EMseq)")+
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.05)) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
ggarrange(p1, p2, nrow = 1, common.legend = T)


####Nerve ####
p2 <- ggplot(phenoData , aes(x=DEPTH, y=Nerve_perc, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  geom_text(aes(colour = DEPTH), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_colour_manual(values=c("#444444", "#A77305", "#0690A1"))+
  ylim(c(0,0.04)) + ggtitle("Nerve ddPCR")+ theme(legend.position="none")+
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less"),
                     comparisons = list(c( "Superficial", "Deep"), c("Deep", "Full")),
                     label = "p.signif")  # Show asterisks
p1 <- ggplot(phenoData , aes(x=DEPTH, y=Nerve, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  geom_text(aes(colour = DEPTH), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_colour_manual(values=c("#444444", "#A77305", "#0690A1"))+
  ylim(c(0.05,0.25)) + ggtitle("Nerve in silico")+ theme(legend.position="none")
p3 <-  ggplot(phenoData , aes(x=DEPTH, y=LAMP_padj, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  geom_text(aes(colour = DEPTH), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_colour_manual(values=c("#444444", "#A77305", "#0690A1"))+
   ggtitle("Nerve ddPCR conc") + theme(legend.position="none")+
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less"),
                     comparisons = list(c( "Superficial", "Deep"), c("Deep", "Full")),
                     label = "p.signif")  # Show asterisks
ggarrange(p1, p2, p3, ncol = 3, nrow = 1)
                    
####SA ####
p2 <- ggplot(phenoData , aes(x=DEPTH, y=SA_perc, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  geom_text(aes(colour = DEPTH), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_colour_manual(values=c("#444444", "#A77305", "#0690A1"))+
  ylim(c(0,0.01)) + ggtitle("SA ddPCR")+ theme(legend.position="none")+
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less"),
                     comparisons = list(c( "Superficial", "Deep"), c("Deep", "Full")),
                     label = "p.signif")  # Show asterisks
p1 <- ggplot(phenoData , aes(x=DEPTH, y=SA, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  geom_text(aes(colour = DEPTH), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_colour_manual(values=c("#444444", "#A77305", "#0690A1"))+
  ylim(c(0,0.04)) + ggtitle("SA in silico")+ theme(legend.position="none")
p3 <-  ggplot(phenoData , aes(x=DEPTH, y=FOXF1_padj, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  geom_text(aes(colour = DEPTH), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_colour_manual(values=c("#444444", "#A77305", "#0690A1"))+
  ggtitle("SA ddPCR conc") + theme(legend.position="none")+
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less"),
                     comparisons = list(c( "Superficial", "Deep"), c("Deep", "Full")),
                     label = "p.signif")  # Show asterisks
ggarrange(p1, p2, p3, 
          #labels = c("A", "B"),
          ncol = 3, nrow = 1)


####VE ####
p2 <- ggplot(phenoData , aes(x=DEPTH, y=Endo_perc, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  geom_text(aes(colour = DEPTH), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_colour_manual(values=c("#444444", "#A77305", "#0690A1"))+
  ylim(c(0,0.01)) + ggtitle("VE ddPCR")+ theme(legend.position="none")+
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less"),
                     comparisons = list(c( "Superficial", "Deep"), c("Deep", "Full")),
                     label = "p.signif")  # Show asterisks
p1 <- ggplot(phenoData , aes(x=DEPTH, y=EC, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  geom_text(aes(colour = DEPTH), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_colour_manual(values=c("#444444", "#A77305", "#0690A1"))+
  ylim(c(0,0.1)) + ggtitle("VE in silico")+ theme(legend.position="none")
p3 <-  ggplot(phenoData , aes(x=DEPTH, y=ZFPM, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  geom_text(aes(colour = DEPTH), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_colour_manual(values=c("#444444", "#A77305", "#0690A1"))+
  ggtitle("VE ddPCR conc") + theme(legend.position="none")+
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less"),
                     comparisons = list(c( "Superficial", "Deep"), c("Deep", "Full")),
                     label = "p.signif")  # Show asterisks
ggarrange(p1, p2, p3, 
          #labels = c("A", "B"),
          ncol = 3, nrow = 1)


####Epidermal ####
p1 <- ggplot(phenoData , aes(x=DEPTH, y=epidermal_perc, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  geom_text(aes(colour = DEPTH), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_colour_manual(values=c("#444444", "#A77305", "#0690A1"))+
  ylim(c(0,1)) + ggtitle("Epidermal ddPCR")+ theme(legend.position="none")
p2 <- ggplot(phenoData , aes(x=DEPTH, y=IFFO1_padj, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  geom_text(aes(colour = DEPTH), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_colour_manual(values=c("#444444", "#A77305", "#0690A1"))+
  ylim(c(0,1)) + ggtitle("Epidermal in silico")+ theme(legend.position="none")
ggarrange(p1, p2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)


p1 <- ggplot(phenoData , aes(x=DEPTH, y=Nerve_perc, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  ylim(c(0,0.04)) + theme(legend.position="none")+
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less"),
                     comparisons = list(c( "Superficial", "Deep"), c("Deep", "Full")),
                     label = "p.signif") + # Show asterisks
  labs(x = "Group",                                              # labelling x axis
       y = "Ratio",                                        # labeling y axis
       title = "Nerve",        # title
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
p2 <- ggplot(phenoData , aes(x=DEPTH, y=SA_perc, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  ylim(c(0,0.01)) + theme(legend.position="none")+
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less"),
                     comparisons = list(c( "Superficial", "Deep"), c("Deep", "Full")),
                     label = "p.signif") + # Show asterisks
  labs(x = "Group",                                              # labelling x axis
       y = "Ratio",                                        # labeling y axis
       title = "SA",        # title
       fill = "Group") +                               # legend
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.01)) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
p3 <- ggplot(phenoData , aes(x=DEPTH, y=Endo_perc, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  ylim(c(0,0.01)) + theme(legend.position="none")+
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less"),
                     comparisons = list(c( "Superficial", "Deep"), c("Deep", "Full")),
                     label = "p.signif") + # Show asterisks
  labs(x = "Group",                                              # labelling x axis
       y = "Ratio",                                        # labeling y axis
       title = "EC",        # title
       fill = "Group") +                               # legend
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.01)) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
ggarrange(p1, p2, p3,
          ncol = 3, nrow = 1)





####Re-epitheliation  ####
p1 <- ggplot(data = phenoData[which(phenoData$DAYS.TO.RE.EP<300),], aes(x = Nerve_perc, y = DAYS.TO.RE.EP)) + 
  geom_point(aes(colour = factor(DEPTH), size = 2)) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  geom_smooth(method = "lm", se = F) + 
   stat_cor(method = "pearson", label.x = 0, label.y = 200)+
  xlab("Nerve Ratio") + ylab("Days to Re-epitheliation")+
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )

p2 <- ggplot(data = phenoData[which(phenoData$DAYS.TO.RE.EP<300),], aes(x = SA_perc, y = DAYS.TO.RE.EP)) + 
  geom_point(aes(colour = factor(DEPTH), size = 2)) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  geom_smooth(method = "lm", se = F) + 
   stat_cor(method = "pearson", label.x = 0, label.y = 200)+
  xlab("SA Ratio") + ylab("Days to Re-epitheliation")+
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )

p3 <- ggplot(data = phenoData[which(phenoData$DAYS.TO.RE.EP<300),], aes(x = Endo_perc, y = DAYS.TO.RE.EP)) + 
  geom_point(aes(colour = factor(DEPTH), size = 2)) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  geom_smooth(method = "lm", se = F) + 
   stat_cor(method = "pearson", label.x = 0, label.y = 200)+
  xlab("EC Ratio") + ylab("Days to Re-epitheliation")+
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
ggplot(data = phenoData[which(phenoData$DAYS.TO.RE.EP<300),], aes(x = Endo_perc, y = DAYS.TO.RE.EP)) + 
  geom_point(aes(colour = factor(DEPTH), size = 2)) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  geom_smooth(method = "lm", se = F) + 
  ggtitle("SA Conc. vs Days to Re-epitheliation") + stat_cor(method = "pearson", label.x = 0, label.y = 200)+
  xlab("EC Ratio") + ylab("Nerve Percentage (EMseq)")+
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
ggarrange(p1, p2, p3,
          ncol = 3, nrow = 1 , common.legend = T)


pal.depth <- c("#999999", "#E69F00", "#56B4E9")
names(pal.depth) <-c( "Superficial", "Deep", "Full")
reepi.df <- phenoData[which(phenoData$DAYS.TO.RE.EP<300),]
p1 <- ggplot(data = reepi.df[grep("Superficial", reepi.df$DEPTH),], aes(x = Nerve_perc, y = DAYS.TO.RE.EP)) + 
  geom_point(aes(colour = factor(DEPTH), size = 2)) + scale_color_manual(values=pal.depth)+
  geom_smooth(method = "lm", se = F) + 
   stat_cor(method = "pearson", label.x = 0, label.y = 50)+
  xlab("Nerve Ratio") + ylab("Days to Re-epitheliation")+
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
p2 <- ggplot(data = reepi.df[grep("Superficial", reepi.df$DEPTH, invert = T),], aes(x = Nerve_perc, y = DAYS.TO.RE.EP)) + 
  geom_point(aes(colour = factor(DEPTH), size = 2)) + scale_color_manual(values=pal.depth)+
  geom_smooth(method = "lm", se = F) + 
   stat_cor(method = "pearson", label.x = 0, label.y = 200)+
  xlab("Nerve Ratio") + ylab("Days to Re-epitheliation")+
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
p3 <- ggplot(data = reepi.df[grep("Superficial", reepi.df$DEPTH),], aes(x = SA_perc, y = DAYS.TO.RE.EP)) + 
  geom_point(aes(colour = factor(DEPTH), size = 2)) + scale_color_manual(values=pal.depth)+
  geom_smooth(method = "lm", se = F) + 
   stat_cor(method = "pearson", label.x = 0, label.y = 50)+
  xlab("SA Ratio") + ylab("Days to Re-epitheliation")+
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
p4 <- ggplot(data = reepi.df[grep("Superficial", reepi.df$DEPTH, invert = T),], aes(x = SA_perc, y = DAYS.TO.RE.EP)) + 
  geom_point(aes(colour = factor(DEPTH), size = 2)) + scale_color_manual(values=pal.depth)+
  geom_smooth(method = "lm", se = F) + 
   stat_cor(method = "pearson", label.x = 0, label.y = 200)+
  xlab("SA Ratio") + ylab("Days to Re-epitheliation")+
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
p5 <- ggplot(data = reepi.df[grep("Superficial", reepi.df$DEPTH),], aes(x = Endo_perc, y = DAYS.TO.RE.EP)) + 
  geom_point(aes(colour = factor(DEPTH), size = 2)) + scale_color_manual(values=pal.depth)+
  geom_smooth(method = "lm", se = F) + 
  stat_cor(method = "pearson", label.x = 0, label.y = 50)+
  xlab("EC Ratio") + ylab("Days to Re-epitheliation")+
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
p6 <- ggplot(data = reepi.df[grep("Superficial", reepi.df$DEPTH, invert = T),], aes(x = Endo_perc, y = DAYS.TO.RE.EP)) + 
  geom_point(aes(colour = factor(DEPTH), size = 2)) + scale_color_manual(values=pal.depth)+
  geom_smooth(method = "lm", se = F) + 
  stat_cor(method = "pearson", label.x = 0, label.y = 200)+
  xlab("EC Ratio") + ylab("Days to Re-epitheliation")+
  theme(
    axis.text.x = element_text(angle = 0,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
ggarrange(p1, p3, p5, p2, p4, p6, ncol = 3, nrow = 2, common.legend = T)


####ROC#####
table(phenoData$DEPTH, phenoData$graft)
table(phenoData$TSBA, phenoData$graft)
library(pROC)
phenoData$outcome <- gsub("full_thickness", 1, phenoData$DEPTH)
phenoData$outcome  <- gsub("deep|superficial", 0,phenoData$outcome)
modelroc <- roc(phenoData$outcome, phenoData$nerve_perc)
plot(modelroc, print.auc=TRUE, auc.polygon=F, grid=c(0.1, 0.2),
     #grid.col=c("green", "red"),
     max.auc.polygon=F,
     auc.polygon.col="skyblue", 
     print.thres=F,
     main = "ROC for prediction of full thickness burn \n by Nerve Perc.")
dev.off()

phenoData$outcome <- gsub("Yes", 1, phenoData$graft)
phenoData$outcome  <- gsub("No", 0,phenoData$outcome)
modelroc <- roc(phenoData$outcome, phenoData$nerve_perc)
plot(modelroc, print.auc=TRUE, auc.polygon=F, grid=c(0.1, 0.2),
     #grid.col=c("green", "red"),
     max.auc.polygon=F,
     auc.polygon.col="skyblue", 
     print.thres=F,
     main = "ROC for prediction of graft \n by Nerve Perc.")
dev.off()

phenoData$tsba_number <- phenoData$TSBA
phenoData$tsba_number <- gsub("<1", 1, phenoData$tsba_number)
phenoData$tsba_number <- as.numeric(phenoData$tsba_number)
modelroc <- roc(phenoData$outcome, phenoData$TSBA)
plot(modelroc, print.auc=TRUE, auc.polygon=F, grid=c(0.1, 0.2),
     #grid.col=c("green", "red"),
     max.auc.polygon=F,
     auc.polygon.col="skyblue", 
     print.thres=F,
     main = "ROC for prediction of graft \n by Nerve Perc.")
dev.off()




###color by graft ####
ggplot(phenoData , aes(x=DEPTH, y=nerve_perc, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  geom_text(aes(colour = graft), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_colour_manual(values=c("#444444", "darkred"))+
  ylim(c(0,0.04)) + ggtitle("Nerve ddPCR")+ theme(legend.position="none")
ggplot(phenoData , aes(x=DEPTH, y=adi_perc, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  geom_text(aes(colour = graft), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_colour_manual(values=c("#444444", "darkred"))+
  ylim(c(0,0.01)) + ggtitle("SA ddPCR")+ theme(legend.position="none")



ggplot(phenoData , aes(x=DEPTH, y=LAMP, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  geom_text(aes(colour = graft), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_colour_manual(values=c("#444444", "darkred"))+
  ggtitle("Nerve ddPCR conc") + theme(legend.position="none")
ggplot(phenoData , aes(x=DEPTH, y=FOXF1, fill = DEPTH, label = SAMPLE_ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=-1, notch=FALSE)+
  geom_text(aes(colour = graft), fontface = "bold",position=position_jitter(width=0.12,height=0))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_colour_manual(values=c("#444444", "darkred"))+
  ggtitle("SA ddPCR conc") + theme(legend.position="none")

