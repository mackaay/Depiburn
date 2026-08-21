library(tidyverse)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(limma)
setwd("./DEpiBurn_EMseq/")

# 1. Load All Sample TSVs#####
file_list <- list.files(path = "/datasets/work/hb-burns-meth/work/Data/level2/motif_ONT/",pattern = "*_motif_frequencies_N.tsv", full.names = T)
names(file_list) <- gsub("_motif_frequencies_N.tsv", "", basename(file_list) )

# Combine into one long dataframe
obs_data <- map_df(file_list, read_tsv, .id = "Sample_ID") %>%
  select(Sample_ID, Motif, Frequency)

# Spread to wide format (Motifs as columns, Samples as rows)
motif_matrix <- obs_data %>%
  pivot_wider(names_from = Motif, values_from = Frequency)


# 2. Load Expected Frequencies (from BBMap k-mer calculation)######
# Assuming a CSV with columns 'Motif' and 'Expected_Freq'
expected_df <- read.csv("./human_genom_4mer_expected.csv", header = F)
expected_df <- read.csv("./human_genome_4mer_expected_full.csv", header = T)
sum(expected_df$Expected_Freq)
expected_df %>% arrange(desc(Expected_Freq)) %>% head()
expected_df <- expected_df[,c(1,3)]
colnames(expected_df) <- c("Motif", "Expected_Freq")


# Calculate mean observed frequency across all samples
mean_obs <- obs_data %>%
  group_by(Motif) %>%
  summarise(Observed_Freq = mean(Frequency)) %>%
  left_join(expected_df, by = "Motif")

# Scatterplot
ggplot(mean_obs, aes(x = Expected_Freq, y = Observed_Freq)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Expected vs. Observed 4-mer Frequencies",
       x = "Genomic Expected Frequency", y = "cfDNA Observed Frequency")

mean_obs[which(mean_obs$Observed_Freq >0.0075 & mean_obs$Expected_Freq <0.01),]
mean_obs[which(mean_obs$Observed_Freq < 0.0075 & mean_obs$Expected_Freq > 0.0075),]


# 3. PCA######
sample_data <- data.frame(
  Sample_ID = names(file_list),
  BAM_File = c(file_list),
  Group = c("Superficial", "Deep", "Deep", "Superficial", "Superficial", "Deep", "Superficial", "Deep", "Superficial", "Deep", "Superficial", "Full"),
  stringsAsFactors = FALSE
)
rownames(sample_data) <- sample_data$Sample_ID
# Prepare matrix for PCA (exclude metadata)
pca_input <- motif_matrix %>% column_to_rownames("Sample_ID")

pal.depth <- c("Superficial"= "#999999","Deep"= "#E69F00","Full"= "#56B4E9")
plotMDS(t(pca_input), pch = 19, col = pal.depth[factor(sample_data$Group)])


