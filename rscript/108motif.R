library(tidyverse)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(limma)
setwd("./DEpiBurn_EMseq/")

# 1. Load All Sample TSVs#####
file_list <- list.files(path = "/datasets/work/hb-burns-meth/work/Data/level2/motif/",pattern = "*_motif_frequencies_N.tsv", full.names = T)
names(file_list) <- gsub("_motif_frequencies_N.tsv", "", basename(file_list) )

# Combine into one long dataframe
obs_data <- map_df(file_list, read_tsv, .id = "Sample_ID") %>%
  select(Sample_ID, Motif, Frequency)

for (i in 1:length(file_list)) {
  read.delim(file = file_list[i], header = T)
}

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


# 3. PCA #######
sample_data <- data.frame(
  Sample_ID = names(file_list),
  Group = c("Full", "Deep", "Deep", "Deep", "Deep", "Full", "Superficial", "Superficial", #207
            "Superficial", "Superficial", "Deep", "Full", 
            "Superficial", "Deep", "Superficial", "Superficial", "Deep", "Superficial", 
            "Deep", "Superficial", "Superficial", "Full"),
  stringsAsFactors = FALSE
)
rownames(sample_data) <- sample_data$Sample_ID
# Prepare matrix for PCA (exclude metadata)
pca_input <- motif_matrix %>% column_to_rownames("Sample_ID")
res.pca <- PCA(pca_input, graph = FALSE)

# Plot PCA
fviz_pca_ind(res.pca, 
             label = "none", 
             habillage = sample_data$Group, # Requires a metadata dataframe
             addEllipses = TRUE, 
             title = "PCA of cfDNA End-Motif Usage")

pal.depth <- c("Superficial"= "#999999","Deep"= "#E69F00","Full"= "#56B4E9")
plotMDS(t(pca_input), pch = 19, col = pal.depth[factor(sample_data$Group)])




# Merge motif data with clinical metadata
library(tidyverse)
library(rstatix) # High-value library for multi-group stats
analysis_df <- motif_matrix %>%
  left_join(sample_data, by = "Sample_ID") # metadata must have 'Group' (Sepsis/Control)

valid_motifs <- analysis_df %>%
  pivot_longer(cols = all_of(colnames(pca_input)), 
               names_to = "Motif", values_to = "Freq") %>%
  group_by(Motif) %>%
  summarise(n_present = sum(Freq > 0), .groups = "drop") %>%
  filter(n_present >= 3) %>% # Adjust this number based on your smallest group size
  pull(Motif)

results_multi <- analysis_df %>%
  pivot_longer(cols = all_of(colnames(valid_motifs)), 
               names_to = "Motif", 
               values_to = "Frequency") %>%
  # Filter out any NAs that might be causing group drops
  filter(!is.na(Frequency), !is.na(Group)) %>%
  filter(n_distinct(Group[Frequency > 0]) > 1) %>%
  group_by(Motif) %>%
  # Perform Kruskal-Wallis test (Is there any difference among the 3 groups?)
  kruskal_test(Frequency ~ Group) %>%
  # Adjust p-values for the 256 motifs
  adjust_pvalue(method = "BH") %>%
  # Filter for motifs that show a significant overall difference (FDR <= 0.05)
  filter(p.adj <= 0.05)

pairwise_results <- analysis_df %>%
  pivot_longer(cols = results_multi$Motif, 
               names_to = "Motif", 
               values_to = "Frequency") %>%
  group_by(Motif) %>%
  pairwise_wilcox_test(Frequency ~ Group, p.adjust.method = "BH")

# Multiple testing correction
results_wilcox <- results_wilcox %>%
  mutate(fdr = p.adjust(p_val, method = "BH")) %>%
  filter(fdr <= 0.05)


