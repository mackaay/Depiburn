setwd("./DEpiBurn_EMseq/")
library(dplyr)
library(GenomicRanges)
library(nnls)
library(pheatmap)
library(RColorBrewer)


library(Rsamtools)
library(GenomicAlignments)
library(tidyverse)

# 1. SETUP AND INITIALIZATION
SHORT_MIN <- 100
SHORT_MAX <- 150
LONG_MIN <- 151
LONG_MAX <- 220
DINUC_MIN <- 300
DINUC_MAX <- 340
BIN_SIZE <- 100000 # 100kb bin size

bamfile_list <- list.files(path = "/datasets/work/hb-burns-meth/work/Data/level2/ONT_5mC/", pattern = "_s.bam$", full.names = T, recursive = T)
sample_data <- data.frame(
  Sample_ID = gsub("_s.bam", "", basename(bamfile_list) ),
  BAM_File = c(bamfile_list),
  Group = c("Superficial", "Deep", "Deep", "Superficial", "Superficial", "Deep", "Superficial", "Deep", "Superficial", "Deep", "Superficial", "Full"),
  stringsAsFactors = FALSE
)
calculate_n50 <- function(lengths) {
  sorted_lengths <- sort(lengths, decreasing = TRUE)
  total_sum <- sum(as.numeric(sorted_lengths))
  cumsum_lengths <- cumsum(as.numeric(sorted_lengths))
  
  # Find the first read that pushes the sum over 50%
  n50 <- sorted_lengths[which(cumsum_lengths >= total_sum / 2)[1]]
  return(n50)
}

# 2. FRAGMENT LENGTH EXTRACTION AND FILTERING
extract_fragment_lengths <- function(bam_file) {
  #gal <- readGAlignmentPairs(bam_file, 
  #                           param = ScanBamParam(what = c("rname", "pos", "qwidth", "isize") ))
  gal <- scanBam(bam_file, 
                 param = ScanBamParam(what = c("rname", "pos", "qwidth", "isize", "mapq", "seq"), 
                                      flag = scanBamFlag(isProperPair = F, 
                                                         isSecondaryAlignment = FALSE, # Usually preferred for ONT
                                                         isSupplementaryAlignment = FALSE # Depending on your analysis
                                      ))
  )[[1]]

  df <- data.frame(
    chr = gal$rname,
    start = gal$pos,
    length = gal$qwidth, 
    mapq = gal$mapq,
    seq = gal$seq,
    stringsAsFactors = FALSE
  ) %>%
    filter(length > 0) # Only keep reads where an insert size was properly calculated
  return(df)
}

n50 <- NULL
for (i in 1:nrow(sample_data)) {
  print(paste("Processing:", sample_data$Sample_ID[i]))
  data <- extract_fragment_lengths(sample_data$BAM_File[i])
  data$Sample_ID <- sample_data$Sample_ID[i]
  data$Group <- sample_data$Group[i]
  n50_value <- calculate_n50(data$length)
  n50 <- c(n50, n50_value)
  print("Saving")
  save(data, file = paste("/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/data/ONTfragment_", sample_data$Sample_ID[i], ".rds", sep = "") )
}

# Combine all data into one large data frame
#all_fragments_df <- bind_rows(fragment_data_list)
all_fragments_df <- NULL
for (i in 1:nrow(sample_data)) {
  load(file = paste("/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/data/ONTfragment_", sample_data$Sample_ID[i], ".rds", sep = "")  )
  all_fragments_df <- rbind(all_fragments_df, data)
  print(paste("Processing:", sample_data$Sample_ID[i], " Done", sep = " "))
  save(all_fragments_df, file ="/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/data/ONTall_fragments_df.rds")
}


# 3. ANALYSIS I) SHORT-TO-LONG RATIO (FRAGMENTATION SCORE)
# Function to calculate GC content
calculate_gc_content <- function(sequence) {
  # Count G's and C's, and total sequence length
  gc_count <- str_count(sequence, "[GCgc]")
  total_length <- nchar(sequence)
  # GC bias correction is often done on the *sequence* of the cfDNA fragment itself
  return(gc_count / total_length)
}
# Calculate GC content for each fragment (This can be computationally intensive)
all_fragments_df <- all_fragments_df %>%
  mutate(gc_content = calculate_gc_content(seq)) %>%
  select(-seq) # Remove sequence data to save memory


library(ggplot2)
# 1. PREPARE THE DATA

for (i in 1:nrow(sample_data)) {
  print(sample_data$Sample_ID[i])
  load(file = paste("/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/data/ONTfragment_", sample_data$Sample_ID[i], ".rds", sep = "")  )
  # We filter to 0-500bp to focus on the most relevant mononucleosomal and dinucleosomal peaks.
plot_data <- data %>%
  filter(length > 0, length <= 1000)
# 2. CREATE THE DENSITY PLOT
# We use density instead of raw counts so that samples with different sequencing depths can be compared directly.

pdf(paste("./plots/109FragmentPlot_ONT_", sample_data$Sample_ID[i],".pdf", sep = ""), width = 8, height = 5)
p1 <- ggplot(plot_data, aes(x = length)) +
  geom_density(alpha = 0.2) +
  # Highlight the regions of interest from your methodology
  annotate("rect", xmin = 100, xmax = 150, ymin = 0, ymax = Inf, alpha = 0.1, fill = "blue") + # Short
  annotate("rect", xmin = 151, xmax = 220, ymin = 0, ymax = Inf, alpha = 0.1, fill = "green") + # Long
  annotate("rect", xmin = 300, xmax = 340, ymin = 0, ymax = Inf, alpha = 0.1, fill = "red") + # Dinucleosome
  # Formatting
  theme_minimal() +
  labs(
    title = sample_data$Sample_ID[i],
    subtitle = "Shaded areas: Blue (Short), Green (Long), Red (Dinucleosome)",
    x = "Fragment Length (bp)",
    y = "Density"
    #color = "Patient Status",
    #fill = "Patient Status"
  ) +
  scale_x_continuous(breaks = seq(0, 1000, by = 50)) +
  theme(legend.position = "top")
print(p1)
dev.off()
}




# Define fragment length boundaries (based on the provided text)
SHORT_MIN <- 100
SHORT_MAX <- 150
LONG_MIN <- 151
LONG_MAX <- 220
DINUC_MIN <- 300
DINUC_MAX <- 340
BIN_SIZE <- 100000 # 100kb bin size

calculate_gc_content <- function(sequence) {
  # Count G's and C's, and total sequence length
  gc_count <- str_count(sequence, "[GCgc]")
  total_length <- nchar(sequence)
  # GC bias correction is often done on the *sequence* of the cfDNA fragment itself
  return(gc_count / total_length)
}

data <- data %>%
  mutate(gc_content = calculate_gc_content(seq)) %>%
  select(-seq) # Remove sequence data to save memory

fragmentation_scores_df <- data %>%
  # Filter fragments only relevant for this analysis (100-220 bp)
  filter(length >= SHORT_MIN, length <= LONG_MAX) %>%
  
  # Create the 100kb non-overlapping bins
  mutate(
    # Ensure chr is factored correctly (Rsamtools usually handles this)
    chr = as.character(chr), 
    bin_start = floor(start / BIN_SIZE) * BIN_SIZE,
    bin_id = paste(chr, bin_start, sep = "_")
  ) %>%
  
  # Classify fragments as Short or Long
  mutate(
    type = case_when(
      length >= SHORT_MIN & length <= SHORT_MAX ~ "Short",
      length >= LONG_MIN & length <= LONG_MAX ~ "Long",
      TRUE ~ NA_character_ # Should not happen after filtering, but good practice
    )
  ) %>%
  filter(!is.na(type))
# Calculate counts per bin, per sample
bin_counts <- fragmentation_scores_df %>%
  group_by(Sample_ID,  bin_id, type, gc_content) %>%
  summarise(
    count = n(),
    .groups = "drop_last"
  )
# Reshape and calculate the raw ratio
raw_ratios <- bin_counts %>%
  pivot_wider(names_from = type, values_from = count, values_fill = list(count = 0)) %>%
  mutate(
    # Add a small pseudo-count to avoid division by zero
    Ratio = (Short + 1) / (Long + 1) 
  ) %>%
  ungroup()
# GC Bias Correction (LOESS)
# This is a critical step. We model the raw ratio as a function of GC content
corrected_ratios <- raw_ratios %>%
  group_by(Sample_ID) %>%
  # The LOESS model is calculated per sample, using the raw Ratio and GC content
  mutate(
    loess_fit = list(loess(Ratio ~ gc_content, data = cur_data())),
    # Calculate the LOESS-smoothed value (the expected ratio for that GC content)
    Expected_Ratio = predict(loess_fit[[1]], newdata = cur_data()),
    # The corrected ratio is the raw ratio divided by the expected ratio
    # Alternatively, residuals (Ratio - Expected_Ratio) can be used.
    # We use division to get a correction factor.
    Corrected_Ratio = Ratio / Expected_Ratio
  ) %>%
  ungroup()
# Calculate the mean fragmentation score for each sample
mean_fragmentation_scores <- corrected_ratios %>%
  group_by(Sample_ID) %>%
  summarise(
    Mean_Fragmentation_Score = mean(Corrected_Ratio, na.rm = TRUE),
    .groups = "drop"
  )
# Display the results 
print(mean_fragmentation_scores)

