#setwd("./DEpiBurn_EMseq/")
#this R code is designed for HPC running. 

library(Rsamtools)
library(GenomicAlignments)
library(tidyverse)

# ----------------------------------------------------------------------
# 1. SETUP AND INITIALIZATION
# ----------------------------------------------------------------------
# Define fragment length boundaries (based on the provided text)
SHORT_MIN <- 100
SHORT_MAX <- 150
LONG_MIN <- 151
LONG_MAX <- 220
DINUC_MIN <- 300
DINUC_MAX <- 340
BIN_SIZE <- 100000 # 100kb bin size

bamfile_list <- list.files(path = "/datasets/work/hb-burns-meth/work/Data/level2/", pattern = "_sd.bam$", full.names = T, recursive = T)
sample_data <- data.frame(
  Sample_ID = substr( basename(bamfile_list), start = 24L, stop = 29L ),
  BAM_File = c(bamfile_list),
  Group = c("Full", "Deep", "Deep", "Deep", "Deep", "Full", "Superficial", "Superficial", #207
            "Superficial", "Superficial", "Deep", "Full", 
            "Superficial", "Deep", "Superficial", "Superficial", "Deep", "Superficial", 
            "Deep", "Superficial", "Superficial", "Full"),
  stringsAsFactors = FALSE
)

# ----------------------------------------------------------------------
# 2. FRAGMENT LENGTH EXTRACTION AND FILTERING
# ----------------------------------------------------------------------
extract_fragment_lengths <- function(bam_file) {
  #gal <- readGAlignmentPairs(bam_file, 
  #                           param = ScanBamParam(what = c("rname", "pos", "qwidth", "isize") ))
  gal <- scanBam(bam_file, 
                 param = ScanBamParam(what = c("rname", "pos", "qwidth", "isize", "mapq", "seq"), 
                                      flag = scanBamFlag(isProperPair = TRUE, isUnmappedQuery = FALSE, hasUnmappedMate = FALSE))
                 )[[1]]
  df <- data.frame(
    chr = gal$rname,
    start = gal$pos,
    length = abs(gal$isize), # Insert size is what we want (always use absolute value)
    mapq = gal$mapq,
    seq = gal$seq,
    stringsAsFactors = FALSE
  ) %>%
    filter(length > 0) # Only keep reads where an insert size was properly calculated
  return(df)
}
# Apply the function to all BAM files and store the results
#sample_data <- data.frame(Sample_ID= c("test1", "test2"), BAM_File = c("./DEpiBurn_EMseq/data/Skin_BlisterFluid_Burn_QUT174_sd.chr22.bam","./DEpiBurn_EMseq/data/Skin_BlisterFluid_Burn_QUT174_sd.chr22.bam"), 
#                          Group = c("Full", "Full"))
#fragment_data_list <- lapply(1:nrow(sample_data), function(i) {
#  message(paste("Processing:", sample_data$Sample_ID[i]))
#  data <- extract_fragment_lengths(sample_data$BAM_File[i])
#  data$Sample_ID <- sample_data$Sample_ID[i]
#  data$Group <- sample_data$Group[i]
#  message(paste(sample_data$Sample_ID[i]), "Done")
#  return(data)
#})
for (i in 23:nrow(sample_data)) {
  print(paste("Processing:", sample_data$Sample_ID[i]))
  data <- extract_fragment_lengths(sample_data$BAM_File[i])
  data$Sample_ID <- sample_data$Sample_ID[i]
  data$Group <- sample_data$Group[i]
  print("Saving")
  save(data, file = paste("/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/data/fragment_", sample_data$Sample_ID[i], ".rds", sep = "") )
}

# Combine all data into one large data frame
#all_fragments_df <- bind_rows(fragment_data_list)
all_fragments_df <- NULL
for (i in 1:nrow(sample_data)) {
  load(file = paste("/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/data/fragment_", sample_data$Sample_ID[i], ".rds", sep = "")  )
  all_fragments_df <- rbind(all_fragments_df, data)
  print(paste("Processing:", sample_data$Sample_ID[i], " Done", sep = " "))
  save(all_fragments_df, file ="/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/data/all_fragments_df.rds")
}

save(all_fragments_df, file ="/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/data/all_fragments_df.rds")



# ----------------------------------------------------------------------
# 3. ANALYSIS I) SHORT-TO-LONG RATIO (FRAGMENTATION SCORE)
# ----------------------------------------------------------------------
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
# We filter to 0-500bp to focus on the most relevant mononucleosomal and dinucleosomal peaks.
plot_data <- all_fragments_df %>%
  filter(length > 0, length <= 500)
# 2. CREATE THE DENSITY PLOT
# We use density instead of raw counts so that samples with different sequencing depths can be compared directly.
ggplot(plot_data, aes(x = length, color = Group, fill = Group)) +
  geom_density(alpha = 0.2) +
  # Highlight the regions of interest from your methodology
  annotate("rect", xmin = 100, xmax = 150, ymin = 0, ymax = Inf, alpha = 0.1, fill = "blue") + # Short
  annotate("rect", xmin = 151, xmax = 220, ymin = 0, ymax = Inf, alpha = 0.1, fill = "green") + # Long
  annotate("rect", xmin = 300, xmax = 340, ymin = 0, ymax = Inf, alpha = 0.1, fill = "red") + # Dinucleosome
  # Formatting
  theme_minimal() +
  labs(
    title = "cfDNA Fragment Size Distribution",
    subtitle = "Shaded areas: Blue (Short), Green (Long), Red (Dinucleosome)",
    x = "Fragment Length (bp)",
    y = "Density",
    color = "Patient Status",
    fill = "Patient Status"
  ) +
  scale_x_continuous(breaks = seq(0, 500, by = 50)) +
  theme(legend.position = "top")



fragmentation_scores_df <- all_fragments_df %>%
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
  group_by(Sample_ID, Group, bin_id, type, gc_content) %>%
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
  group_by(Sample_ID, Group) %>%
  summarise(
    Mean_Fragmentation_Score = mean(Corrected_Ratio, na.rm = TRUE),
    .groups = "drop"
  )
# Display the results 
print(mean_fragmentation_scores)


# Perform non-parametric Wilcoxon test (Mann-Whitney U test)
wilcox_test_i <- wilcox.test(
  Mean_Fragmentation_Score ~ Group,
  data = mean_fragmentation_scores
)
# Output the test result
message("\n--- Fragmentation Score (i) Wilcoxon Test Result ---")
print(wilcox_test_i)


# ----------------------------------------------------------------------
# 4. ANALYSIS II) DINUCLEOSOME FRACTION
# ----------------------------------------------------------------------
# Calculate the total number of fragments per sample
total_reads <- all_fragments_df %>%
  group_by(Sample_ID) %>%
  summarise(Total_Fragments = n(), .groups = "drop")
# Calculate the number of dinucleosome fragments (300-340 bp)
dinucleosome_counts <- all_fragments_df %>%
  filter(length >= DINUC_MIN, length <= DINUC_MAX) %>%
  group_by(Sample_ID) %>%
  summarise(Dinucleosome_Fragments = n(), .groups = "drop")
# Combine and calculate the fraction
dinucleosome_fractions <- total_reads %>%
  left_join(dinucleosome_counts, by = "Sample_ID") %>%
  replace_na(list(Dinucleosome_Fragments = 0)) %>%
  left_join(sample_data %>% select(Sample_ID, Group), by = "Sample_ID") %>%
  mutate(
    Fraction_Dinucleosome = Dinucleosome_Fragments / Total_Fragments
  )
# Display the results 
print(dinucleosome_fractions)
# Statistical Test ii)
wilcox_test_ii <- wilcox.test(
  Fraction_Dinucleosome ~ Group,
  data = dinucleosome_fractions
)
# Output the test result
message("\n--- Dinucleosome Fraction (ii) Wilcoxon Test Result ---")
print(wilcox_test_ii)


