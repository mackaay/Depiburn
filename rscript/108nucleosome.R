library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(dplyr)

load(file = "/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/data/TSS_tissue.RData")
setwd("/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/")

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


# Calculate Windowed Protection Score (WPS)###########
#' @param bam_file Path to the BAM file
#' @param tss_gr GenomicRanges object of TSS regions (centered, e.g., 5kb)
#' @param window_size Integer, usually 120bp
calc_wps <- function(bam_file, tss_gr, window_size = 120) {
  # 1. Filter Mononucleosomal Fragments (120 - 200 bp)
  param <- ScanBamParam(
    flag = scanBamFlag(isProperPair = TRUE, isDuplicate = FALSE),
    what = c("pos", "isize", "qwidth")
  )
  # Read alignments (simplifying to mononucleosomal range)
  # Note: isize is the template length (fragment size)
  raw_reads <- readGAlignmentPairs(bam_file, param = param)
  reads_gr <- granges(raw_reads)
  print("Reading bam")
  # Filter by length: 120 <= length <= 200
  reads_filtered <- reads_gr[width(reads_gr) >= 120 & width(reads_gr) <= 200]
  print("Filtering Mononucleosome")
  # 2. Sliding Window Calculation
  # For each TSS region, calculate WPS
  results <- lapply(seq_along(tss_gr), function(i) {
    region <- tss_gr[i]

    # Get fragments overlapping this specific TSS region
    local_reads <- subsetByOverlaps(reads_filtered, region)
    
    # Define sliding windows (1bp steps)
    starts <- seq(start(region), end(region) - window_size, by = 1)
    ends <- starts + window_size - 1
    windows <- GRanges(seqnames(region), IRanges(starts, ends))
    
    # 3. Compute Overlaps (Fully Spanning vs. Partially Spanning)
    # WPS = (Full Spanning) - (Non-Full Spanning)
    
    # We use findOverlaps to identify spatial relationships
    hits <- findOverlaps(windows, local_reads)
    
    # Convert hits to a logical matrix or list to check boundaries
    # A fragment fully spans if: read_start <= window_start AND read_end >= window_end
    read_starts <- start(local_reads)[subjectHits(hits)]
    read_ends <- end(local_reads)[subjectHits(hits)]
    win_starts <- start(windows)[queryHits(hits)]
    win_ends <- end(windows)[queryHits(hits)]
    
    is_spanning <- (read_starts <= win_starts) & (read_ends >= win_ends)
    
    # Aggregate counts per window
    #spanning_counts <- aggregate(is_spanning, list(queryHits(hits)), sum)
    #non_spanning_counts <- aggregate(!is_spanning, list(queryHits(hits)), sum)
    res_df <- data.frame(win_idx = queryHits(hits), is_spanning = is_spanning) %>%
      group_by(win_idx) %>%
      summarize(
        spanning = sum(is_spanning),
        non_spanning = sum(!is_spanning),
        .groups = "drop"
      )
    # Calculate raw WPS
    #wps_raw <- rep(0, length(windows))
    #wps_raw[spanning_counts$Group.1] <- spanning_counts$x - non_spanning_counts$x
    wps_raw <- rep(0, length(windows))
    wps_raw[res_df$win_idx] <- res_df$spanning - res_df$non_spanning
    
    # 4. Adjustment and Smoothing
    # Subtract mean WPS
    wps_adj <- wps_raw - mean(wps_raw)
    
    # Smooth with Gaussian kernel (bandwidth = 30bp)
    # stats::ksmooth is efficient for this
    wps_smooth <- ksmooth(x = 1:length(wps_adj), y = wps_adj, 
                          kernel = "normal", bandwidth = 30)$y
    
    return(wps_smooth)
  })
  
  return(results)
}

calc_wps_fast <- function(bam_file, tss_gr, window_size = 120) {
  # 1. Read and Filter BAM (One-time operation)
  param <- ScanBamParam(
    flag = scanBamFlag(isProperPair = TRUE, isDuplicate = FALSE),
    what = c("pos", "isize")
  )
  raw_reads <- readGAlignmentPairs(bam_file, param = param)
  reads_gr <- granges(raw_reads)
  
  # Filter for mononucleosomal fragments (120-200bp)
  reads_filtered <- reads_gr[width(reads_gr) >= 120 & width(reads_gr) <= 200]
  
  # 2. Vectorized Counting using Coverage
  # 'cov_all' counts fragments covering each base
  cov_all <- coverage(reads_filtered)
  
  # To find "Fully Spanning", we look at fragments where (End - Start) >= window_size
  # and the fragment covers the entire window.
  # A trick: Fragments spanning window [i, i+k] are those that cover 
  # both the start AND the end of the window AND are long enough.
  
  # Get coverage of fragment STARTS and ENDS
  cov_starts <- coverage(resize(reads_filtered, width = 1, fix = "start"))
  cov_ends <- coverage(resize(reads_filtered, width = 1, fix = "end"))
  
  # 3. Process each TSS region efficiently
  # We use 'Views' to extract data from the Rle coverage objects instantly
  results <- lapply(seq_along(tss_gr), function(i) {
    region <- tss_gr[i]
    chr <- as.character(seqnames(region))
    # Extract coverage for the specific chromosome/region
    # view_cov is the number of fragments covering each base in the TSS region
    view_cov <- as.numeric(cov_all[[chr]][ranges(region)])
    
    # Calculate WPS using a rolling window approach
    # We use a moving sum/mean or just indexing
    # n_spanning = fragments covering every base in window
    # Because we filtered fragments to be > 120bp, we can approximate 
    # spanning as the minimum coverage within the window.
    # Rollapply equivalent using RcppRoll for extreme speed
    # Or simple zoo::rollapply
    library(RcppRoll)
    
    # WPS formula simplified: 
    # (Count of fragments covering entire window) - (Count of fragments covering part of window)
    # count_any: fragments overlapping the window at all
    # count_full: fragments where start <= win_start AND end >= win_end
    
    # Optimized approximation:
    win_starts <- 1:(length(view_cov) - window_size + 1)
    win_ends <- win_starts + window_size - 1
    # For a fragment to span, it must cover both start and end of window
    # This is a common shortcut for WPS:
    cov_vec <- view_cov[win_starts : win_ends] # This needs matrix handling or roll_min
    # Using RcppRoll to get min coverage in window (Spanning) 
    # and max coverage (Any overlap)
    spanning <- roll_min(view_cov, n = window_size)
    any_overlap <- roll_max(view_cov, n = window_size)
    non_spanning <- any_overlap - spanning
    wps_raw <- spanning - non_spanning
    # 4. Adjustment and Smoothing
    wps_adj <- wps_raw - mean(wps_raw)
    wps_smooth <- ksmooth(x = 1:length(wps_adj), y = wps_adj, 
                          kernel = "normal", bandwidth = 30)$y
    return(wps_smooth)
  })
  return(results)
}

library(tidyr)
# Assume 'active_wps_list' and 'inactive_wps_list' are results from calc_wps_ont()
# These are lists of numeric vectors (each vector is the 5kb signal for one gene)
prepare_profile <- function(wps_list, label) {
  # Calculate the mean signal across all genes at each position
  # Each entry in the list must be the same length (e.g., ~4880 bps after windowing)
  mat <- do.call(rbind, wps_list)
  mean_signal <- colMeans(mat, na.rm = TRUE)
  
  data.frame(
    position = seq(-2440, 2440), # Adjusted for 5kb - 120bp window
    WPS = mean_signal,
    Group = label
  )
}


for (i in 1:nrow(sample_data)) {
  print(sample_data$BAM_File[i])
  SA_footprint <- calc_wps(sample_data$BAM_File[i], SA_tss)
blood_footprint <- calc_wps(sample_data$BAM_File[i], blood_tss)
nerve_footprint <- calc_wps(sample_data$BAM_File[i], nerve_tss)
skin_footprint <- calc_wps(sample_data$BAM_File[i], skin_tss)
#inactive_footprint <- calc_wps(sample_data$BAM_File[i], inactive_tss)
print("tissue footprint done")

plot_data <- rbind(
  prepare_profile(SA_footprint, "SA"),
  prepare_profile(blood_footprint, "Blood"), 
  prepare_profile(nerve_footprint, "Nerve"),
  prepare_profile(skin_footprint, "Skin")
  #prepare_profile(inactive_footprint, "inactive")
)
  plot_data$Sample <- basename(sample_data$BAM_File[i])
  print("Done")
save(plot_data, file =  paste("/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/data/nucleosome_mono_", basename(sample_data$BAM_File[i]), ".rds", sep = ""))
}



#Plot Nucleosome analysis ####
plot_data_merge <- NULL
for (i in  c(1:22)) {
  load(file = paste("/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/data/nucleosome_mono_", basename(sample_data$BAM_File[i]), ".rds", sep = ""))
  plot_data$burn_group <- sample_data$Group[i]
  plot_data_merge <- rbind(plot_data_merge, plot_data)
  print(basename( sample_data$BAM_File[i]))
}
plot_data_merge$burn_group <- factor(plot_data_merge$burn_group, levels = c("Superficial", "Deep", "Full"))

pal.depth <- c("Superficial"= "#999999","Deep"= "#E69F00","Full"= "#56B4E9")
library(ggplot2)
ggplot(plot_data_merge, aes(x = position, y = WPS, color = burn_group, fill = burn_group)) +
  stat_summary(geom = "ribbon", 
               fun.data = mean_cl_boot, 
               alpha = 0.1, 
               color = NA) +
  stat_summary(geom = "line", 
               fun = mean, 
               size = 1) +
  scale_color_manual(values = pal.depth) +
  scale_fill_manual(values = pal.depth) +         # changing the y axis nber format
  theme_minimal(base_size = 15  )+
  labs(
    title = "Nucleosome Footprinting (WPS) at TSS",
    x = "Distance from TSS (bp)",
    y = "Mean Adjusted WPS \n(Aggregate)",
    subtitle = "Mean signal with 95% CI"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")

ggplot(plot_data_merge, aes(x = position, y = WPS, color = burn_group, fill = burn_group)) +
  stat_summary(geom = "ribbon", 
               fun.data = mean_cl_boot, 
               alpha = 0.1, 
               color = NA) +
  stat_summary(geom = "line", 
               fun = mean, 
               size = 1) +
  scale_color_manual(values = pal.depth) +
  scale_fill_manual(values = pal.depth) +         # changing the y axis nber format
  theme_minimal(base_size = 15  )+
  labs(
    title = "Nucleosome Footprinting (WPS) at TSS",
    x = "Distance from TSS (bp)",
    y = "Mean Adjusted WPS \n(Aggregate)",
    subtitle = "Mean signal with 95% CI"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + xlim(-1200, 1200)

ggplot(plot_data_merge, aes(x = position, y = WPS, color = burn_group, fill = burn_group)) +
  stat_summary(geom = "ribbon", 
               fun.data = mean_cl_boot, 
               alpha = 0.1, 
               color = NA) +
  stat_summary(geom = "line", 
               fun = mean, 
               size = 0.8) +
  scale_color_manual(values = pal.depth) +
  facet_wrap(~Group, nrow = 1) +
  scale_fill_manual(values = pal.depth) +         # changing the y axis nber format
  theme_minimal(base_size = 15  )+
  labs(
    title = "Nucleosome Footprinting (WPS) at TSS",
    x = "Distance from TSS (bp)",
    y = "Mean Adjusted WPS \n(Aggregate)",
    subtitle = "Mean signal with 95% CI"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") 


