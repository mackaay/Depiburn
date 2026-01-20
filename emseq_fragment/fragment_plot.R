library(Rsamtools)
library(GenomicAlignments)
library(tidyverse)



load(file ="/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/data/fragment_QUT187.rds")
print("Reading Done")
pdf("/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/plots/108FragmentPlot_QUT187.pdf", width = 8, height=5)
ggplot(data, aes(x = length)) +
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
    y = "Density"
    #color = "Patient Status",
    #fill = "Patient Status"
  ) +
  scale_x_continuous(breaks = seq(0, 1000, by = 50)) +
  theme(legend.position = "top")
dev.off()
print("Plot Done")

load(file ="/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/data/fragment_QUT197.rds")
print("Reading Done")
pdf("/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/plots/108FragmentPlot_QUT197.pdf", width = 8, height=5)
ggplot(data, aes(x = length)) +
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
    y = "Density"
    #color = "Patient Status",
    #fill = "Patient Status"
  ) +
  scale_x_continuous(breaks = seq(0, 1000, by = 50)) +
  theme(legend.position = "top")
dev.off()
print("Plot Done")

load(file ="/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/data/fragment_QUT198.rds")
print("Reading Done")
pdf("/datasets/work/hb-diab-cfdna/work/Users/chenkai/Meth_decon/DEpiBurn_EMseq/plots/108FragmentPlot_QUT198.pdf", width = 8, height=5)
ggplot(data, aes(x = length)) +
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
    y = "Density"
    #color = "Patient Status",
    #fill = "Patient Status"
  ) +
  scale_x_continuous(breaks = seq(0, 1000, by = 50)) +
  theme(legend.position = "top")
dev.off()
print("Plot Done")
