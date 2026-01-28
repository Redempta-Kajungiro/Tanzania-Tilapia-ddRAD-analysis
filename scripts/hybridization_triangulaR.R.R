# Set working directory
setwd("/home/redempta/work_2025/tilapia/new_analysis")

# Install and load packages
devtools::install_github("omys-omics/triangulaR")
install.packages("ggplot2")
library(triangulaR)
library(ggplot2)
library(vcfR)

# Read VCF and population map
data <- read.vcfR("populations.snps.vcf", verbose = FALSE)

popmap <- read.table(
  file = "populations_info_lines_.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

colnames(popmap) <- c("id", "pop")

# Subset populations

keep_pops <- c("Kunduchi", "Ruhila", "Wami", "Mindu",
               "Karanga", "Pangani_Nile", "Utete",
               "Victoria", "Nyamisati")

keep_samples <- popmap$id[popmap$pop %in% keep_pops]
popmap_sub <- popmap[popmap$pop %in% keep_pops, ]

# Subset VCF
data_sub <- data[, c("FORMAT", keep_samples)]

# Allele frequency difference filtering
vcfR.diff <- alleleFreqDiff(
  vcfR = data_sub,
  pm = popmap_sub,
  p1 = "Victoria",
  p2 = "Nyamisati",
  difference = 1
)

# Calculate hybrid index and heterozygosity
hi.het <- hybridIndex(
  vcfR = vcfR.diff,
  pm = popmap_sub,
  p1 = "Victoria",
  p2 = "Nyamisati"
)

# Define colors for populations
cols9 <- c(
  Karanga       = "#E69F00",  
  Kunduchi      = "#7B3294",  
  Mindu         = "#F781BF",  
  Nyamisati     = "#1B9E77",  
  Pangani_Nile  = "#F0E442",  
  Ruhila        = "#377EB8",  
  Utete         = "#D95F02",  
  Victoria      = "#B2182B",  
  Wami          = "#999999"   
)

# Generate triangle plot
p1 <- triangle.plot(
  hi.het,
  col = adjustcolor(cols9, alpha.f = 0.6),
  cex = 1
)
p1
ggsave("tri.tiff", plot = p1, width = 10, height = 6, dpi = 300)

# Combine multiple plots
# Assuming p5, p75, p09 are defined similarly to p1
combined <- plot_grid(
  p5, p75, p09, p1,
  labels = c("A", "B", "C", "D"),
  label_size = 10,
  nrow = 2
)

# Save combined plot
ggsave("combinedtrianew.png", combined, width = 10, height = 6, dpi = 300)
