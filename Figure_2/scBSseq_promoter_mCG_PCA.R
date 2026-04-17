###--- run in bash ---###
#!/bin/bash
##-- activate methscan
source /g/data/vn68/aa3618/tools/methscan/methscan_env/bin/activate

##-- prepare data
INPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/bismark/meth_bedGraph/"
OUTPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/bismark/meth_bedGraph/methscan/"
mkdir -p ${OUTPUT_DIR}

methscan prepare ${INPUT_DIR}*.cov ${OUTPUT_DIR}

##-- filter low quality cells
INPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/bismark/meth_bedGraph/methscan/"
OUTPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/bismark/meth_bedGraph/methscan/filter/"
mkdir -p ${OUTPUT_DIR}

methscan filter --min-sites 1000000 --min-meth 10 --max-meth 70 ${INPUT_DIR} ${OUTPUT_DIR}

##-- run methscan smooth 
# this command treats all cells as a pseudo-bulk sample and calculates smoothed mean methylation along the genome
methscan smooth ${OUTPUT_DIR}

##-- intersect matrix with promoters
REGION_DIR="/g/data/vn68/aa3618/genomes/Sminthopsis_crassicaudata/T2T_ONT_genome/gtf/"
INPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/bismark/meth_bedGraph/methscan/filter/"
OUTPUT_DIR="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/bismark/meth_bedGraph/methscan/promoter_matrix/"
mkdir -p ${OUTPUT_DIR}

methscan matrix --threads 4 ${REGION_DIR}dunnart_male_T2T_promoter_2kb.bed ${INPUT_DIR} ${OUTPUT_DIR}

###--- run in R ---###
library(tidyverse)
library(irlba)
library(data.table)
library(rtracklayer)
library(plyr)
library(ggplot2)
library(circlize)
library(GenomicRanges)

##-- get PCA of promoter mCG
# read in promoter methylation matrix
meth_mtx <- read.csv("~/Desktop/dunnart_reprogramming/scBSseq/methscan/mean_shrunken_residuals.csv.gz", row.names=1) %>%
  as.matrix()

# PCA that iteratively imputes missing values
prcomp_iterative <- function(x, n=10, n_iter=50, min_gain=0.001, ...) {
  mse <- rep(NA, n_iter)
  na_loc <- is.na(x)
  x[na_loc] = 0  # zero is our first guess
  
  for (i in 1:n_iter) {
    prev_imp <- x[na_loc]  # what we imputed in the previous round
    # PCA on the imputed matrix
    pr <- prcomp_irlba(x, center = F, scale. = F, n = n, ...)
    # impute missing values with PCA
    new_imp <- (pr$x %*% t(pr$rotation))[na_loc]
    x[na_loc] <- new_imp
    # compare our new imputed values to the ones from the previous round
    mse[i] = mean((prev_imp - new_imp) ^ 2)
    # if the values didn't change a lot, terminate the iteration
    gain <- mse[i] / max(mse, na.rm = T)
    if (gain < min_gain) {
      message(paste(c("\n\nTerminated after ", i, " iterations.")))
      break
    }
  }
  pr$mse_iter <- mse[1:i]
  pr
}

# run modified PCA on the centered methylation matrix
pca <- meth_mtx %>%
  scale(center = T, scale = F) %>%
  prcomp_iterative(n = 15)  # increase this value to e.g. 15 for real data sets

# calculate percentage of variance explained
variance_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100

# create labels with variance explained
pc1_label <- paste0("PC1 (", round(variance_explained[1], 1), "%)")
pc2_label <- paste0("PC2 (", round(variance_explained[2], 1), "%)")

# plot PCA
pdf("PCA_promoter_mCG_methscan.pdf", useDingbats = FALSE, width = 4, height = 4)
pca_tbl <- as_tibble(pca$x) %>% 
  add_column(cell=rownames(meth_mtx))

print(
  ggplot(pca_tbl, aes(x = PC1, y = PC2)) +
    geom_point(size = 3, 
               color = "black",
               stroke = 0.8,
               shape = 21,
               fill = "#2E86AB") +
    coord_fixed() +
    labs(title = "PCA of promoter mCG - E9 dunnart",
         x = pc1_label,
         y = pc2_label) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, 
                                margin = margin(b = 5)),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40",
                                   margin = margin(b = 15)),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "gray90", size = 0.5),
      panel.grid.minor = element_line(color = "gray95", size = 0.3),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "gray80", fill = NA, size = 0.8)
    ) +
    guides(color = "none", fill = "none") +
    scale_x_continuous(limits = c(-20, 30), expand = c(0, 0))
) 

dev.off()

##-- plot genome mCG
genome_mCG <- read.csv("average_methylation_per_sample.csv")
colnames(genome_mCG)[1] <- "sample"
colnames(genome_mCG)[2] <- "genome_methylation"

# extract sample name
pca_tbl <- pca_tbl %>%
  mutate(sample = str_replace(cell, ".*-(P\\d+)-([A-H]\\d{2})\\.bismark", "\\1.\\2"))

print(sum(pca_tbl$sample %in% genome_mCG$sample))  # Should be > 0

# join with genome_mCG
pca_merged <- pca_tbl %>%
  left_join(genome_mCG, by = "sample")

# plot PCA coloured by genome methylation
pdf("PCA_promoter_mCG_coloured_by_genome.pdf", useDingbats = FALSE, width = 8, height = 6)

print(
  ggplot(pca_merged, aes(x = PC1, y = PC2)) +
    geom_point(aes(fill = genome_methylation), 
               shape = 21,           # Filled circle with border
               color = "black",      # Border color
               size = 3, 
               stroke = 0.8) +       # Border thickness
    scale_fill_gradientn(
      colours = c("#00008b", "white", "#D20A2E"),
      limits = c(0, 1),
      na.value = "grey"
    ) +
    coord_fixed() +
    labs(
      title = "PCA of promoter mCG - E9 dunnart",
      x = pc1_label,
      y = pc2_label,
      fill = "genome mCG"  # Note: matches 'fill' aesthetic
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, 
                                margin = margin(b = 5)),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40",
                                   margin = margin(b = 15)),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "gray90", size = 0.5),
      panel.grid.minor = element_line(color = "gray95", size = 0.3),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "gray80", fill = NA, size = 0.8)
    ) +
    scale_x_continuous(limits = c(-20, 30), expand = c(0, 0))
)
dev.off()

##-- plot elf5 mCG
# define functions
overlapSums <- function(x, y, z, na.rm=FALSE) {
  stopifnot(class(x)=="GRanges")
  stopifnot(class(y)=="GRanges")
  stopifnot(is.numeric(z) || is.integer(z))
  stopifnot(length(x)==length(z))
  
  ov <- as.matrix(findOverlaps(y, x))
  ovSums <- numeric(length(y))
  ovSums[] <- NA_real_
  ovSums[unique(ov[,1])] <- viewSums(Views(z[ov[,2]], ranges(Rle(ov[,1]))), na.rm=na.rm)
  ovSums
}

overlapRatios <- function(x, y, C, cov, minCov, na.rm=FALSE) {
  C.sum <- overlapSums(x, y, values(x)[[C]], na.rm=na.rm)
  cov.sum <- overlapSums(x, y, values(x)[[cov]], na.rm=na.rm)
  
  ov_all <- as.data.frame(findOverlaps(y, x)) 
  number_of_CpGs <- rle(ov_all[,1])$lengths
  
  CpG_ov_regions <- x[ov_all$subjectHits]
  ov_all$cov <- as.numeric(values(CpG_ov_regions)[[cov]])
  ov_all$cov_binary <- as.numeric(ov_all$cov >= minCov)
  number_of_CpGs_minCov <- as.numeric(tapply(ov_all$cov_binary, ov_all$queryHits, sum))
  
  CpGs_covered <- number_of_CpGs_minCov / number_of_CpGs
  rat <- as.data.frame(C.sum / cov.sum)
  rat_cov <- cbind(rat, CpGs_covered)
  
  colnames(rat_cov) <- c("rat", "CpGs_covered")
  rat_cov$rat[rat_cov$CpGs_covered < 0.5] <- NA
  
  rat_cov[,1]
}

# load bigtable
load("dunnart_E9_scBSseq_bigtable.RData")

# convert to GRanges object and extract relevant columns
CpG <- GRanges(bigtable, ignore.mcols=TRUE)
values(CpG) <- as.data.frame(bigtable)[,-c(1:5)]  # keep only methylation columns

# load ELF5 promoter region from bed file
# chr6  36949684  36951684
elf5_promoter <- import.bed("elf5.bed")

# identify column pairs
col <- seq(1, ncol(mcols(CpG)), by=2)
sample_names <- sub("\\.(C|cov)$", "", colnames(mcols(CpG))[col])
col <- setNames(col, sample_names)

# set minimum coverage threshold (requiring at least 1 read per CpG)
minCov <- 1  

# calculate ELF5 promoter methylation per sample
elf5_methylation <- numeric(length(col))
names(elf5_methylation) <- sample_names

for(i in seq_along(col)) {
  C_col <- col[i]
  cov_col <- C_col + 1
  
  # use overlapRatios function to get methylation ratio for ELF5 promoter
  elf5_meth <- overlapRatios(
    x = CpG, 
    y = elf5_promoter, 
    C = colnames(mcols(CpG))[C_col], 
    cov = colnames(mcols(CpG))[cov_col], 
    minCov = minCov, 
    na.rm = TRUE
  )
  
  elf5_methylation[i] <- elf5_meth[1]  # First (and only) element
  
  if(i %% 20 == 0) {
    cat("Processed", i, "of", length(col), "samples\n")
  }
}

# create results data frame for ELF5 promoter methylation
elf5_results <- data.frame(
  sample = names(elf5_methylation),
  elf5_methylation = elf5_methylation,
  stringsAsFactors = FALSE
)

# add plate and row information from sample names
elf5_results$Plate <- sub("\\.(A|B|C|D).*", "", elf5_results$sample)
elf5_results$Row <- sub(".*\\.", "", sub("\\.[0-9]+.*", "", elf5_results$sample))

# save results
print(head(elf5_results))
write.csv(elf5_results, "ELF5_promoter_methylation_per_sample.csv", row.names = FALSE)
saveRDS(elf5_results, "ELF5_promoter_methylation_per_sample.rds")

elf5_mCG <- read.csv("ELF5_promoter_methylation_per_sample.csv")
# extract sample name
pca_tbl <- pca_tbl %>%
  mutate(sample = str_replace(cell, ".*-(P\\d+)-([A-H]\\d{2})\\.bismark", "\\1.\\2"))
print(sum(pca_tbl$sample %in% elf5_mCG$sample))  # Should be > 0
# join with elf5_mCG
pca_merged <- pca_tbl %>%
  left_join(elf5_mCG, by = "sample")

# plot PCA coloured by elf5 methylation
pdf("PCA_promoter_mCG_coloured_by_elf5.pdf", useDingbats = FALSE, width = 8, height = 6)
print(
  ggplot(pca_merged, aes(x = PC1, y = PC2)) +
    geom_point(data = pca_merged %>% filter(is.na(elf5_methylation)),
               shape = 21,           
               color = "black",      
               fill = "grey",
               size = 3, 
               stroke = 0.8) +       
    geom_point(data = pca_merged %>% filter(!is.na(elf5_methylation)),
               aes(fill = elf5_methylation), 
               shape = 21,           
               color = "black",      
               size = 3, 
               stroke = 0.8) +       
    scale_fill_gradientn(
      colours = c("#00008b", "white", "#D20A2E"),
      limits = c(0, 1),
      na.value = "grey"
    ) +
    coord_fixed() +
    labs(
      title = "PCA of promoter mCG - E9 dunnart",
      x = pc1_label,
      y = pc2_label,
      fill = "elf5 mCG"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, 
                                margin = margin(b = 5)),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40",
                                   margin = margin(b = 15)),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "gray90", size = 0.5),
      panel.grid.minor = element_line(color = "gray95", size = 0.3),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "gray80", fill = NA, size = 0.8)
    ) +
    scale_x_continuous(limits = c(-20, 30), expand = c(0, 0))
)
dev.off()

##-- plot nanog mCG
# load bigtable
load("dunnart_E9_scBSseq_bigtable.RData")

# convert to GRanges object and extract relevant columns
CpG <- GRanges(bigtable, ignore.mcols=TRUE)
values(CpG) <- as.data.frame(bigtable)[,-c(1:5)]  # keep only methylation columns

# load NANOG promoter region from bed file
# chr5  187330636 187332636
nanog_promoter <- import.bed("nanog.bed")

# identify column pairs
col <- seq(1, ncol(mcols(CpG)), by=2)
sample_names <- sub("\\.(C|cov)$", "", colnames(mcols(CpG))[col])
col <- setNames(col, sample_names)

# set minimum coverage threshold (requiring at least 1 read per CpG)
minCov <- 1  

# calculate NANOG promoter methylation per sample
nanog_methylation <- numeric(length(col))
names(nanog_methylation) <- sample_names

for(i in seq_along(col)) {
  C_col <- col[i]
  cov_col <- C_col + 1
  
  # use overlapRatios function to get methylation ratio for NANOG promoter
  nanog_meth <- overlapRatios(
    x = CpG, 
    y = nanog_promoter, 
    C = colnames(mcols(CpG))[C_col], 
    cov = colnames(mcols(CpG))[cov_col], 
    minCov = minCov, 
    na.rm = TRUE
  )
  
  nanog_methylation[i] <- nanog_meth[1]  # First (and only) element
  
  if(i %% 20 == 0) {
    cat("Processed", i, "of", length(col), "samples\n")
  }
}

# create results data frame for NANOG promoter methylation
nanog_results <- data.frame(
  sample = names(nanog_methylation),
  nanog_methylation = nanog_methylation,
  stringsAsFactors = FALSE
)

# add plate and row information from sample names
nanog_results$Plate <- sub("\\.(A|B|C|D).*", "", nanog_results$sample)
nanog_results$Row <- sub(".*\\.", "", sub("\\.[0-9]+.*", "", nanog_results$sample))

# save results
print(head(nanog_results))
write.csv(nanog_results, "NANOG_promoter_methylation_per_sample.csv", row.names = FALSE)
saveRDS(nanog_results, "NANOG_promoter_methylation_per_sample.rds")

nanog_mCG <- read.csv("NANOG_promoter_methylation_per_sample.csv")
# extract sample name
pca_tbl <- pca_tbl %>%
  mutate(sample = str_replace(cell, ".*-(P\\d+)-([A-H]\\d{2})\\.bismark", "\\1.\\2"))
print(sum(pca_tbl$sample %in% nanog_mCG$sample))  # Should be > 0
# join with nanog_mCG
pca_merged <- pca_tbl %>%
  left_join(nanog_mCG, by = "sample")

# plot PCA coloured by nanog methylation
pdf("PCA_promoter_mCG_coloured_by_nanog.pdf", useDingbats = FALSE, width = 8, height = 6)
print(
  ggplot(pca_merged, aes(x = PC1, y = PC2)) +
    geom_point(data = pca_merged %>% filter(is.na(nanog_methylation)),
               shape = 21,           
               color = "black",      
               fill = "grey",
               size = 3, 
               stroke = 0.8) +       
    geom_point(data = pca_merged %>% filter(!is.na(nanog_methylation)),
               aes(fill = nanog_methylation), 
               shape = 21,           
               color = "black",      
               size = 3, 
               stroke = 0.8) +       
    scale_fill_gradientn(
      colours = c("#00008b", "white", "#D20A2E"),
      limits = c(0, 0.15),
      na.value = "grey"
    ) +
    coord_fixed() +
    labs(
      title = "PCA of promoter mCG - E9 dunnart",
      x = pc1_label,
      y = pc2_label,
      fill = "nanog mCG"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, 
                                margin = margin(b = 5)),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40",
                                   margin = margin(b = 15)),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "gray90", size = 0.5),
      panel.grid.minor = element_line(color = "gray95", size = 0.3),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "gray80", fill = NA, size = 0.8)
    ) +
    scale_x_continuous(limits = c(-20, 30), expand = c(0, 0))
)
dev.off()