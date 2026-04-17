library(data.table)
library(rtracklayer)
library(plyr)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)

# load bigtable 
load("dunnart_E9_scBSseq_bigtable.RData")

# convert to GRanges object and extract relevant columns
CpG <- GRanges(bigtable, ignore.mcols=TRUE)
values(CpG) <- as.data.frame(bigtable)[,-c(1:5)]  # keep only methylation columns

# identify column pairs
col <- seq(1, ncol(mcols(CpG)), by=2)
sample_names <- sub("\\.(C|cov)$", "", colnames(mcols(CpG))[col])
col <- setNames(col, sample_names)

# prepare metadata for mean methylation calculation
CpG_mcols <- mcols(CpG)
mean_methylation <- numeric(length(col))
names(mean_methylation) <- sample_names

# calculate mean methylation per sample
for(i in seq_along(col)) {
  C_col <- col[i]
  cov_col <- C_col + 1
  
  C_vals <- CpG_mcols[[C_col]]
  cov_vals <- CpG_mcols[[cov_col]]
  
  ratio <- C_vals / cov_vals
  mean_methylation[i] <- mean(ratio[is.finite(ratio)], na.rm = TRUE)
  
  if(i %% 10 == 0) {
    cat("Processed", i, "of", length(col), "samples\n")
  }
}

# create results data frame
mean_results <- data.frame(
  variable = names(mean_methylation),
  rating.mean = mean_methylation,
  stringsAsFactors = FALSE
)

# save results
print(mean_results)
write.csv(mean_results, "average_methylation_per_sample.csv", row.names = FALSE)
saveRDS(mean_results, "average_methylation_per_sample.rds")

# plot average methylation per sample
p1 <- ggplot(mean_results, aes(x = variable, y = rating.mean)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Average Methylation Per Cell",
       x = "Cell",
       y = "Mean Methylation Ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("average_methylation_barplot.pdf", p1, width = 20, height = 6)

# add plate and row information from sample names
mean_results$Plate <- sub("\\.(A|B|C|D).*", "", mean_results$variable)
mean_results$Row <- sub(".*\\.", "", sub("\\.[0-9]+.*", "", mean_results$variable))
saveRDS(mean_results, "average_methylation_per_cell.rds")

# plot grouped by plate
p2 <- ggplot(mean_results, aes(x = variable, y = rating.mean, fill = Plate)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylim(0, 1) +
  labs(title = "Average 5mC per Cell",
       x = "Cell",
       y = "Mean 5mCG/CG") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Plate, scales = "free_x")

ggsave("average_methylation_by_plate.pdf", p2, width = 20, height = 8)

# mean mCG percentage
mean_percentage <- mean_results
mean_percentage$rating.mean <- mean_percentage$rating.mean * 100
colnames(mean_percentage) <- c("sample", "mCG/CG", "Plate", "Row")
saveRDS(mean_percentage, "Dunnart_SC_mean_mCG.rds")

# plot
av_mc <- readRDS("average_methylation_per_cell.rds")
mean(av_mc$rating.mean)

heatmap_data <- matrix(av_mc$rating.mean, ncol = 1)
rownames(heatmap_data) <- av_mc$variable
colnames(heatmap_data) <- "mC"

row_annotation <- data.frame(
  Plate = av_mc$Plate,
  Row = substr(av_mc$Row, 1, 1)
)
rownames(row_annotation) <- av_mc$variable

annotation_colors <- list(
  Plate = c("P3" = "#2166ac", "P4" = "#762a83", "P5" = "#5aae61"),
  Row = c("A" = "#d73027", "B" = "#f46d43", "C" = "#74add1", "D" = "#4575b4")
)

# create row annotations
row_ha <- rowAnnotation(
  df = row_annotation,
  col = annotation_colors,
  width = unit(1, "cm")
)

# create barplot annotation
barplot_ha <- rowAnnotation(
  "mean mC" = anno_barplot(
    av_mc$rating.mean,
    ylim = c(0, 1),
    bar_width = 0.8,
    gp = gpar(fill = "steelblue"),
    width = unit(3, "cm")
  ),
  annotation_name_side = "top"
)

# create the heatmap with barplot
ht <- Heatmap(
  heatmap_data,
  name = "mCG/CG",
  col = colorRamp2(seq(0, 1, length.out = 101), colorRampPalette(c("blue", "white", "red"))(101)),
  rect_gp = gpar(col = "black", lwd = 1),  # Black border around cells
  heatmap_legend_param = list(
    at = seq(0, 1, by = 0.25),
    color_bar = "continuous",
    legend_height = unit(4, "cm")
  ),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  left_annotation = row_ha,
  right_annotation = barplot_ha,
  column_title = "single cell average mC",
  row_title = "cells"
)
pdf("dunnart_d9_scBSseq_mean_mC_heatmap.pdf", width = 5, height = 12)
draw(ht)
dev.off()