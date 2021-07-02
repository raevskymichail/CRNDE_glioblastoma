# setwd('C:/Users/raevs/Desktop/Oncobox/QC_for_preprint')
setwd('/mnt/c/Users/raevs/Desktop/Oncobox/QC_for_preprint')

source('bin/utils_gallow_plot.R', chdir = TRUE)
source('bin/preprocessing.R', chdir = TRUE)

if (!require("pacman")) install.packages("pacman")
pacman::p_load("magrittr", "stringr", "readr", "readxl", "FactoMineR", "factoextra", "ggplot2", "ggpubr", "gridExtra", "ggdendro", "dendroextras")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("preprocessCore")) BiocManager::install("preprocessCore")

# Read the expression matrix
exp_matrix <- load_exp_matrix(exp_matrix_path = "data/exp_matrix.txt")

# Find Num. uniquely mapped reads
mapped_reads_vec <- collect_QC_stat(exp_matrix)

# Obtain Sample source
clin_data <- process_clinical_info(clin_info_path = "clinical info.xlsx", exclude_control = FALSE)
sample_source <- obtain_samples_source(clin_data)

# Filter low-quality samples (unique_mapped_reads < 2.5e+6)
samples_data_list <- filter_low_QC_samples(exp_matrix, mapped_reads_vec, sample_source, qc_threshold = 2.5e+6)
exp_data_sample_names <- colnames(samples_data_list[["exp_matrix"]])

# Get response per samples
response_per_sample <- get_response_per_sample(response_info_path = "data/Suppl1_clear_MS.xlsx",
                                               exp_data_sample_names = colnames(samples_data_list[["exp_matrix"]]))

# Subset samples without response info
samples_data_list[["exp_matrix"]] <- samples_data_list[["exp_matrix"]][, names(response_per_sample)]

# Run PCA
data.pca <- PCA(t(samples_data_list[["exp_matrix"]]), graph = FALSE)

pcs_variance <- as.data.frame(data.pca$eig)$`percentage of variance`
pc1_lab <- paste0("PC1 (", round(pcs_variance[1], 2), "%)")
pc2_lab <- paste0("PC2 (", round(pcs_variance[2], 2), "%)")

p <- fviz_pca_ind(data.pca,
                  geom.ind = "point", # show points only (nbut not "text")
                  col.ind = response_per_sample, # color by groups
                  addEllipses = FALSE, # Concentration ellipses
                  legend.title = NULL
                  ) +
                  labs(x = pc1_lab, y = pc2_lab) +
theme(legend.title = element_blank(),
                        plot.title = element_blank())

ggsave(
  filename = "pca_plot_by_response.pdf",
  plot = p,
  device = "pdf",
  path = "plots/pca_plot_by_response",
  scale = 1,
  width = 6,
  height = 4,
  units = "in",
  dpi = 300)

ggsave(
  filename = "pca_plot_by_response.tiff",
  plot = p,
  device = "tiff",
  path = "plots/pca_plot_by_response",
  scale = 1,
  width = 6,
  height = 4,
  units = "in",
  compression = "lzw",
  dpi = 1200)

# Prepare plot for A4 panel

# ggsave(
#   filename = "pca_plot_by_response_for_panel.pdf",
#   plot = p,
#   device = "pdf",
#   path = "plots/panel_plot",
#   scale = 1,
#   width = 210 / 2,
#   height = 297 / 3,
#   units = "mm",
#   dpi = 300)
