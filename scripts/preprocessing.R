#!/usr/bin/env Rscript

#' Combine source expression files into a matrix
#'
#' @import tidyverse
#' @import purrr
#' @import dplyr
#' @import readr
#' @export
source_data_to_matrix <- function(exp_data_dir_path = "data/raw", output_matrix_path = "data/exp_matrix.txt", delim = ",") {
  exp_data <- list.files(exp_data_dir_path, full.names = TRUE)
  list_of_files <- lapply(exp_data, read_delim, delim = delim, escape_double = FALSE, trim_ws = TRUE)
  exp_matrix <- purrr::reduce(list_of_files, dplyr::full_join, by = "SYMBOL")
  write.table(exp_matrix, output_matrix_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

#' Load preprocessed expression matrix
#'
#' @param exp_matrix_path path to processed expression matrix obtained from \code{source_data_to_matrix}
#' @param drop_dupl_genes whether to remove duplicated gene names (ex. MIR*, LINC*). By default is \code{FALSE}
#' @param sample_list_path optional. Path to `.txt` file containing a list of TCGA patient codes to subset from the expression matrix.
#' @param quantile_normalization bool. Whether to conduct Quantile Normalization of expression data
#' @param log10_scalling bool. Wherher to convert counts to log10 values
#' @return expression matrix where cols are sample names and rows are SYMBOL genes
#'
#' @import readr
#' @import preprocessCore
#' @export
load_exp_matrix <- function(exp_matrix_path = "data/exp_matrix.txt",
                            drop_dupl_genes = FALSE, sample_list_path = NULL,
                            quantile_normalization = TRUE, log10_scalling = TRUE) {
  exp_matrix <- read_delim(exp_matrix_path,
        "\t", escape_double = FALSE, trim_ws = TRUE)
  exp_matrix <- as.data.frame(exp_matrix)

  if (drop_dupl_genes) {
    exp_matrix <- exp_matrix[!duplicated(exp_matrix$SYMBOL),]
  }

  exp_matrix <- exp_matrix[complete.cases(exp_matrix),]
  rownames(exp_matrix) <- exp_matrix$SYMBOL
  exp_matrix$SYMBOL <- NULL

  colnames(exp_matrix) <- samples_names_to_case_ids(sample_names = colnames(exp_matrix))

  if (!is.null(sample_list_path)) {
    sample_list <- read_delim(sample_list_path,
        "\t", col_names = FALSE, escape_double = FALSE, trim_ws = TRUE)
    common_samples_tumour <- intersect(colnames(exp_matrix), sample_list$X1)
    common_samples_norm <- intersect(colnames(exp_matrix), paste0("Norm_", sample_list$X1))
    common_samples <- unique(c(common_samples_tumour, common_samples_norm))
    exp_matrix <- exp_matrix[, common_samples]
  }

  if (quantile_normalization) {
    exp_matrix_genes <- rownames(exp_matrix)
    exp_matrix_samples <- colnames(exp_matrix)
    exp_matrix <- normalize.quantiles(as.matrix(exp_matrix + 1))
    exp_matrix <- as.data.frame(exp_matrix)
    rownames(exp_matrix) <- exp_matrix_genes
    colnames(exp_matrix) <- exp_matrix_samples
  }

  if (log10_scalling) {
    exp_matrix <- log10(exp_matrix)
  }

  return(exp_matrix)
}

#' Convert TCGA sample names into case IDs
#'
#' @param sample_names a vector, containing TCGA sample names
#' @return a vector, contating TCGA case IDs
#' @export
samples_names_to_case_ids <- function(sample_names) {
  case_ids <- gsub("(TCGA-.+-.+)-.*", "\\1", sample_names)
  return(case_ids)
}

#' Load annotations for expression samples
#'
#' @param sample_ann_table_path path to sample annotations table
#' @param check_gene_mutated HGNC gene name. If provided return annotation
#' for selected genes with either "Mutated" or "Control" values
#' @param sample_list_path optional. Path to `.txt` file containing a list of TCGA patient codes to subset from the expression matrix.
#' @param set_control_wo_any_genes_mut bool. Sets controls as samples without any mutations in `check_gene_mutated` AND without any mutations in other genes (EGFR, BRAF, *RAS)
#' @param set_control_as_norms bool. Sets controls as adjanced normal samples, that should be marked with "Norm_" prefix in sample annotation table
#' @return a named vector, where values are mutations (AA changes) and names are TCGA case IDs
#'
#' @import readxl
#' @import tidyverse
#' @import tidyr
#' @import dplyr
#' @import data.table
#' @export
get_mut_sample_annotations <- function(sample_ann_table_path = "data/TCGA_VCFs_found_AA_changes.xlsx",
                                       check_gene_mutated, check_aa_change, sample_list_path = NULL,
                                       set_control_wo_any_genes_mut = FALSE,
                                       set_control_as_norms = FALSE) {
  sample_ann_table <- read_excel(sample_ann_table_path, col_types = c("skip", "text", "text"))
  sample_ann_table$`Case ID` <- samples_names_to_case_ids(sample_names = sample_ann_table$`Case ID`)
  sample_ann_table <- separate_rows(sample_ann_table, Mutations, sep = "; ")
  sample_ann_table <- sample_ann_table[!duplicated(sample_ann_table),]
  sample_ann_table <- setDT(sample_ann_table)[, if (any(!is.na(Mutations))) .SD[!is.na(Mutations)] else .SD, by = `Case ID`]
  sample_ann_table <- as.data.frame(sample_ann_table)

  if (set_control_as_norms == FALSE) {
    sample_ann_table <- sample_ann_table %>% filter((Mutations != "normal_sample") %>% replace_na(TRUE))
  }

  if (!missing(check_gene_mutated) & missing(check_aa_change)) {
    if (set_control_wo_any_genes_mut == TRUE) {
      sample_ann_table <- sample_ann_table[is.na(sample_ann_table$Mutations) | grepl(check_gene_mutated, sample_ann_table$Mutations),]
    }
    if (set_control_as_norms == FALSE) {
      sample_ann_table$Mutations[is.na(sample_ann_table$Mutations)] <- "Control"
      sample_ann_table$Mutations[sample_ann_table$Mutations == "NA"] <- "Control"
    } else {
      sample_ann_table$Mutations[sample_ann_table$Mutations == "NA"] <- NA
      sample_ann_table <- sample_ann_table[!is.na(sample_ann_table$Mutations),]
      sample_ann_table <- sample_ann_table[grepl("normal_sample", sample_ann_table$Mutations) | grepl(check_gene_mutated, sample_ann_table$Mutations),]
      sample_ann_table$Mutations <- gsub("normal_sample", "Control", sample_ann_table$Mutations)
    }
    sample_ann_table$Mutations <- ifelse(grepl(check_gene_mutated, sample_ann_table$Mutations), "Mutated", "Control")
  }

  if (!missing(check_gene_mutated) & !missing(check_aa_change)) {
    if (set_control_wo_any_genes_mut == TRUE) {
      sample_ann_table <- sample_ann_table[is.na(sample_ann_table$Mutations) | grepl(paste(check_gene_mutated, check_aa_change, sep = " "), sample_ann_table$Mutations),]
    }
    if (set_control_as_norms == FALSE) {
      sample_ann_table$Mutations[is.na(sample_ann_table$Mutations)] <- "Control"
      sample_ann_table$Mutations[sample_ann_table$Mutations == "NA"] <- "Control"
    } else {
      sample_ann_table$Mutations[sample_ann_table$Mutations == "NA"] <- NA
      sample_ann_table <- sample_ann_table[!is.na(sample_ann_table$Mutations),]
      sample_ann_table <- sample_ann_table[grepl("normal_sample", sample_ann_table$Mutations) | grepl(paste(check_gene_mutated, check_aa_change, sep = " "), sample_ann_table$Mutations),]
      sample_ann_table$Mutations <- gsub("normal_sample", "Control", sample_ann_table$Mutations)
    }
    sample_ann_table$Mutations <- ifelse(grepl(paste(check_gene_mutated, check_aa_change, sep = " "), sample_ann_table$Mutations), "Mutated", "Control")
  }

  sample_ann_table <- setDT(sample_ann_table)[, if (all(Mutations == "Control")) .SD else Mutations <- "Mutated", by = `Case ID`]
  sample_ann_table <- sample_ann_table[!duplicated(sample_ann_table),]

  sample_ann <- sample_ann_table$Mutations
  names(sample_ann) <- sample_ann_table$`Case ID`

  if (!is.null(sample_list_path)) {
    sample_list <- read_delim(sample_list_path,
        "\t", col_names = FALSE, escape_double = FALSE, trim_ws = TRUE)
    common_samples_tumour <- intersect(names(sample_ann), sample_list$X1)
    common_samples_norm <- intersect(names(sample_ann), paste0("Norm_", sample_list$X1))
    common_samples <- unique(c(common_samples_tumour, common_samples_norm))
    sample_ann <- sample_ann[common_samples]
  }

  return(sample_ann)
}

#' Preprocess Drug Scores database into sample annotation
#'
#' @param sample_list_path optional. Path to `.txt` file containing a list of TCGA patient codes to subset from the expression matrix.
#'
#' @import readr
#' @export
get_drug_scores_sample_annotations <- function(drug_scores_db_path = "data/drug_scores/ds_balanced.csv",
                                               drug_name = "Cetuximab",
                                               return_ranks = FALSE,
                                               sample_list_path = NULL) {
  drug_scores_db <- read_csv(drug_scores_db_path)
  colnames(drug_scores_db) <- samples_names_to_case_ids(sample_names = colnames(drug_scores_db))
  sample_cols <- !grepl("Norm_TCGA_pseudo|Drug|Database|Drug Type``", colnames(drug_scores_db), ignore.case = TRUE)
  if (return_ranks == TRUE) {
    drug_scores_db[, sample_cols] <- apply(drug_scores_db[, sample_cols], 2, rank)
  }

  drug_scores_db <- drug_scores_db[drug_scores_db$Drug == drug_name, sample_cols]

  sample_ann <- as.vector(t(drug_scores_db))
  names(sample_ann) <- gsub("Tumour_(TCGA\\.[0-9a-zA-Z]+\\.[0-9a-zA-Z]+)\\..*", "\\1", colnames(drug_scores_db))
  names(sample_ann) <- gsub("\\.", "-", names(sample_ann))
  names(sample_ann) <- gsub("Tumour_", "", names(sample_ann), ignore.case = TRUE) # if Drug Scores DB came from BES calculations

  if (!is.null(sample_list_path)) {
    sample_list <- read_delim(sample_list_path,
        "\t", col_names = FALSE, escape_double = FALSE, trim_ws = TRUE)
    common_samples_tumour <- intersect(names(sample_ann), sample_list$X1)
    common_samples_norm <- intersect(names(sample_ann), paste0("Norm_", sample_list$X1))
    common_samples <- unique(c(common_samples_tumour, common_samples_norm))
    sample_ann <- sample_ann[common_samples]
  }

  sample_names <- names(sample_ann)
  sample_ann <- as.numeric(sample_ann)
  names(sample_ann) <- sample_names
  return(sample_ann)
}

#' Preprocess Pathways Activation database into sample annotation
#'
#' @import readxl
#' @return a list, which elements are named as pathway names and
#' each element itself correspond to a named vector, where values are
#' activations of a given pathway in samples (TCGA case IDs) and names are TCGA case IDs
#' @export
get_pathways_activation_sample_annotations <- function(pathways_activation_db_path = "data/pathways_activation_scores/EGF_PWs.xlsx",
                                                       return_dataframe = FALSE,
                                                       return_ranks = FALSE) {
  pathways_activation_db <- read_excel(pathways_activation_db_path)
  colnames(pathways_activation_db) <- samples_names_to_case_ids(sample_names = colnames(pathways_activation_db))
  pathways_names <- paste(pathways_activation_db$Pathway, pathways_activation_db$Database, sep = " - ")
  sample_cols <- !grepl("Norm_TCGA_pseudo", colnames(pathways_activation_db), ignore.case = TRUE)

  if (return_ranks == TRUE) {
    pathways_activation_db[, sample_cols] <- apply(pathways_activation_db[, sample_cols], 2, rank)
  }

  pathways_activation_db <- pathways_activation_db[, sample_cols]
  colnames(pathways_activation_db) <- gsub("Tumour_(TCGA\\.[0-9a-zA-Z]+\\.[0-9a-zA-Z]+)\\..*", "\\1", colnames(pathways_activation_db))
  colnames(pathways_activation_db) <- gsub("\\.", "-", colnames(pathways_activation_db))
  colnames(pathways_activation_db) <- gsub("Tumour_", "", colnames(pathways_activation_db), ignore.case = TRUE) # if Drug Scores DB came from BES calculations
  sample_names <- colnames(pathways_activation_db)

  if (return_dataframe == FALSE) {
    sample_ann_list <- split(pathways_activation_db, pathways_names)
    sample_ann_list <- lapply(sample_ann_list, as.numeric)
    sample_ann_list <- lapply(sample_ann_list, setNames, sample_names)
    sample_ann <- sample_ann_list
  } else {
    sample_ann <- pathways_activation_db
    rownames(sample_ann) <- pathways_names
  }

  return(sample_ann)
}

#' Preprocess DESeq2 results into sample annotation for gene signature of significantly differentially expressed genes
#'
#' @param dif_exp_results_path path to .RDS object corresponding to the DESeq2 results data frame
#' @param candidate_genes character vector, containing unique HUGO (HGNC) name for genes,
#' that will be used to form a signature if their differential expression passes \code{pCutoff} and \code{pCutoff} thresholds
#' @param expression_matrix a data frame containing preprocessed values for gene expression. See \code{load_exp_matrix} function
#' @param n_genes_per_signature numeric. If specified it will generate multiple signatures from a set of passed `candidate_genes`,
#' where each signature will have a size of specified number.
#' @return a list, which elements are named as signature name and
#' each element itself correspond to a named vector, where values are 
#' scores (log2FC weighted expression) of a given signature in samples (TCGA case IDs) and names are TCGA case IDs
#'
#' @import DESeq2
#' @import readxl
#' @export
get_gene_signature_sample_annotations <- function(dif_exp_results_path,
                                                  candidate_genes,
                                                  expression_matrix,
                                                  pCutoff = 0.05,
                                                  FCcutoff = 0.5,
                                                  n_genes_per_signature = NULL) {
  dif_exp_results <- readRDS(dif_exp_results_path)
  dif_exp_results <- as.data.frame(dif_exp_results)
  dif_exp_results <- dif_exp_results[candidate_genes,]

  genes_passed <- rownames(dif_exp_results)[abs(dif_exp_results["log2FoldChange"]) > FCcutoff & dif_exp_results["padj"] < pCutoff]
  genes_passed <- genes_passed[!is.na(genes_passed)]

  # Calculate gene signature scores
  signature_sample_ann <- list()

  signature_all_passed_genes <- calculate_gene_signature_scores(expression_matrix, dif_exp_results,
                                                                signature_genes = genes_passed)

  signature_sample_ann <- c(signature_sample_ann, signature_all_passed_genes)

  if (!is.null(n_genes_per_signature)) {
    genes_combinations <- combn(genes_passed, n_genes_per_signature)
    for (i in 1:ncol(genes_combinations)) {
      signature_i_genes_comb <- calculate_gene_signature_scores(expression_matrix, dif_exp_results,
                                                                signature_genes = genes_combinations[, i])
      signature_sample_ann <- c(signature_sample_ann, signature_i_genes_comb)
    }
  }

  return(signature_sample_ann)
}

#' Calculates scores for gene signature
#'
#' Calculate gene signature score as an average weighed expression for
#' selected genes, where weights are corresponding gene log2FC
#'
#' @return a list, which elements are named as signature name and
#' each element itself correspond to a named vector, where values are
#' scores (log2FC weighted expression) of a given signature in samples (TCGA case IDs) and names are TCGA case IDs
#' @export
calculate_gene_signature_scores <- function(exp_matrix, dif_exp_results, signature_genes) {
  signature_name <- paste(signature_genes, collapse = " + ")
  signature_genes_exp <- exp_matrix[signature_genes,]
  signature_genes_log2FC_coef <- dif_exp_results[signature_genes,]$log2FoldChange
  signature_scores <- apply(signature_genes_log2FC_coef * signature_genes_exp, 2, mean)
  signature <- list(signature_scores)
  names(signature) <- signature_name
  return(signature)
}

#' Filter non-common samples
#'
#' Update in a global environment \code{exp_matrix} and \code{sample_ann}
#' by versions having matching sample names
#'
#' @param exp_data a data frame, where rows are HGNC genes (HUGO)
#' and cols are TCGA sample names (case IDs)
#' @param sample_ann a named vector, where names are TCGA sample names (case IDs)
#' @export
filter_non_common_samples <- function(exp_matrix, sample_ann) {
  exp_matrix_glob_name <- deparse(substitute(exp_matrix))
  sample_ann_glob_name <- deparse(substitute(sample_ann))
  common_samples <- intersect(colnames(exp_matrix), names(sample_ann))
  exp_matrix <- exp_matrix[, common_samples]
  sample_ann <- sample_ann[common_samples]

  data <- list(exp_matrix, sample_ann)
  names(data) <- c(exp_matrix_glob_name, sample_ann_glob_name)
  list2env(data, envir = .GlobalEnv)
}

#' Perform Differential Expression Analysis
#'
#' @param exp_data a data frame, where rows are HGNC genes (HUGO)
#' and cols are TCGA sample names (case IDs)
#' @param sample_ann a named vector, where names are TCGA sample names (case IDs)
#' @return a DESeq2 results object
#'
#' @import DESeq2
#' @import apeglm
#' @import EnhancedVolcano
#' @export
dif_exp_analysis <- function(exp_matrix, sample_ann, shring_results = FALSE) {
  dds <- DESeqDataSetFromMatrix(countData = exp_matrix,
                                colData = data.frame(condition = as.factor(sample_ann),
                                                    row.names = names(sample_ann)),
                                design = ~condition)
  dds <- DESeq(dds)
  print("Calculating Differential Expression")
  res <- results(dds, name = "condition_Mutated_vs_Control")

  if (shring_results == TRUE) {
    print("Shrinking log fold changes association with condition")
    res <- lfcShrink(dds, coef = "condition_Mutated_vs_Control", type = "apeglm")
  }

  return(res)
}

#' Save plot of provided analysis step
#'
#' Export ggplot2 object to desired location in vector (`.pdf`) and raster (`.tiff`) formats
#' @import ggplot2
#' @import ggpubr
#' @import extrafont
#' @export
export_analysis_plot <- function(filename, plot, path,
                                 scale, width, height, units, load_fonts = FALSE) {
  if (load_fonts == TRUE) {
    extrafont::font_import(prompt = FALSE)
  }

  extrafont::loadfonts(device = "pdf")

  ggsave(
    filename = paste0(filename, ".pdf"),
    plot = plot,
    device = cairo_pdf,
    path = path,
    scale = scale,
    width = width,
    height = height,
    units = units,
    dpi = 300)

  ggsave(
    filename = paste0(filename, ".tiff"),
    plot = plot,
    device = "tiff",
    path = path,
    scale = scale,
    width = width,
    height = height,
    units = units,
    compression = "lzw",
    dpi = 1200)
}

#' Calculate AUC ROC per each pathway based on pathways activation levels (PAL) and return ROC curve chart
#'
#' @param sample_names a character vector, containing TCGA patient case id codes
#' @param is_mutated a character vector, containing either "Control" or "Mutated" values
#' @param pathways_names a character vector, contatining pathways names per each `sample_name`
#' @param pw_activation_score a numeric vector, containing PAL scores per each `sample_name` and `pathways_names` pair
#' @param gene SYMBOL of analysed gene
#'
#' @return list of arranged ROC plots for each of pathways
#' @import pROC
#' @import dplyr
#' @import rlist
#' @export
calculate_roc_per_pw <- function(sample_names, is_mutated, pathways_names, pw_activation_score,
                                 gene, aa_change = NULL) {
  data <- data.frame(sample_names = as.character(sample_names),
                     is_mutated = as.factor(is_mutated),
                     pathways_names = as.character(gsub(" - ", "\n", pathways_names)),
                     pw_activation_score = as.numeric(pw_activation_score))

  plots_list <- list()

  for (pathway_to_analyse in unique(data$pathways_names)) {
    pw_data <- data[data$pathways_names %in% pathway_to_analyse,]

    res.roc <- roc(pw_data$is_mutated, pw_data$pw_activation_score)

    # Collect ROC data
    roc.data <- data.frame(thresholds = res.roc$thresholds,
                           sensitivity = res.roc$sensitivities,
                           specificity = res.roc$specificities)

    roc.data <- roc.data %>%
    filter(thresholds != -Inf) %>%
    mutate(pathways_names = pw_data$pathways_names)

    # Plot ROC Curve per each pathway
    auc_roc_label <- paste("AUC:", round(res.roc[["auc"]], 3))
    roc_curve <- plot_roc_curve(roc.data, gene = gene,
                                          aa_change = aa_change,
                                          auc_roc_label = auc_roc_label,
                                          title = pathway_to_analyse)
    plots_list <- list.append(plots_list, roc_curve)
  }
  return(plots_list)
}

#' Plot ROC Curve
#'
#' @import ggplot2
#' @import ggpubr
#' @import grid
#' @import gridExtra
#' @import extrafont
#' @export
plot_roc_curve <- function(roc.data, gene, aa_change = NULL, auc_roc_label = NULL, title = "Pathways Activation") {
  if (is.null(aa_change)) {
    name_label <- paste("Is", gene, "mutated")
  } else {
    name_label <- paste("Is", gene, aa_change, "mutated")
  }

  auc_roc_plot <- ggplot(roc.data, aes(specificity, sensitivity)) +
                  geom_path(aes(color = "red"), size = 1.25) +
                  scale_x_reverse(expand = c(0, 0), limits = c(1, 0)) +
                  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
                  geom_abline(size = 1, intercept = 1, slope = 1, linetype = "dashed") +
                  annotate("text", x = 0.3, y = 0.3, size = 5, label = auc_roc_label) +
                  labs(title = title,
                      x = "1 - Specificity (FPR)",
                      y = "Sensitivity (TPR)") +
                      theme(text = element_text(size = 14)) +
                      theme_Publication() +
                      scale_fill_Publication() +
                      theme(legend.position = "none")
  return(auc_roc_plot)
}

#' Calculate AUC ROC per each signature based on signature scores and return ROC curve chart
#'
#' @param sample_names a character vector, containing TCGA patient case id codes
#' @param is_mutated a character vector, containing either "Control" or "Mutated" values
#' @param signature_names a character vector, contatining signature names per each `sample_name`
#' @param signature_scores a numeric vector, containing PAL scores per each `sample_name` and `signature_names` pair
#' @param gene SYMBOL of analysed gene
#'
#' @return list of arranged ROC plots for each of pathways
#' @import pROC
#' @import dplyr
#' @import rlist
#' @export
calculate_roc_per_signature <- function(sample_names, is_mutated, signature_names, signature_scores,
                                        gene, aa_change = NULL) {
  data <- data.frame(sample_names = as.character(sample_names),
                     is_mutated = as.factor(is_mutated),
                     signature_names = as.character(gsub(" - ", "\n", signature_names)),
                     signature_scores = as.numeric(signature_scores))

  plots_list <- list()

  for (signature_to_analyse in unique(data$signature_names)) {
    signature_data <- data[data$signature_names %in% signature_to_analyse,]

    res.roc <- roc(signature_data$is_mutated, signature_data$signature_scores)

    # Collect ROC data
    roc.data <- data.frame(thresholds = res.roc$thresholds,
                           sensitivity = res.roc$sensitivities,
                           specificity = res.roc$specificities)

    roc.data <- roc.data %>%
    filter(thresholds != -Inf) %>%
    mutate(signature_names = signature_data$signature_names)

    # Plot ROC Curve per each pathway
    auc_roc_label <- paste("AUC:", round(res.roc[["auc"]], 3))
    roc_curve <- plot_roc_curve(roc.data, gene = gene,
                                          aa_change = aa_change,
                                          auc_roc_label = auc_roc_label,
                                          title = signature_to_analyse)
    plots_list <- list.append(plots_list, roc_curve)
  }
  return(plots_list)
}

#' Maps ENSEMBL gene names into HUGO (HGNC / SYMBOL) format and vice versa
#'
#' @param gene_names a vector, containing gene names in ENSEMBL format
#' @param gene_format a str, specifying either `ensembl_gene_id` or `hgnc_symbol` input gene format
#' @param counts a vector, containig raw counts for gene expression
#' @param hsmart_ann a R object, containig HSMART ENSEMBL genes annotation
#' @return a data frame with columns: `ensembl_gene_id`, `SYMBOL`, `counts`
#'
#' @import biomaRt
#' @export
convert_ensembl_to_symbol_gene_names <- function(gene_names, gene_format = c("ensembl_gene_id", "hgnc_symbol"),
                                                 counts, hsmart_ann, ...) {
  if (missing(hsmart_ann)) {
    hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
  }

  if (gene_format == "ensembl_gene_id") {
    gene_names <- sapply(gene_names, function(x) { gsub("\\.[0-9]*", "", x) })
  }

  mapping <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    values = gene_names,
    mart = hsmart
  )

  data <- data.frame(gene = gene_names, counts = counts)
  colnames(data)[1] <- match.arg(gene_format)

  mapping <- mapping[!duplicated(mapping$hgnc_symbol),]
  data <- merge(data, mapping, by = match.arg(gene_format))
  data <- data[c("ensembl_gene_id", "hgnc_symbol", "counts")]
  return(data)
}

#' Process an expression sample from Oncobox FASTQ STAR pipeline into need format
#'
#' @param exp_sample_path a str, specifying path to input expression sample
#' @param hsmart_ann a R object, containig HSMART ENSEMBL genes annotation
#' @param output_dir a str, specifying path to output dir
#' @return a data frame with columns: `ensembl_gene_id`, `SYMBOL`, `<sample_code>`
#'
#' @import readr
#' @export
preprocess_fastq_star_exp_sample <- function(exp_sample_path, hsmart_ann, output_dir) {
  sample_name <- gsub("_paired_endReadsPerGene.out.tab", "", basename(exp_sample_path))
  exp_sample <- read_delim(exp_sample_path, delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  colnames(exp_sample) <- c("ENSEMBL", "R1_R2_counts", "R1_counts", "R2_counts")
  exp_sample <- exp_sample[c("ENSEMBL", "R1_R2_counts")][grepl("ENSG[0-9]+.*", exp_sample$ENSEMBL),]

  exp_sample <- convert_ensembl_to_symbol_gene_names(gene_names = exp_sample$ENSEMBL,
                                                     gene_format = "ensembl_gene_id",
                                                     counts = exp_sample$R1_R2_counts,
                                                     hsmart_ann = hsmart_ann)
  exp_sample <- exp_sample[c("hgnc_symbol", "counts")]
  colnames(exp_sample) <- c("SYMBOL", sample_name)

  output_file_path <- file.path(output_dir, basename(exp_sample_path))
  write.table(exp_sample, output_file_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

#' Process in a given directory all expression samples from Oncobox FASTQ STAR pipeline into need format and saves output to dir
#'
#' @import biomaRt
#' @export
preprocess_fastq_star_samples_dir <- function(input_exp_data_dir_path = "data/literature_validation_samples_data/exp_data/raw",
                                              output_exp_data_dir_path = "data/literature_validation_samples_data/exp_data/processed") {
  exp_data <- list.files(input_exp_data_dir_path, full.names = TRUE)
  print("Loading HSMART gene annotations...")
  hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
  print("Starting preprocessing...")
  lapply(exp_data, preprocess_fastq_star_exp_sample, hsmart_ann = hsmart, output_dir = output_exp_data_dir_path)
}
