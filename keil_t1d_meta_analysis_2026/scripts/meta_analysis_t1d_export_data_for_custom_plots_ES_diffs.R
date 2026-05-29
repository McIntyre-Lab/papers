# -----------------------------
# Setup
# -----------------------------
# install.packages("metafor")
# install.packages("dplyr")
# install.packages("tidyr")

library(metafor)
library(dplyr)
library(tidyr)

# -----------------------------
# Set directories
# -----------------------------
workdir <- "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType"

outData <- file.path(workdir, "quantify_t1d_pacbio_transcripts/data_4_forest_plots_ES_diff")
dir.create(outData, showWarnings = FALSE)

# -----------------------------
# Load data
# -----------------------------
# T1D candidate genes
t1d_df <- read.csv(file.path(workdir, "candidate_gene_lists/T1D_candidate_genes_robertson_2021_ensembl.txt"))
t1d_genes <- t1d_df$ensembl_geneID

# Annotation data
annotFile <- file.path(workdir, "allChr_isoseq_FSM_ISM_NIC_event_analysis_ef.csv")
annotIn <- read.csv(annotFile)
annot_subset <- annotIn[, c("ef_id", "er_id", "gene_id", "ef_ir_flag", "ef_chr", "ef_start", "ef_end", "ef_strand")]

# -----------------------------
# Meta-analysis per cell type
# -----------------------------
cell_list <- c("CD4", "CD8")

for (cell in cell_list) {
  message("Processing cell type: ", cell)
  
  # Load input data
  inputFile <- paste0("ES_SV_", cell, "_pdiff_4_ma.csv")
  datafile <- file.path(workdir, "quantify_t1d_pacbio_transcripts", inputFile)
  datain <- read.csv(datafile)
  datain$geneID <- sapply(strsplit(datain$featureID, ":"), `[`, 1)
  
  # Compute effect size differences (M - F) and variance
  effect_df <- datain %>%
    select(featureID, geneID, sex, effect_size) %>%
    pivot_wider(names_from = sex, values_from = effect_size)
  
  variance_df <- datain %>%
    select(featureID, sex, variance) %>%
    pivot_wider(names_from = sex, values_from = variance) %>%
    rename(F_var = F, M_var = M)
  
  diff_df <- left_join(effect_df, variance_df, by = "featureID") %>%
    filter(!is.na(F) & !is.na(M)) %>%
    mutate(
      diff_ES = M - F,
      variance_diff = (M_var + F_var) / 2
    ) %>%
    select(featureID, geneID, diff_ES, variance_diff)
  
  # Filter only T1D genes
  diff_df <- diff_df %>% filter(geneID %in% t1d_genes)
  
  for (gene in unique(diff_df$geneID)) {
    genedf <- diff_df %>% filter(geneID == gene)
    num_features <- nrow(genedf)
    
    # Join annotations
    gene_annot_df <- merge(genedf, annot_subset, by.x = "featureID", by.y = "ef_id", all.x = TRUE)
    
    if (num_features > 1) {
      message("Running meta-analysis for T1D gene: ", gene)
      
      # Meta-analysis (fixed effect)
      dataout <- rma(
        yi = diff_ES,
        vi = variance_diff,
        data = gene_annot_df,
        slab = gene_annot_df$featureID,
        method = "FE"
      )
      
      # Add meta-analysis statistics
      gene_annot_df$weight <- weights(dataout)
      gene_annot_df$proportional_weight <- gene_annot_df$weight / sum(gene_annot_df$weight)
      gene_annot_df$ci_lower <- gene_annot_df$diff_ES - 1.96 * sqrt(gene_annot_df$variance_diff)
      gene_annot_df$ci_upper <- gene_annot_df$diff_ES + 1.96 * sqrt(gene_annot_df$variance_diff)
      
      # Export annotated per-feature data
      feature_plot_data <- gene_annot_df %>%
        select(
          featureID, geneID, diff_ES, variance_diff,
          weight, proportional_weight, ci_lower, ci_upper,
          er_id, ef_ir_flag, ef_chr, ef_start, ef_end, ef_strand
        ) %>%
        rename(
          effect_size = diff_ES,
          variance = variance_diff
        )
      
      feature_file <- file.path(outData, paste0("plot_data_ES_diff_", gene, "_", cell, "_features_FE.csv"))
      write.csv(feature_plot_data, feature_file, row.names = FALSE)
      
      # Save overall estimate
      overall_data <- data.frame(
        featureID = "Overall",
        geneID = gene,
        effect_size = dataout$b[1],
        variance = dataout$se[1]^2,
        weight = 1 / (dataout$se[1]^2),
        proportional_weight = 1,
        ci_lower = dataout$ci.lb[1],
        ci_upper = dataout$ci.ub[1],
      )
      
      overall_file <- file.path(outData, paste0("plot_data_ES_diff_", gene, "_", cell, "_overall_FE.csv"))
      write.csv(overall_data, overall_file, row.names = FALSE)
    } else {
      message("Skipping gene ", gene, " in ", cell, ": not enough features")
    }
  }
}
