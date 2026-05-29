# install.packages("metafor")
library(metafor)


#Set up directories
workdir <- "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType"

outSum  <- file.path(workdir, "quantify_t1d_pacbio_transcripts/meta_analysis_DE_pdiffs_summary_bysex")
dir.create(outSum, showWarnings = FALSE, recursive = TRUE)

# outPlot <- file.path(workdir, "quantify_t1d_pacbio_transcripts/meta_analysis_DE_pdiffs_plots")
# dir.create(outPlot, showWarnings = FALSE, recursive = TRUE)

#Input files
annotFile <- file.path(workdir, "allChr_isoseq_FSM_ISM_NIC_event_analysis_ef.csv")
annotIn   <- read.csv(annotFile)

annot_subset <- annotIn[, c("ef_id", "ef_ir_flag", "ef_start")]

t1d_df    <- read.csv(file.path(workdir, "candidate_gene_lists/T1D_candidate_genes_robertson_2021_ensembl.txt"))
t1d_genes <- t1d_df$ensembl_geneID 



cell_list <- c("CD4", "CD8")

# -----------------------------------------
# Helper to safely run rma and write summary
# -----------------------------------------
run_and_write <- function(df, gene, cell, sex_label, out_dir) {
  # Require at least 2 features to estimate FE meta-analysis
  if (nrow(df) < 2) {
    message(sprintf("%s in %s (%s) has <2 features; skipping FE meta-analysis.", gene, cell, sex_label))
    return(invisible(NULL))
  }
  
  # Ensure needed columns are present
  needed <- c("effect_size", "variance")
  if (!all(needed %in% names(df))) {
    stop(sprintf("Missing required columns (%s) for %s %s %s",
                 paste(needed, collapse=", "), gene, cell, sex_label))
  }
  
  # Sort for nicer, reproducible slab order
  df <- df[order(df$ef_start, df$featureID), ]
  
  # Fixed-effects meta-analysis (no moderator), MD
  fit <- tryCatch(
    rma(effect_size, variance,
        data = df,
        slab = df$featureID,
        measure = "MD",
        method  = "FE"),
    error = function(e) e
  )
  
  outfile <- file.path(out_dir, sprintf("summary_%s_%s_%s_FE.txt", gene, cell, sex_label))
  
  if (inherits(fit, "error")) {
    msg <- sprintf("ERROR for %s %s %s: %s", gene, cell, sex_label, fit$message)
    writeLines(msg, con = outfile)
    message(msg)
  } else {
    sink(outfile)
    print(summary(fit))
    sink()
  }
}

# -----------------------------------------
# Main loop
# -----------------------------------------
for (cell in cell_list) {
  inputFile <- paste0("ES_SV_", cell, "_pdiff_4_ma.csv")
  datafile  <- file.path(workdir, "quantify_t1d_pacbio_transcripts", inputFile)
  
  datain <- read.csv(datafile)
  datain$geneID <- sapply(strsplit(datain$featureID, ":"), `[`, 1)
  
  # Merge in annotations
  datain <- merge(datain, annot_subset, by.x = "featureID", by.y = "ef_id", all.x = TRUE)
  
  # Iterate genes
  for (gene in unique(datain$geneID)) {
    genedf <- datain[datain$geneID == gene, , drop = FALSE]
    
    # Split by sex
    for (sex_now in c("M", "F")) {
      sexdf <- genedf[genedf$sex == sex_now, , drop = FALSE]
      # Only attempt if that sex exists in this gene
      if (nrow(sexdf) > 0) {
        run_and_write(sexdf, gene, cell, sex_now, outSum)
      } else {
        message(sprintf("%s in %s has no features for sex %s.", gene, cell, sex_now))
      }
    }
  }
}

message("Done. Summaries written to: ", outSum)