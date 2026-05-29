install.packages("metafor")
install.packages("broom")
# Load libraries
library(dplyr)
library(metafor)
library(broom)
library(tidyr)

#-----------------------------
# Meta-analysis visualization guide
#-----------------------------
# - Vertical line = null effect
# - X-axis = effect size (case vs control)
# - Lines = individual features for a gene
# - Black box = point estimate (larger = more data)
# - Line = 95% confidence interval
# - Black diamond = fixed effect estimate (no moderator)
# - Gray diamond = estimate with moderator (e.g., group)

#-----------------------------
# Set directories
#-----------------------------
workdir <- "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType"

outPlot <- file.path(workdir, "quantify_t1d_pacbio_transcripts/meta_analysis_sex_diff_plots")
dir.create(outPlot, showWarnings = FALSE)

outSum <- file.path(workdir, "quantify_t1d_pacbio_transcripts/meta_analysis_sex_diff_summary")
dir.create(outSum, showWarnings = FALSE)

# Load T1D candidate gene list
t1d_df <- read.csv(file.path(workdir, "candidate_gene_lists/T1D_candidate_genes_robertson_2021_ensembl.txt"))
t1d_genes <- t1d_df$ensembl_geneID

#-----------------------------
# Loop through cell types
#-----------------------------
cell_list <- c("CD4", "CD8")

for (cell in cell_list) {
  
  inputFile <- paste0("ES_SV_", cell, "_pdiff_4_ma.csv")
  datafile <- file.path(workdir, "quantify_t1d_pacbio_transcripts", inputFile)
  datain <- read.csv(datafile)
  
  datain$geneID <- sapply(strsplit(datain$featureID, ":"), `[`, 1)
  
  # Compute effect size difference and average variance
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
  
  geneCol <- datain$geneID
  
  for (gene in unique(geneCol)) {
    
    genedf <- diff_df[diff_df$geneID == gene, ]
    
    num_features <- nrow(genedf)
    
    print(paste(cell, gene))
    
    if (num_features > 1 ) {
      
      # (1) Fixed-effects meta-analysis - no moderator
      dataout <- rma(
        yi = diff_ES,
        vi = variance_diff,
        data = genedf,
        slab = genedf$featureID,
        measure = "MD",
        method = "FE"
      )
      
      # Write summary output
      fileSum1 <- file.path(outSum, paste0("summary_ES_diff_", gene, "_", cell, "_FE.txt"))
      capture.output(summary(dataout), file = fileSum1)
      
      # Plot forest only if gene is a T1D candidate
      if (gene %in% t1d_genes) {
        outPNG <- paste0("forestPlot_ES_diff_", gene, "_", cell, "_FE.png")
        png(file.path(outPlot, outPNG), width = 960, height = 640)
        forest(
          dataout,
          header = TRUE,
          slab = genedf$featureID,
          xlab = "Male - Female ES, (MD, FE)"
        )
        dev.off()
      }
    }
  }
}
