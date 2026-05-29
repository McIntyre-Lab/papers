install.packages("metafor")
install.packages("broom")

library(metafor)
library(broom)
library(dplyr)

# Setup work directories
workdir <- "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType"

outPlot <- file.path(workdir, "quantify_t1d_pacbio_transcripts/meta_analysis_DE_pdiffs_plots_means")
dir.create(outPlot, showWarnings = FALSE)
outSum <- file.path(workdir, "quantify_t1d_pacbio_transcripts/meta_analysis_DE_pdiffs_summary_means")
dir.create(outSum, showWarnings = FALSE)
outData <- file.path(workdir, "quantify_t1d_pacbio_transcripts/data_4_forest_plots_means")
dir.create(outData, showWarnings = FALSE)

# Import annotations of fragments
annotFile <- file.path(workdir, "allChr_isoseq_FSM_ISM_NIC_event_analysis_ef.csv")
annotIn <- read.csv(annotFile)
annot_subset <- annotIn[, c("ef_id", "er_id", "gene_id", "ef_ir_flag", "ef_chr", "ef_start", "ef_end", "ef_strand")]

#Inpor t1d list
t1d_df <- read.csv(file.path(workdir, "candidate_gene_lists/T1D_candidate_genes_robertson_2021_ensembl.txt"))
t1d_genes <- t1d_df$ensembl_geneID

cell_list <- c("CD4", "CD8")

for (cell in cell_list) {
  inputFile <- paste("ES_SV_", cell, "_pdiff_4_ma.csv", sep = "")
  datafile <- file.path(workdir, "quantify_t1d_pacbio_transcripts", inputFile)
  datain <- read.csv(datafile)
  
  datain_means <- datain %>%
    group_by(featureID) %>%
    summarize(
      # Weighted Mean Effect Size
      # Weighted Mean ES = (ES_m / Var_m + ES_f / Var_f) / (1 / Var_m + 1 / Var_f)
      wtd_mean_effect_size = sum(effect_size / variance) / sum(1 / variance),
      
      # Weighted Variance
      # Weighted Var = 1 / (1 / Var_m + 1 / Var_f)
      wtd_variance = 1 / sum(1 / variance),
      
      # Unweighted Mean Effect Size
      # Unweighted Mean ES = (ES_m + ES_f) / 2
      mean_effect_size = mean(effect_size),
      
      # Unweighted Variance
      # Unweighted Var = (Var_m + Var_f) / 4
      mean_variance = sum(variance) / n()^2,
      
      .groups = 'drop'
    )
  
  datain_means$geneID <- sapply(strsplit(datain_means$featureID, ":"), `[`, 1)
  
  geneCol <- datain_means$geneID
  for(gene in unique(geneCol)){
    genedf <- datain_means[datain_means$geneID == gene, ]
    gene_annot_df <- merge(genedf, annot_subset, by.x = "featureID", by.y = "ef_id", all.x = TRUE)
    
    num_features <- nrow(gene_annot_df)
    
    print(paste(cell, gene))
    
    if (num_features > 1) {
      # Meta-analysis using weighted means
      dataout_wtd <- rma(wtd_mean_effect_size, wtd_variance, data = gene_annot_df, method = "FE")
      out_wtd <- summary(dataout_wtd)
      fileSum_wtd <- file.path(outSum, paste("summary_wtd_", gene, "_", cell, "_FE.txt", sep = ""))
      capture.output(out_wtd, file = fileSum_wtd)
      
      # Meta-analysis using unweighted means
      dataout_unwtd <- rma(mean_effect_size, vi = mean_variance, data = gene_annot_df, method = "FE")
      out_unwtd <- summary(dataout_unwtd)
      fileSum_unwtd <- file.path(outSum, paste("summary_unwtd_", gene, "_", cell, "_FE.txt", sep = ""))
      capture.output(out_unwtd, file = fileSum_unwtd)
      
      if (gene %in% t1d_genes) {
        # Forest plot for weighted means
        outPNG_wtd <- paste("forestPlot_wtd_means_", gene, "_", cell, "_FE.png", sep = "")
        png(file.path(outPlot, outPNG_wtd), width = 960, height = 640)
        forest(dataout_wtd, header = TRUE, slab = gene_annot_df$featureID, xlab = "Case vs Control (MD, FE)")
        dev.off()
        
        # Forest plot for unweighted means
        outPNG_unwtd <- paste("forestPlot_unwtd_means_", gene, "_", cell, "_FE.png", sep = "")
        png(file.path(outPlot, outPNG_unwtd), width = 960, height = 640)
        forest(dataout_unwtd, header = TRUE, slab = gene_annot_df$featureID, xlab = "Case vs Control (MD, FE)")
        dev.off()
        
        # Extract weights and confidence intervals for individual studies (weighted)
        weights_wtd <- weights(dataout_wtd)
        proportional_weights_wtd <- weights_wtd / sum(weights_wtd)
        ci_lower_wtd <- gene_annot_df$wtd_mean_effect_size - 1.96 * sqrt(gene_annot_df$wtd_variance)
        ci_upper_wtd <- gene_annot_df$wtd_mean_effect_size + 1.96 * sqrt(gene_annot_df$wtd_variance)
        
        # Save the individual study data for custom plotting (weighted)
        plot_data_wtd <- data.frame(
          featureID = gene_annot_df$featureID,
          effect_size = gene_annot_df$wtd_mean_effect_size,
          variance = gene_annot_df$wtd_variance,
          weight = weights_wtd,
          proportional_weight = proportional_weights_wtd,
          ci_lower = ci_lower_wtd,
          ci_upper = ci_upper_wtd,
          ef_chr = gene_annot_df$ef_chr,
          ef_start = gene_annot_df$ef_start,
          ef_end = gene_annot_df$ef_end,
          ef_strand = gene_annot_df$ef_strand,
          ef_ir_flag = gene_annot_df$ef_ir_flag
        )
        
        # Output file path for individual study data (weighted)
        output_file_study_wtd <- file.path(outData, paste("plot_data_wtd_", gene, "_", cell, "_FE_study_wtd_means.csv", sep = ""))
        
        # Write individual study data to CSV (weighted)
        write.csv(plot_data_wtd, output_file_study_wtd, row.names = FALSE)
        
        # Extract weights and confidence intervals for individual studies (unweighted)
        weights_unwtd <- weights(dataout_unwtd)
        proportional_weights_unwtd <- weights_unwtd / sum(weights_unwtd)
        ci_lower_unwtd <- gene_annot_df$mean_effect_size - 1.96 * sqrt(gene_annot_df$mean_variance)
        ci_upper_unwtd <- gene_annot_df$mean_effect_size + 1.96 * sqrt(gene_annot_df$mean_variance)
        
        # Save the individual study data for custom plotting (unweighted)
        plot_data_unwtd <- data.frame(
          featureID = gene_annot_df$featureID,
          effect_size = gene_annot_df$mean_effect_size,
          variance = gene_annot_df$mean_variance,
          weight = weights_unwtd,
          proportional_weight = proportional_weights_unwtd,
          ci_lower = ci_lower_unwtd,
          ci_upper = ci_upper_unwtd,
          ef_chr = gene_annot_df$ef_chr,
          ef_start = gene_annot_df$ef_start,
          ef_end = gene_annot_df$ef_end,
          ef_strand = gene_annot_df$ef_strand,
          ef_ir_flag = gene_annot_df$ef_ir_flag
        )
        
        # Output file path for individual study data (unweighted)
        output_file_study_unwtd <- file.path(outData, paste("plot_data_unwtd_", gene, "_", cell, "_FE_study_unwtd_means.csv", sep = ""))
        
        # Write individual study data to CSV (unweighted)
        write.csv(plot_data_unwtd, output_file_study_unwtd, row.names = FALSE)
      }
    } else {
      print(paste(gene, "in", cell, "does not have enough features. Cannot estimate parameters"))
    }
  }
}
