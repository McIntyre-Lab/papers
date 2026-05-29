install.packages("metafor")
install.packages("broom")

  # vertical line = "line of null effect"
  # horizontal axis = effect size, case vs control
  # for each gene, study lines = features for a gene
    # black box = point estimate for feature, bigger box = more data points for that feature
    # line = 95% confidence interval
  # black diamond (no moderator) = point estimate and confidence intervals when combine and average all features together
  # grey diamond (moderator = group) = estimate with 
workdir <- "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType"

outPlot <- file.path(workdir, "quantify_t1d_pacbio_transcripts/meta_analysis_DE_pdiffs_plots_exonsOnly")
dir.create(outPlot)
outSum  <- file.path(workdir, "quantify_t1d_pacbio_transcripts/meta_analysis_DE_pdiffs_summary_exonsOnly")
dir.create(outSum)
outData  <- file.path(workdir, "quantify_t1d_pacbio_transcripts/data_4_forest_plots")
dir.create(outData)

library(metafor)
library(broom)

#Import annotations of fragments
annotFile=file.path(workdir,"allChr_isoseq_FSM_ISM_NIC_event_analysis_ef.csv")
annotIn=read.csv(annotFile)

annot_subset=annotIn[, c("ef_id","er_id","gene_id","ef_ir_flag","ef_chr","ef_start", "ef_end","ef_strand")]
  
t1d_df=read.csv(paste(workdir,"/candidate_gene_lists/T1D_candidate_genes_robertson_2021_ensembl.txt", sep=""))
t1d_genes=t1d_df$ensembl_geneID

cell_list <- c("CD4", "CD8")   #create character vector

for (cell in cell_list) {
  geneLst <- c('ENSG00000072110','ENSG00000112182',
               'ENSG00000134242','ENSG00000118503','ENSG00000122484')
  
  for (gene in t1d_genes) {
    inputFile <- paste("ES_SV_", cell, "_pdiff_4_ma.csv", sep = "")
    
    datafile <- file.path(workdir, "quantify_t1d_pacbio_transcripts", inputFile)
    datain <- read.csv(datafile)
    datain$geneID <- sapply(strsplit(datain$featureID, ":"), `[`, 1) 
    
    genedf <- datain[datain$geneID == gene,]

    gene_annot_df <- merge(genedf, annot_subset, by.x = "featureID", by.y = "ef_id", all.x = TRUE)
    
    num_features <- nrow(gene_annot_df)
    
    print(paste(cell, gene))
    
    if (num_features > 1 && "M" %in% genedf$sex && "F" %in% genedf$sex) {
      ## Mixed-effects meta-analysis - sex (categorical) as moderator, model = FE
      dataout2 <- rma(effect_size, variance, data=gene_annot_df, mods = ~ sex, 
          slab=gene_annot_df$featureID,
          measure="MD", method="FE", )
      
      out2 <- summary(dataout2)
      
      # Extract moderator effect estimates for each level of sex
      moderator_effect_estimates <- coef(dataout2)
      
      # Extract the variance-covariance matrix of the coefficients
      moderator_variances <- diag(vcov(dataout2))
      
      # Calculate weights for the moderator effects as the inverse of the variances
      moderator_weights <- 1 / moderator_variances
      
      # Extract the specific weights for the moderator effects
      intercept_weight <- moderator_weights["intrcpt"]
      sexM_weight <- moderator_weights["sexM"]
      
      # Ensure the levels of sex are consistent with the moderator estimates
      intercept <- moderator_effect_estimates["intrcpt"]
      sexM_effect <- moderator_effect_estimates["sexM"]
      
      # Create a new column in gene_annot_df to store the moderator effects
      gene_annot_df$moderator_effect <- ifelse(gene_annot_df$sex == "M", intercept + sexM_effect, intercept)
      
      # Extract weights and confidence intervals for individual studies
      weights <- weights(dataout2)
      proportional_weights <- weights / sum(weights)
      ci_lower <- gene_annot_df$effect_size - 1.96 * sqrt(gene_annot_df$variance)
      ci_upper <- gene_annot_df$effect_size + 1.96 * sqrt(gene_annot_df$variance)
      
      # Save the individual study data for custom plotting
      plot_data <- data.frame(
        featureID = gene_annot_df$featureID,
        effect_size = gene_annot_df$effect_size,
        variance = gene_annot_df$variance,
        sex = gene_annot_df$sex,
        moderator_effect = gene_annot_df$moderator_effect,
        weight = weights,
        proportional_weight = proportional_weights,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        ef_chr = gene_annot_df$ef_chr,
        ef_start = gene_annot_df$ef_start,
        ef_end = gene_annot_df$ef_end,
        ef_strand = gene_annot_df$ef_strand,
        ef_ir_flag = gene_annot_df$ef_ir_flag
      )
      
      # Output file path for individual study data
      output_file_study <- file.path(outData, paste("plot_data_", gene, "_", cell, "_allFeatures_FE_study.csv", sep = ""))
      
      # Write individual study data to CSV
      write.csv(plot_data, output_file_study, row.names = FALSE)
      
      # Save the overall model estimate details for custom plotting
      overall_estimate <- dataout2$b[1]
      overall_ci_lower <- dataout2$ci.lb[1]
      overall_ci_upper <- dataout2$ci.ub[1]
      overall_se <- dataout2$se[1]
      overall_variance <- overall_se^2
      overall_weight <- 1 / overall_variance
      
      overall_data <- data.frame(
        featureID = "Overall",
        effect_size = overall_estimate,
        variance = overall_variance,
        sex = NA,
        ef_ir_flag = NA,
        moderator_effect = NA,
        intercept_weight = intercept_weight,
        sexM_weight = sexM_weight,
        weight = overall_weight,
        ci_lower = overall_ci_lower,
        ci_upper = overall_ci_upper
      )
      
      # Output file path for overall model data
      output_file_overall <- file.path(outData, paste("plot_data_", gene, "_", cell, "_allFeatures_FE_overall.csv", sep = ""))
      
      # Write overall model data to CSV
      write.csv(overall_data, output_file_overall, row.names = FALSE)
      
    } else if (num_features < 1) {
      print(paste(gene, " in ", cell, "does not have enough features. Cannot estimate parameters"))
    } else if (!("M" %in% genedf$sex)) {
      print(paste(gene, " in ", cell, "does not have any features with sex male"))
    } else if (!("F" %in% genedf$sex)) {
      print(paste(gene, " in ", cell, "does not have any features with sex female"))
    }
  }
}