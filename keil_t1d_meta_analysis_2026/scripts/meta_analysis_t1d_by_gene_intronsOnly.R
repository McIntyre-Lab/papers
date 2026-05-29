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

outPlot <- file.path(workdir, "quantify_t1d_pacbio_transcripts/meta_analysis_DE_pdiffs_plots_intronsOnly")
dir.create(outPlot)
outSum  <- file.path(workdir, "quantify_t1d_pacbio_transcripts/meta_analysis_DE_pdiffs_summary_intronsOnly")
dir.create(outSum)

library(metafor)
library(broom)

#Import annotations of fragments
annotFile=file.path(workdir,"allChr_isoseq_FSM_ISM_NIC_event_analysis_ef.csv")
annotIn=read.csv(annotFile)

annot_subset=annotIn[, c("ef_id", "ef_ir_flag","ef_start")]

t1d_df=read.csv(paste(workdir,"/candidate_gene_lists/T1D_candidate_genes_robertson_2021_ensembl.txt", sep=""))
t1d_genes=t1d_df$ensembl_geneID

cell_list <- c("CD4", "CD8")   #create character vector

for(cell in cell_list) {
  inputFile <- paste("ES_SV_", cell, "_pdiff_4_ma.csv", sep="")
  
  datafile <- file.path(workdir, "quantify_t1d_pacbio_transcripts", inputFile)
  datain <- read.csv(datafile)
  datain$geneID <- sapply(strsplit(datain$featureID, ":"), `[`, 1) 
  
  geneCol <- datain$geneID
  for(gene in unique(geneCol)) {
    
    genedf = datain[datain$geneID == gene,]
    
    gene_annot_df=merge(genedf, annot_subset, by.x = "featureID", by.y ="ef_id", all.x = TRUE)
    
    introndf= subset(gene_annot_df, ef_ir_flag == 1)
    
    num_features <- nrow(introndf)
    
    print(paste(cell,gene))
    
    if (num_features > 1 && "M" %in% genedf$sex && "F" %in% genedf$sex) {
      ## (1) fixed-effects meta-analysis - no moderator
      dataout = rma(effect_size, variance, data=introndf,
                    slab=introndf$featureID,
                    measure="MD", method="FE", )
      
      fileSum1 <- file.path(outSum, paste("summary_", gene, "_", cell, "_intronsOnly_FE.txt", sep=""))
      out <- summary(dataout)
      capture.output(out, file = fileSum1)
      
      ## (2) mixed-effects meta-analysis - sex (categorical) as moderator, model = FE
      dataout2 = rma(effect_size, variance, data=introndf, mods = ~ sex, 
                     slab=introndf$featureID,
                     measure="MD", method="FE", )
      
      fileSum2 <- file.path(outSum, paste("summary_sexMod_", gene, "_", cell, "_intronsOnly_FE.txt", sep=""))
      out2 <- summary(dataout2)
      capture.output(out2, file = fileSum2)
      
      
      ## (3) mixed-effects meta-analysis - sex (categorical) as moderator, model = REML
      #dataout3 = rma(effect_size, variance, data=introndf, mods = ~ sex, 
      #               slab=introndf$featureID,
      #               measure="MD", method="REML", )
      #
      #fileSum3 <- file.path(outSum, paste("summary_sexMod_", gene, "_", cell, "_REML.txt", sep=""))
      #out3 <- summary(dataout3)
      #capture.output(out3, file = fileSum3)
      
      if (gene %in% t1d_genes) {
        
        ## (1A) forest plot from fixed effects model no moderator
        outPNG <- paste("forestPlot_", gene, "_", cell, "_intronsOnly_FE.png", sep="")
        png(file.path(outPlot, outPNG), width=960, height=640)
        forest(dataout, header=TRUE, slab=introndf$featureID, order=sex, ilab=sex,
               xlab="For each sex, Case vs Control (MD, FE)")
        dev.off()
        
        ## (2A) forest plot from fixed effects model - sex moderator, FE
        outPNG2 <- paste("forestPlot_sexMod_", gene, "_", cell, "_intronsOnly_FE.png", sep="")
        png(file.path(outPlot, outPNG2), width=960, height=640)
        forest(dataout2, header=TRUE, slab=introndf$featureID, 
               ilab=sex, ilab.xpos=-3.5, order=sex, 
               xlab="For each sex, Case vs Control (MD, FE)")    
        dev.off()
        
        ## (3A) forest plot from mixed effects model - sex moderator, REML
        #outPNG3 <- paste("forestPlot_sexMod_", gene, "_", cell, "_REML.png", sep="")
        #png(file.path(outPlot, outPNG3), width=960, height=640)
        #forest(dataout3, header=TRUE, slab=introndf$featureID, 
        #       ilab=sex, ilab.xpos=-3.5, order=sex, 
        #       xlab="For each sex, Case vs Control (MD, REML)")    
        #dev.off()
      }
      
    } else if (num_features < 1) {
      print(paste(gene," in ",cell, "does not have enough features. Cannot estimate parameters"))
    } else if (!("M" %in% genedf$sex)) {
      print(paste(gene," in ",cell, "does not have any features with sex male"))
    } else if (!("F" %in% genedf$sex)) {
      print(paste(gene," in ",cell, "does not have any features with sex female"))
    }
  }
}
