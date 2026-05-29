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

outPlot <- file.path(workdir, "quantify_t1d_pacbio_transcripts/meta_analysis_DE_pdiffs_plots")
dir.create(outPlot)
outSum  <- file.path(workdir, "quantify_t1d_pacbio_transcripts/meta_analysis_DE_pdiffs_summary")
dir.create(outSum)

#Import annotations of fragments
annotFile=file.path(workdir,"allChr_isoseq_FSM_ISM_NIC_event_analysis_ef.csv")
annotIn=read.csv(annotFile)

annot_subset=annotIn[, c("ef_id", "ef_ir_flag","ef_start")]

t1d_df=read.csv(paste(workdir,"/candidate_gene_lists/T1D_candidate_genes_robertson_2021_ensembl.txt", sep=""))
t1d_genes=t1d_df$ensembl_geneID

library(metafor)
library(broom)

cell_list <- c("CD4", "CD8")   #create character vector

for(cell in cell_list) {
  inputFile <- paste("ES_SV_", cell, "_pdiff_4_ma.csv", sep="")
 
  datafile <- file.path(workdir, "quantify_t1d_pacbio_transcripts", inputFile)
  datain <- read.csv(datafile)
  datain$geneID <- sapply(strsplit(datain$featureID, ":"), `[`, 1)    

      geneCol <- datain$geneID
      for(gene in unique(geneCol)) {
        print(gene)
    
        genedf = datain[datain$geneID == gene,]
        gene_annot_df=merge(genedf, annot_subset, by.x = "featureID", by.y ="ef_id", all.x = TRUE)
        sorted_genedf=gene_annot_df[order(gene_annot_df$ef_start,gene_annot_df$sex),]
        
        num_features <- nrow(genedf)
        
        if (num_features > 1 && "M" %in% genedf$sex && "F" %in% genedf$sex) {
        ## (1) fixed-effects meta-analysis - no moderator
        dataout = rma(effect_size, variance, data=sorted_genedf,
                      slab=sorted_genedf$featureID,
                      measure="MD", method="FE", )
      
        fileSum1 <- file.path(outSum, paste("summary_", gene, "_", cell, "_FE.txt", sep=""))
        out <- summary(dataout)
        capture.output(out, file = fileSum1)
        
        ## (2) mixed-effects meta-analysis - sex (categorical) as moderator, model = FE
        dataout2 = rma(effect_size, variance, data=sorted_genedf, mods = ~ sex, 
                       slab=sorted_genedf$featureID, 
                       measure="MD", method="FE", )
        
        fileSum2 <- file.path(outSum, paste("summary_sexMod_", gene, "_", cell, "_FE.txt", sep=""))
        out2 <- summary(dataout2)
        capture.output(out2, file = fileSum2)
     
        ## (3) mixed-effects meta-analysis - sex (categorical) as moderator, model = REML
        #dataout3 = rma(effect_size, variance, data=sorted_genedf, mods = ~ sex, 
        #               slab=sorted_genedf$featureID,
        #               measure="MD", method="REML", )
        
        #fileSum3 <- file.path(outSum, paste("summary_sexMod_", gene, "_", cell, "_REML.txt", sep=""))
        #out3 <- summary(dataout3)
        #capture.output(out3, file = fileSum3)
        
        ## (4) mixed-effects meta-analysis - cellType (categorical) as moderator, model = FE
        #dataout4 = rma(effect_size, variance, data=sorted_genedf, mods = ~ cellType, 
        #               slab=sorted_genedf$featureID,
        #               measure="MD", method="FE", )
        
        #fileSum4 <- file.path(outSum, paste("summary_cellTypeMod_", gene, "_", cell, "_FE.txt", sep=""))
        #out4 <- summary(dataout4)
        #capture.output(out4, file = fileSum4)
        
        ## (5) mixed-effects meta-analysis - cellType (categorical) as moderator, model = REML
        #dataout5 = rma(effect_size, variance, data=genedf, mods = ~ cellType, 
        #               slab=genedf$featureID,
        #               measure="MD", method="REML", )
        
        #fileSum5 <- file.path(outSum, paste("summary_cellTypeMod_", gene, "_", cell, "_REML.txt", sep=""))
        #out5 <- summary(dataout5)
        #capture.output(out5, file = fileSum5)
        
        
        if (gene %in% t1d_genes) {
        ## (1A) forest plot from fixed effects model no moderator
        outPNG <- paste("forestPlot_", gene, "_", cell, "_FE.png", sep="")
        png(file.path(outPlot, outPNG), width=960, height=640)
        forest(dataout, header=TRUE, slab=sorted_genedf$featureID, order=sex, ilab=sex,
               xlab="For each sex, Case vs Control (MD, FE)")
        dev.off()
 
        ## (2A) forest plot from fixed effects model - sex moderator, FE
        outPNG2 <- paste("forestPlot_sexMod_", gene, "_", cell, "_FE.png", sep="")
        png(file.path(outPlot, outPNG2), width=960, height=640)
        forest(dataout2, header=TRUE, slab=sorted_genedf$featureID, 
               ilab=sex, ilab.xpos=-3.5, order=sex, 
               xlab="For each sex, Case vs Control (MD, FE)")    
        dev.off()
       
        ## (3A) forest plot from mixed effects model - sex moderator, REML
        #outPNG3 <- paste("forestPlot_sexMod_", gene, "_", cell, "_REML.png", sep="")
        #png(file.path(outPlot, outPNG3), width=960, height=640)
        #forest(dataout3, header=TRUE, slab=genedf$featureID, 
        #       ilab=sex, ilab.xpos=-3.5, order=sex, 
        #       xlab="For each sex, Case vs Control (MD, REML)")    
        #dev.off()
        
        
        ## (4A) forest plot from mixed effects model - sex moderator, REML
        #outPNG4 <- paste("forestPlot_cellTypeMod_", gene, "_", cell, "_FE.png", sep="")
        #png(file.path(outPlot, outPNG4), width=960, height=640)
        #forest(dataout4, header=TRUE, slab=genedf$featureID, 
        #       ilab=sex, ilab.xpos=-3.5, order=cellType, 
        #       xlab="For each cell type, Case vs Control (MD, FE)")    
        #dev.off()
        
        ## (5A) forest plot from mixed effects model - sex moderator, REML
        #outPNG5 <- paste("forestPlot_cellTypeMod_", gene, "_", cell, "_REML.png", sep="")
        #png(file.path(outPlot, outPNG5), width=960, height=640)
        #forest(dataout5, header=TRUE, slab=genedf$featureID, 
        #       ilab=sex, ilab.xpos=-3.5, order=cellType, 
        #       xlab="For each cell type, Case vs Control (MD, REML)")    
        #dev.off()
        } else if (num_features < 1) {
          print(paste(gene," in ",cell, "does not have enough features. Cannot estimate parameters"))
        } else if (!("M" %in% genedf$sex)) {
          print(paste(gene," in ",cell, "does not have any features with sex male"))
        } else if (!("F" %in% genedf$sex)) {
          print(paste(gene," in ",cell, "does not have any features with sex female"))
        }
        
        }
    }
}
