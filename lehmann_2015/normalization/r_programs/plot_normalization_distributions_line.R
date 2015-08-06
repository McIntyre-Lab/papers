library(ggplot2)
library(grid)
library(gridExtra)

mclab <- Sys.getenv("MCLAB")

# Parse command line arguments 1 = filename, 2 = mated or virgin
args <- commandArgs(TRUE)
mname <- args[1]
vname <- args[2]

# Read in data
mated <- read.csv(mname)
virgin <- read.csv(vname)


# Create Plot of LOG UQ APN
pm <- ggplot(mated,aes(sample_id,log_uq_apn)) + geom_boxplot(outlier.shape = NA) + theme(axis.text.x = element_text(angle=90,vjust=0.5)) + ylim(0,10) + ggtitle("Boxplot of Log UQ APN by Genotype\nMated")
pv <- ggplot(virgin,aes(sample_id,log_uq_apn)) + geom_boxplot(outlier.shape = NA) + theme(axis.text.x = element_text(angle=90,vjust=0.5)) + ylim(0,10) + ggtitle("Boxplot of Log UQ APN by Genotype\nVirgin")


pdf(paste0(mclab,"/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/OE_normalization/report/Mated_Virgin_log_uq_apn_boxplot_line.pdf"),width=11,height=8.5)
    grid.arrange(pm, pv, ncol=1)
dev.off()
