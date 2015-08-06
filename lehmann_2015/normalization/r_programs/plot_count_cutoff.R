library(ggplot2)
library(grid)
library(gridExtra)

mclab <- Sys.getenv("MCLAB")

# Parse command line arguments 1 = filename, 2 = mated or virgin
args <- commandArgs(TRUE)
sname <- args[1]

# Read in data
sums <- read.csv(sname)

# Scatter plot looking at relationship between APN > 0 and APN >= 5
pdf(paste0(mclab,"/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/OE_normalization/report/exon_expression_counts.pdf"),width=8.5,height=8.5)
    cuts <- data.frame(Thresh = c("2k","5k","7.5k"), vals = c(2000,5000,7500))
    p <- ggplot(sums,aes(x=cnt_apn_gt_0,y=cnt_apn_gt_5)) + geom_point() + ylab("Number exonic regions with APN >5") + xlab("Number exonic regions with APN >0")
    p + geom_vline(xintercept=29300) + geom_hline(data=cuts,aes(yintercept=vals,linetype=Thresh),show_guide=TRUE) + scale_color_manual(values=c("blue","red")) + geom_text(data=NULL, x=25000, y=18000, label="Vertical Cutoff 29,300", color="Purple")
dev.off()

