library(ggplot2)
mclab <- Sys.getenv("MCLAB")

# Parse command line arguments 1 = filename, 2 = mated or virgin
args <- commandArgs(TRUE)
fname <- args[1]
tname <- args[2]

# Read in data
mydata <- read.csv(fname)

# Create Plot of LOG UQ APN
png(paste0(mclab,"/cegs_sergey/reports/line_normalization/",tname,"_log_uq_apn_boxplot_line.png"),width=2500,heigh=800)
ggplot(mydata,aes(sample_id,log_uq_apn)) + geom_boxplot(outlier.shape = NA) + theme(axis.text.x = element_text(angle=90,vjust=0.5)) + ylim(0,10) + ggtitle(paste0("Boxplot of Log UQ APN by Genotype\n",tname))
dev.off()

