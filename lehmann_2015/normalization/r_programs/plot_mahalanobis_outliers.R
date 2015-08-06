library(ggplot2)
mclab <- Sys.getenv("MCLAB")

# Parse command line arguments

args <- commandArgs(TRUE)
fname <- args[1]
fname <- '/tmp/outlier.csv'

# Read data
mydata <- read.csv(fname)

png('/home/jfear/mclab/cegs_sergey/reports/mahalanobis_outlier.png')
qplot(factor(flag_mahalanobis_outlier), uq_log_uq_center, data=mydata, geom="boxplot", xlab="Flag Outlier", main="Distirubtion of expression for\nMahalanobis Distance outliers")
dev.off()
