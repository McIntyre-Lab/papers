# sript to import and plot the distribution of the number of duplicates.

plotfun <- function(fname){
    odir <- "/home/jfear/mclab/cegs_sergey/reports_internal/duplicate_distributions/"
    name <- sub("\\.txt","",basename(fname))
    mydata <- read.table(fname,header=FALSE,colClasses=c('integer','NULL'))

    png(paste(odir,name,".png",sep=""))
    par(mar=c(5,5,5,5))
    plot(table(mydata$V1),type="h",main=fname,xlab='number of duplicates', ylab='frequency',cex.lab=2)
    dev.off()
}

setwd("/home/jfear/mclab/cegs_sergey/qc_original_data/alison_duplicates/files/")
design <- list.files()

sapply(design,plotfun)
