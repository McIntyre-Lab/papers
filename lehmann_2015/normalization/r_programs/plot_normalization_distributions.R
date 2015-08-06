library(ggplot2)
mclab <- Sys.getenv("MCLAB")

# MATED
    mydata <- read.csv('/home/jfear/tmp/mated.csv')
    mydata$flag_raleigh <- factor(mydata$flag_raleigh)
    mydata$line <- factor(mydata$line)
    mydata$rep <- factor(mydata$rep)
    mydata$sample <- factor(mydata$sample)
    mydata$log_rpkm <- log(mydata$rpkm,2)
    mydata$log_apn <- log(mydata$apn,2)

    # APN
        png(paste(mclab,"/cegs_sergey/reports/line_normalization/mated_apn_density.png",sep=""),width=800,height=800)
        ggplot(mydata,aes(apn,color=flag_raleigh,group=sample)) + geom_density() + xlim(0,100) + ylim(0,.4) + ggtitle("Density Plot of APN by Genotype\nMated")
        dev.off()

        png(paste(mclab,"/cegs_sergey/reports/line_normalization/mated_apn_boxplot.png",sep=""),width=2500,heigh=800)
        ggplot(mydata,aes(sample,apn,fill=flag_raleigh)) + geom_boxplot(outlier.shape = NA) + theme(axis.text.x = element_text(angle=90,vjust=.5)) + ylim(0,50) + ggtitle("Boxplot of APN by Genotype\nMated")
        dev.off()

    # RPKM 
        png(paste(mclab,"/cegs_sergey/reports/line_normalization/mated_rpkm_boxplot.png",sep=""),width=2500,heigh=800)
        ggplot(mydata,aes(sample,rpkm,fill=flag_raleigh)) + geom_boxplot(outlier.shape = NA) + theme(axis.text.x = element_text(angle=90,vjust=.5)) + ylim(0,50) + ggtitle("Boxplot of RPKM by Genotype\nMated")
        dev.off()

    # Log RPKM 
        png(paste(mclab,"/cegs_sergey/reports/line_normalization/mated_log_rpkm_boxplot.png",sep=""),width=2500,heigh=800)
        ggplot(mydata,aes(sample,log_rpkm,fill=flag_raleigh)) + geom_boxplot(outlier.shape = NA) + theme(axis.text.x = element_text(angle=90,vjust=.5)) + ylim(0,10) + ggtitle("Boxplot of Log RPKM by Genotype\nMated")
        dev.off()

    # Log APN
        png(paste(mclab,"/cegs_sergey/reports/line_normalization/mated_log_apn_boxplot.png",sep=""),width=2500,heigh=800)
        ggplot(mydata,aes(sample,log_apn,fill=flag_raleigh)) + geom_boxplot(outlier.shape = NA) + theme(axis.text.x = element_text(angle=90,vjust=.5)) + ylim(0,10) + ggtitle("Boxplot of Log APN by Genotype\nMated")
        dev.off()

    # LOG UQ APN
        png(paste(mclab,"/cegs_sergey/reports/line_normalization/mated_log_uq_apn_density.png",sep=""),width=2500,height=800)
        ggplot(mydata,aes(log_uq_apn,color=flag_raleigh,group=sample)) + geom_density() + xlim(0,7.5) + ylim(0,.4) + ggtitle("Density Plot of Log UQ APN by Genotype\nMated")
        dev.off()

        png(paste(mclab,"/cegs_sergey/reports/line_normalization/mated_log_uq_apn_boxplot.png",sep=""),width=2500,heigh=800)
        ggplot(mydata,aes(sample,log_uq_apn,fill=flag_raleigh)) + geom_boxplot(outlier.shape = NA) + theme(axis.text.x = element_text(angle=90,vjust=0.5)) + ylim(0,10) + ggtitle("Boxplot of Log UQ APN by Genotype\nMated")
        dev.off()

    # UQ FF 
        uq_ff <- mydata[c("line", "rep", "uq_ff", "flag_raleigh")]
        uniq_ff <- unique(uq_ff)
        png(paste(mclab,"/cegs_sergey/reports/line_normalization/mated_uq_ff_boxplot.png",sep=""),width=2500,heigh=800)
        ggplot(uniq_ff,aes(line,uq_ff)) + geom_boxplot(aes(fill=flag_raleigh)) + theme(axis.text.x = element_text(angle=90,vjust=0.5)) + ggtitle("Boxplot of Log UQ Fudge Factor by Genotype\nMated") + geom_point(aes(color=rep))
        dev.off()

    # Plot Centered data
        png(paste(mclab,"/cegs_sergey/reports/line_normalization/mated_mean_center_boxplot.png",sep=""),width=2500,heigh=800)
        ggplot(mydata,aes(sample,mean_log_uq_center)) + geom_boxplot(aes(fill=flag_raleigh),outlier.shape=NA) + theme(axis.text.x = element_text(angle=90,vjust=0.5)) + ylim(0,10) + ggtitle("Boxplot of Mean Center Log UQ APN\nMated")
        dev.off()

        png(paste(mclab,"/cegs_sergey/reports/line_normalization/mated_median_center_boxplot.png",sep=""),width=2500,heigh=800)
        ggplot(mydata,aes(sample,median_log_uq_center)) + geom_boxplot(aes(fill=flag_raleigh),outlier.shape=NA) + theme(axis.text.x = element_text(angle=90,vjust=0.5)) + ylim(0,10) + ggtitle("Boxplot of Median Center Log UQ APN\nMated")
        dev.off()

        png(paste(mclab,"/cegs_sergey/reports/line_normalization/mated_uq_center_boxplot.png",sep=""),width=2500,heigh=800)
        ggplot(mydata,aes(sample,uq_log_uq_center)) + geom_boxplot(aes(fill=flag_raleigh),outlier.shape=NA) + theme(axis.text.x = element_text(angle=90,vjust=0.5)) + ylim(0,10) + ggtitle("Boxplot of UQ Center Log UQ APN\nMated")
        dev.off()

# VIRGIN
    mydata <- read.csv('/home/jfear/tmp/virgin.csv')
    mydata$flag_raleigh <- factor(mydata$flag_raleigh)
    mydata$line <- factor(mydata$line)
    mydata$rep <- factor(mydata$rep)
    mydata$sample <- factor(mydata$sample)
    mydata$log_rpkm <- log(mydata$rpkm,2)
    mydata$log_apn <- log(mydata$apn,2)

    # APN
        png(paste(mclab,"/cegs_sergey/reports/line_normalization/virgin_apn_density.png",sep=""),width=800,height=800)
        ggplot(mydata,aes(apn,color=flag_raleigh,group=sample)) + geom_density() + xlim(0,100) + ylim(0,.4) + ggtitle("Density Plot of APN by Genotype\nVirgin")
        dev.off()

        png(paste(mclab,"/cegs_sergey/reports/line_normalization/virgin_apn_boxplot.png",sep=""),width=2500,heigh=800)
        ggplot(mydata,aes(sample,apn,fill=flag_raleigh)) + geom_boxplot(outlier.shape = NA) + theme(axis.text.x = element_text(angle=90,vjust=.5)) + ylim(0,50) + ggtitle("Boxplot of APN by Genotype\nVirgin")
        dev.off()

    # RPKM 
        png(paste(mclab,"/cegs_sergey/reports/line_normalization/virgin_rpkm_boxplot.png",sep=""),width=2500,heigh=800)
        ggplot(mydata,aes(sample,rpkm,fill=flag_raleigh)) + geom_boxplot(outlier.shape = NA) + theme(axis.text.x = element_text(angle=90,vjust=.5)) + ylim(0,50) + ggtitle("Boxplot of RPKM by Genotype\nVirgin")
        dev.off()

    # Log RPKM 
        png(paste(mclab,"/cegs_sergey/reports/line_normalization/virgin_log_rpkm_boxplot.png",sep=""),width=2500,heigh=800)
        ggplot(mydata,aes(sample,log_rpkm,fill=flag_raleigh)) + geom_boxplot(outlier.shape = NA) + theme(axis.text.x = element_text(angle=90,vjust=.5)) + ylim(0,10) + ggtitle("Boxplot of Log RPKM by Genotype\nVirgin")
        dev.off()

    # Log APN
        png(paste(mclab,"/cegs_sergey/reports/line_normalization/virgin_log_apn_boxplot.png",sep=""),width=2500,heigh=800)
        ggplot(mydata,aes(sample,log_apn,fill=flag_raleigh)) + geom_boxplot(outlier.shape = NA) + theme(axis.text.x = element_text(angle=90,vjust=.5)) + ylim(0,10) + ggtitle("Boxplot of Log APN by Genotype\nVirgin")
        dev.off()

    # LOG UQ APN
        png(paste(mclab,"/cegs_sergey/reports/line_normalization/virgin_log_uq_apn_density.png",sep=""),width=800,height=800)
        ggplot(mydata,aes(log_uq_apn,color=flag_raleigh,group=sample)) + geom_density() + xlim(0,100) + ylim(0,.4) + ggtitle("Density Plot of Log UQ APN by Genotype\nVirgin")
        dev.off()

        png(paste(mclab,"/cegs_sergey/reports/line_normalization/virgin_log_uq_apn_boxplot.png",sep=""),width=2500,heigh=800)
        ggplot(mydata,aes(sample,log_uq_apn,fill=flag_raleigh)) + geom_boxplot(outlier.shape=NA) + theme(axis.text.x = element_text(angle=90,vjust=0.5)) + ylim(0,15) + ggtitle("Boxplot of Log UQ APN by Genotype\nVirgin")
        dev.off()

    # UQ FF 
        uq_ff <- mydata[c("line", "rep", "uq_ff", "flag_raleigh")]
        uniq_ff <- unique(uq_ff)
        png(paste(mclab,"/cegs_sergey/reports/line_normalization/virgin_uq_ff_boxplot.png",sep=""),width=2500,heigh=800)
        ggplot(uniq_ff,aes(line,uq_ff)) + geom_boxplot(aes(fill=flag_raleigh)) + theme(axis.text.x = element_text(angle=90,vjust=0.5)) + ggtitle("Boxplot of Log UQ Fudge Factor by Genotype\nVirgin") + geom_point(aes(color=rep))
        dev.off()

    # Plot Centered data
        png(paste(mclab,"/cegs_sergey/reports/line_normalization/virgin_mean_center_boxplot.png",sep=""),width=2500,heigh=800)
        ggplot(mydata,aes(sample,mean_log_uq_center)) + geom_boxplot(aes(fill=flag_raleigh),outlier.shape=NA) + theme(axis.text.x = element_text(angle=90, vjust=0.5)) + ylim(0,10) + ggtitle("Boxplot of Mean Center Log UQ APN\nVirgin")
        dev.off()

        png(paste(mclab,"/cegs_sergey/reports/line_normalization/virgin_median_center_boxplot.png",sep=""),width=2500,heigh=800)
        ggplot(mydata,aes(sample,median_log_uq_center)) + geom_boxplot(aes(fill=flag_raleigh),outlier.shape=NA) + theme(axis.text.x = element_text(angle=90, vjust=0.5)) + ylim(0,10) + ggtitle("Boxplot of Median Center Log UQ APN\nVirgin")
        dev.off()

        png(paste(mclab,"/cegs_sergey/reports/line_normalization/virgin_uq_center_boxplot.png",sep=""),width=2500,heigh=800)
        ggplot(mydata,aes(sample,uq_log_uq_center)) + geom_boxplot(aes(fill=flag_raleigh),outlier.shape=NA) + theme(axis.text.x = element_text(angle=90, vjust=0.5)) + ylim(0,10) + ggtitle("Boxplot of UQ Center Log UQ APN\nVirgin")
        dev.off()
