#===============================================================================
#
#          FILE:  wiggleplots_example.R
# 
#         USAGE:  R CMD wiggleplots_example.R 
# 
#   DESCRIPTION:  This is an example script to create wiggle plots using R.
# 
#  REQUIREMENTS:  If you are going to output as an SVG you need to install
#  "RSvgDevice" by running from inside R:
#
#                 > install.packages("RSvgDevice")
#
#        AUTHOR:  Justin Fear
# 
#===============================================================================
library(ggplot2)

## Get Command Line Arguments
    args <- commandArgs(TRUE);
    my.gene <- args[2]
    my.sex <- args[3]
    bga <- as.numeric(args[4])
    bgb <- as.numeric(args[5])
    bgc <- as.numeric(args[6])

## Get Environmental Variables
    MCLAB <- Sys.getenv("MCLAB");

## load wiggle plot functions
    source(file=paste(MCLAB,"/arbeitman_fru_network/r_programs/wiggleplot_functions.R",sep=""));

## set working direcotry to where files are
    setwd(args[1])

## Get list of files in directory
    files <- list.files(pattern='*.csv')

## Create Wiggle Tracks
    ### Read in track information
        frua.file <- files[1];
        frub.file <- files[2];
        fruc.file <- files[3];

        frua <- read.csv(frua.file,header=FALSE);
        frub <- read.csv(frub.file,header=FALSE);
        fruc <- read.csv(fruc.file,header=FALSE);

        colnames(frua) <- c("chrom", "coord", "count")
        colnames(frub) <- c("chrom", "coord", "count")
        colnames(fruc) <- c("chrom", "coord", "count")


    ### Remove Background
        frua$count <- frua$count - bga
        frua$count[frua$count < 0] <- 0

        frub$count <- frub$count - bgb
        frub$count[frub$count < 0] <- 0

        fruc$count <- fruc$count - bgc
        fruc$count[fruc$count < 0] <- 0

    ### Combine all datasets into a single data frame
        frua$FruM <- 'A'
        frub$FruM <- 'B'
        fruc$FruM <- 'C'

        fru.all <- rbind(frua,frub,fruc)

    ### Create a ggplot object for Wiggle plots
        wiggle <- ggplot() +  geom_ribbon(data=fru.all, aes(x=coord,y=count,ymin=0,ymax=count,fill=FruM,alpha=0.2)) + scale_fill_manual(values=c("red","orange","blue"))+ ggtitle(paste("AH",my.sex,my.gene,sep=" ")) +
        theme(panel.background=element_blank(),panel.grid=element_blank(),plot.background=element_blank())


## Create Gene Model Tracks
    ### Read in Exon information
        exon2symbol.file <- paste(MCLAB,"/useful_dmel_data/flybase530/exon2symbol.csv",sep="")
        exon2symbol  <- read.csv(exon2symbol.file, stringsAsFactors=F, header=T,colClasses=c(rep("NULL",3),"integer","integer","NULL","character",rep("NULL",3),"character","NULL"))

    ### Create A subset dataset based on my current gene of interest
        exon2symbol.subset <- subset(exon2symbol,exon2symbol$symbol == my.gene)
        exon2symbol.subset$y <- 0

    ### Create A uniq list of FBtr and sort by start site
        uniq.fbtr <- unique(sort(exon2symbol.subset$FBtr))
        df = data.frame(FBtr = uniq.fbtr, min.start = rep(0, length(uniq.fbtr)))

        j <- 1
        while(j <= length(uniq.fbtr)){
            
            min.start <- min(exon2symbol.subset[exon2symbol.subset$FBtr == uniq.fbtr[j],]$start_exon);
            max.end <- max(exon2symbol.subset[exon2symbol.subset$FBtr == uniq.fbtr[j],]$end_exon);
            df[df$FBtr == uniq.fbtr[j],]$min.start <- min.start

            j <- j + 1
        }

        df.sort <- df[order(df$min.start),]

    ### Create A list of Y axis points so that Gene models will be plotted below each other
        j <- 1
        y <- -5
        while(j <= length(uniq.fbtr)){
            exon2symbol.subset[exon2symbol.subset$FBtr == df.sort[j,1],]$y <- y;
            y <- y - 5
            j <- j + 1
        }

## Plot Wiggle Plots and Gene models
    svg(file=paste(my.gene,"_AH_",my.sex,"_wig_v2.svg",sep=""), width=28, height=7)
    plot(wiggle)    
    dev.off()

    svg(file=paste(my.gene,"_AH_",my.sex,"_wig_models_v2.svg",sep=""), width=28, height=7)
    plot(wiggle + geom_rect(data=exon2symbol.subset, aes(xmin=start_exon, xmax=end_exon, ymin = y, ymax=y+3),fill="black"))
    dev.off()

