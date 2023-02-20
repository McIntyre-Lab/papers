#!/usr/bin/env Rscript --vanilla

# Generate stacked bar graph of histone detection groups
# (Green Figure)

# Make plots with and without "None" (aka no detection of any marks in either sex)

library(readr)
library(ggplot2)
library(scales)

args = commandArgs(trailingOnly = TRUE)

if(length(args)!=2) {
  stop("3 arguments required in this order: /path/to/counts.csv /path/to/output", call.=FALSE)
} else if (length(args)==2){
  if(endsWith(args[1],".csv")){
    inFILE = args[1]
  } else {
    stop("First argument must be /path/to/counts.csv", call.=FALSE)
  }
  outFILE = args[2]
}

counts <- read_csv(inFILE)

## Set detection flag string
counts$flag_string <- paste(counts$fK4, counts$fK27, counts$mK4, counts$mK27, sep="")

## Drop intergenic regions
counts <- counts[counts$featureType != "intergenic",]

## Rename fragment to exonic and intron to intronic
counts$featureType[counts$featureType == 'fragment'] <- 'Exonic'
counts$featureType[counts$featureType == 'intron'] <- 'Intronic'

## Set group names based on detection flag strings
counts$group <- ifelse(counts$flag_string=="0001", "H3K27me2me3 Male Only",
                       ifelse(counts$flag_string=="0010", "H3K4me3 Male Only",
                              ifelse(counts$flag_string=="0100", "H3K27me2me3 Female Only",
                                     ifelse(counts$flag_string=="0101", "H3K27me2me3 Both Sexes",
                                            ifelse(counts$flag_string=="1000","H3K4me3 Female Only",
                                                   ifelse(counts$flag_string=="1010","H3K4me3 Both Sexes",
                                                          ifelse(counts$flag_string=="0000","None",
                                                                 ifelse(counts$flag_string %in% c("0011","0110","0111","1001","1011","1100","1101","1110","1111"),"Mixed","OOPS"))))))))
## Reorder by feature type and group
counts$group <- factor(counts$group, levels=c("None","Mixed","H3K27me2me3 Male Only","H3K27me2me3 Female Only","H3K27me2me3 Both Sexes","H3K4me3 Male Only","H3K4me3 Female Only","H3K4me3 Both Sexes"))
counts$featureType <- factor(counts$featureType, levels=c("TSS","5UTR","Exonic","Intronic","3UTR"))

## Plot with None inculded
for(fileType in c("png","tiff","eps")){
  tempPlot = ggplot(counts, aes(x = featureType, y = freq, fill = group)) +
    geom_col(position = "fill", width = 0.5) +
    labs(x="", y="", fill="Detected Sex Bias") +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = c("black","gray90","gray75","gray50","gray40","#00B050","#385723","#C5E0B4")) +
    theme_minimal()
  ggsave(plot=tempPlot,filename=paste(outFILE,"_with_None.",fileType,sep=""),device=fileType,
         height=4, width=8)
}

## Drop counts for features undetected in all samples
counts <- counts[counts$flag_string != "0000",]

## Plot with None inculded
for(fileType in c("png","tiff","eps")){
  tempPlot = ggplot(counts, aes(x = featureType, y = freq, fill = group)) +
    geom_col(position = "fill", width = 0.5) +
    labs(x="", y="", fill="Detected Sex Bias") +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = c("gray90","gray75","gray50","gray40","#00B050","#385723","#C5E0B4")) +
    theme_minimal()
  ggsave(plot=tempPlot,filename=paste(outFILE,"_detected_only.",fileType,sep=""),device=fileType,
         height=4, width=8)
}