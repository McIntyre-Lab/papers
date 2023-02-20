#!/usr/bin/env Rscript --vanilla

# Generate bar graph of histone detection agreement (Cohen's Kappa)

# Input agreement CSV file must contain column names "featureType", "K4_kappa", and K27_kappa"

library(readr)
library(ggplot2)

args = commandArgs(trailingOnly = TRUE)

if(length(args)!=5) {
  stop("3 arguments required in this order: /path/to/species1_agreement.csv /path/to/species2_agreement.csv /path/to/output species1 species2", call.=FALSE)
} else if (length(args)==5){
  if(endsWith(args[1],".csv")){
    inFILE1 = args[1]
  } else {
    stop("First argument must be /path/to/species1_agreement.csv", call.=FALSE)
  }
  if(endsWith(args[2],".csv")){
    inFILE2 = args[2]
  } else {
    stop("First argument must be /path/to/species2_agreement.csv", call.=FALSE)
  }
  outFILE = args[3]
  name1 = args[4]
  name2 = args[5]
}

if(name1 == "mel"){
  fullName1 = "D. melanogaster"
} else if(name1 == "sim"){
  fullName1 = "D. simulans"
}
if(name2 == "mel"){
  fullName2 = "D. melanogaster"
} else if(name2 == "sim"){
  fullName2 = "D. simulans"
}

kappaDF1 <- read_csv(inFILE1)
kappaDF2 <- read_csv(inFILE2)

## Subset kappa values for male vs. female H3K4me3 and H3K27me2me3
kappaDF1 <- kappaDF1[,c("featureType","K4_kappa","K27_kappa")]
kappaDF2 <- kappaDF2[,c("featureType","K4_kappa","K27_kappa")]

## Reshape dataframes for each species to set histone mark groups
kappaDF1 <- reshape(as.data.frame(kappaDF1),idvar = "featureType", varying = c("K4_kappa","K27_kappa"),
        v.names = "kappa", direction = "long",
        times = c(paste(fullName1," Male vs. Female H3K4me3", sep=""),
                  paste(fullName1," Male vs. Female H3K27me2me3", sep="")),
        new.row.names = 1:1000)
names(kappaDF1)[names(kappaDF1) == "time"] <- "group"
kappaDF2 <- reshape(as.data.frame(kappaDF2),idvar = "featureType", varying = c("K4_kappa","K27_kappa"),
        v.names = "kappa", direction = "long",
        times = c(paste(fullName2," Male vs. Female H3K4me3", sep=""),
                  paste(fullName2," Male vs. Female H3K27me2me3", sep="")),
        new.row.names = 1:1000)
names(kappaDF2)[names(kappaDF2) == "time"] <- "group"

## Combine species
combKappa <- rbind(kappaDF1,kappaDF2)

## Drop fusions and TSS1kbWindows
combKappa <- combKappa[!combKappa$featureType %in% c("fusion","TSS1kbWindow"),]

## Rename fragment to exonic, intron to intronic, TSS300bpwindow to TSS
combKappa$featureType[combKappa$featureType == 'fragment'] <- 'Exonic'
combKappa$featureType[combKappa$featureType == 'intron'] <- 'Intronic'
combKappa$featureType[combKappa$featureType == 'TSS300bpWindow'] <- 'TSS'
combKappa$featureType[combKappa$featureType == 'intergenic'] <- 'Intergenic'

## Set feature type order
combKappa$featureType <- factor(combKappa$featureType, levels=c("TSS","5UTR","Exonic","Intronic","3UTR","Intergenic"))

## Plot
for(fileType in c("png","tiff","eps")){
  tempPlot = ggplot(combKappa, aes(x = featureType, y = kappa, fill = group)) +
  geom_col(position = position_dodge(), width = 0.5) +
  labs(x="", y="", fill="Kappa Value") +
  scale_fill_manual(values = c("#4472C4","#A5A5A5","#5B9BD5","#264478")) +
  theme_minimal() + theme(legend.position="bottom", legend.direction = "vertical")
  ggsave(plot=tempPlot,filename=paste(outFILE,".",fileType,sep=""),device=fileType,
         height=4, width=8)
}
