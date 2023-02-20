#!/usr/bin/env Rscript --vanilla

# Plot FRiP values of features (2 species together)
# Requires 5 arguments: /path/to/species1_summary_frip.csv species2_summary_frip.csv /path/to/output.png species1 species2

# FRiP summary files must contain at least "sampleID" and "frip" columns

library(readr)
library(ggplot2)

args = commandArgs(trailingOnly = TRUE)

if(length(args)!=5) {
  stop("3 arguments required in this order: /path/to/species1_summary_frip.csv species2_summary_frip.csv /path/to/output.png species1 species2", call.=FALSE)
} else if (length(args)==5){
  if(endsWith(args[1],".csv")){
    inFILE1 = args[1]
  } else {
    stop("First argument must be /path/to/species1_summary_frip.csv", call.=FALSE)
  }
  if(endsWith(args[2],".csv")){
    inFILE2 = args[2]
  } else {
    stop("Second argument must be /path/to/species2_summary_frip.csv", call.=FALSE)
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

frip1 <- read_csv(inFILE1)
frip2 <- read_csv(inFILE2)

## Set sampleType (including species)
frip1$sampleType <- ifelse(grepl("I_",frip1$sampleID),paste(fullName1,"IgG Control"),ifelse(grepl("K4_",frip1$sampleID),paste(fullName1,"H3K4me3"),paste(fullName1,"H3K27me2me3")))
frip2$sampleType <- ifelse(grepl("I_",frip2$sampleID),paste(fullName2,"IgG Control"),ifelse(grepl("K4_",frip2$sampleID),paste(fullName2,"H3K4me3"),paste(fullName2,"H3K27me2me3")))

## Combine species data frames
fripFull <- rbind(frip1, frip2)

## Rename TSS300bpWindow to just TSS
fripFull$feature_type[fripFull$feature_type == 'TSS300bpWindow'] <- 'TSS'
#fripFull$feature_type[fripFull$feature_type == 'TSS300bpWindowPB'] <- 'TSS'

## Rename fragments to just exons
fripFull$feature_type[fripFull$feature_type == 'fragments'] <- 'exons'
#fripFull$feature_type[fripFull$feature_type == 'fragmentsPB'] <- 'exons'

## Rename fragments to just exons
fripFull$feature_type[fripFull$feature_type == 'fusions'] <- 'exons(fusion)'

## Reorder feature types for graph
fripFull$feature_type_reorder = factor(fripFull$feature_type, levels=c('TSS','5UTR','3UTR','exons','introns','intergenic','TSS1kbWindow','exons(fusion)'))

## Reorder sample types for graph
fripFull$sampleType_reorder = factor(fripFull$sampleType, levels=c("D. melanogaster H3K27me2me3","D. simulans H3K27me2me3","D. melanogaster H3K4me3","D. simulans H3K4me3","D. melanogaster IgG Control","D. simulans IgG Control"))

## Remove unwanted features
fripFull = fripFull[fripFull$feature_type %in% c("TSS","5UTR","3UTR","exons","introns","intergenic"),]

## Remove outlier sim sample
fripFull <- fripFull[fripFull$sampleID != 'K4_sim_sz11_m_noEtoh_rep3',]


## Plot
## Colors can be changed... look at http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
for(fileType in c("png","tiff","eps")){
  tempPlot = ggplot(fripFull, aes(x=feature_type, y=frip)) +
    geom_boxplot(aes(fill=sampleType_reorder)) +
    labs(title="FRiP Per Sample Type", y="FRiP", fill = "ChIP Antibody") +
    scale_fill_manual(values = c("red4", "red", "green4", "green", "dodgerblue4", "dodgerblue")) +
    theme_bw() + facet_grid(. ~ feature_type_reorder) +
    facet_wrap(~feature_type_reorder,scales="free",ncol=3) +
    theme(plot.title=element_text(hjust = 0.5),axis.title.x=element_blank(),
          axis.text.x=element_blank(),axis.ticks.x=element_blank())
  ggsave(plot=tempPlot,filename=paste(outFILE,".",fileType,sep=""),device=fileType,
         height=8, width=10)
}