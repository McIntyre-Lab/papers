#!/usr/bin/env Rscript --vanilla

# Plot FRiP values of MACS2 consensus peaks
# Requires 4 arguments: /path/to/summary_frip.csv  /path/to/output.png species1 species2

# FRiP summary file must contain at least "sampleID" and "frip" columns

library(readr)
library(ggplot2)

args = commandArgs(trailingOnly = TRUE)

if(length(args)!=4) {
	stop("4 arguments required in this order: /path/to/summary_frip.csv /path/to/output_prefix species1 species2", call.=FALSE)
} else if (length(args)==4){
	if(endsWith(args[1],".csv")){
		inFILE = args[1]
	} else {
		stop("First argument must be /path/to/summary_frip.csv", call.=FALSE)
	}
	outFILE = args[2]
	name1 = args[3]
	name2 = args[4]
}

summary_full <- read_csv(inFILE)

frip <- summary_full[c('sampleID','frip')]

frip$species <- ifelse(grepl(name1,frip$sampleID),name1,name2)

frip$sampleType <- ifelse(grepl("K27_cvrg_cnts_I_",frip$sampleID),"IgG Control (H3K27me2me3 Peaks)",ifelse(grepl("K4_cvrg_cnts_I_",frip$sampleID),"IgG Control (H3K4me3 peaks)",ifelse(grepl("^K4_",frip$sampleID),"H3K4me3","H3K27me2me3")))

frip.noEtoh <- subset(frip,grepl("noEtoh",sampleID))

for(chip in c("H3K4me3","H3K27me2me3")){
  for(fileType in c("png","tiff","eps")){
    fileName = paste(outFILE,"_",chip,".",fileType,sep="")

    tempPlot = ggplot(subset(frip.noEtoh,grepl(chip,sampleType)),aes(x=sampleType, y=frip)) +
                  geom_boxplot(aes(fill=species)) +
                  labs(title=paste(chip," FRiP",sep=""), y="FRiP", colour = "Species") +
                  theme_bw() + facet_grid(. ~ species) + facet_wrap(~species,scales="free",ncol=2) +
                  theme(plot.title=element_text(hjust = 0.5),axis.title.x=element_blank(),
                      axis.text.x=element_text(angle=45, hjust=1),
                      strip.background=element_blank(), strip.text = element_blank()) +
                  scale_fill_discrete(name="Species",breaks=c("mel","sim"),
                      labels=c(expression(italic("D. melanogaster")),expression(italic("D. simulans"))))

    ggsave(plot=tempPlot,file=fileName,device=fileType)
  }
}
