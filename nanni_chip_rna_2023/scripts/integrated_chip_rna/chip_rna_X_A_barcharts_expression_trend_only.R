#!/usr/bin/env Rscript --vanilla

# Plot propotion of sex-biased expression trend on X and autosomes (2 species counted separately but together in graph)
# Requires 5 arguments: /path/to/species1_ttest_gene_annot.csv /path/to/species2_ttest_gene_annot.csv /path/to/output species1 species2

# Ttest gene annotation files must contain at least "gene_trend_ttest" and "xsome" columns

library(readr)
library(ggplot2)
library(scales)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)

if(length(args)!=6) {
  stop("7 arguments required in this order: /path/to/species1_gene_annot.csv /path/to/species2_gene_annot.csv /path/to/ortholog_gene_annot.csv /path/to/output species1 species2", call.=FALSE)
} else if (length(args)==6){
  if(endsWith(args[1],".csv")){
    inFile1 = args[1]
  } else {
    stop("First argument must be /path/to/species1_ttest_gene_annot.csv", call.=FALSE)
  }
  if(endsWith(args[2],".csv")){
    inFile2 = args[2]
  } else {
    stop("Second argument must be /path/to/species2_ttest_gene_annot.csv", call.=FALSE)
  }
  if(endsWith(args[3],".csv")){
    inFile3 = args[3]
  } else {
    stop("Third argument must be /path/to/ortholog_ttest_gene_annot.csv", call.=FALSE)
  }
  outFILE = args[4]
  name1 = args[5]
  name2 = args[6]
}

if(name1 == "mel"){
  fullName1 = "D. melanogaster"
  geneCol1 = "FBgn"
} else if(name1 == "sim"){
  fullName1 = "D. simulans"
  geneCol1 = "fbgn"
}
if(name2 == "mel"){
  fullName2 = "D. melanogaster"
  geneCol2 = "FBgn"
} else if(name2 == "sim"){
  fullName2 = "D. simulans"
  geneCol2 = "fbgn"
}

geneDF1 <- read_csv(inFile1)
geneDF2 <- read_csv(inFile2)
orthoDF <- read_csv(inFile3)

# Get set of one-to-one orthologous genes from ortholog file (for ortho subset plots)
# Drop undetected and sex-limited genes
one2one <- orthoDF[(orthoDF$flag_one2one_ortholog == 1),]
one2oneGene1 <- unlist(unique(one2one[,paste(name1,"_geneID", sep="")]), use.names = FALSE)
one2oneGene2 <- unlist(unique(one2one[,paste(name2,"_geneID", sep="")]), use.names = FALSE)

## Drop undetected and sex-limited genes
expressDF1 <- geneDF1[geneDF1$flag_expressed==1 & geneDF1$flag_sex_limited==0,]
expressDF2 <- geneDF2[geneDF2$flag_expressed==1 & geneDF2$flag_sex_limited==0,]

## Set species chromosome groups
expressDF1$chromGroup <- paste(fullName1, expressDF1$xsome)
expressDF2$chromGroup <- paste(fullName2, expressDF2$xsome)

## Set gene column names to be the same
names(expressDF2)[names(expressDF2) == geneCol2] <- geneCol1

## Set sex-bias group
expressDF1 <- expressDF1 %>%
    mutate(., sexGroup = with(., case_when(
        expressDF1$flag_M==1 ~ "Male-bias (Overcompensated)",
        expressDF1$flag_F==1 ~ "Female-bias (Undercompensated)",
        expressDF1$flag_MF==1 ~ "Male and Female-bias (Overcompensated and Undercompensated)",
        expressDF1$flag_U==1 ~ "Unbiased (Compensated)",
        TRUE ~ "OOPS"
    )))
expressDF2 <- expressDF2 %>%
    mutate(., sexGroup = with(., case_when(
        expressDF2$flag_M==1 ~ "Male-bias (Overcompensated)",
        expressDF2$flag_F==1 ~ "Female-bias (Undercompensated)",
        expressDF2$flag_MF==1 ~ "Male and Female-bias (Overcompensated and Undercompensated)",
        expressDF2$flag_U==1 ~ "Unbiased (Compensated)",
        TRUE ~ "OOPS"
    )))


## Combine both species sexGroup and chromGroup
## Drop rows with chromGroup that is not X or A
combTTESTDF <- rbind(expressDF1[expressDF1$xsome %in% c("X","A"),
                              c(geneCol1, "sexGroup","chromGroup")],
                     expressDF2[expressDF2$xsome %in% c("X","A"),
                              c(geneCol1, "sexGroup","chromGroup")])
combTTESTDFortho <- combTTESTDF[(unlist(combTTESTDF[,paste(geneCol1)], use.names = FALSE) %in% one2oneGene1)
                                | (unlist(combTTESTDF[,paste(geneCol1)], use.names = FALSE) %in% one2oneGene2),]

###### Plot proportion of sex-biased expression on X and A for all expressed genes (including unbiased)

# Frequencies
allExpDF <- as.data.frame(table(combTTESTDF[,c("chromGroup","sexGroup")]))
allExpDFortho <- as.data.frame(table(combTTESTDFortho[,c("chromGroup","sexGroup")]))

# Percentages
allExpPercDF <- as.data.frame(prop.table(table(combTTESTDF[,c("chromGroup","sexGroup")]),1))
allExpPercDF$Percent <- allExpPercDF$Freq*100
allExpPercDFortho <- as.data.frame(prop.table(table(combTTESTDFortho[,c("chromGroup","sexGroup")]),1))
allExpPercDFortho$Percent <- allExpPercDFortho$Freq*100
# Set frequency in percent table
allExpPercDF$Freq <- allExpDF$Freq
allExpPercDFortho$Freq <- allExpDFortho$Freq

## Reorder sexGroup
allExpDF$sexGroup_reorder = factor(allExpDF$sexGroup,levels=c(
  'Unbiased (Compensated)',
  'Male and Female-bias (Overcompensated and Undercompensated)',
  'Female-bias (Undercompensated)',
  'Male-bias (Overcompensated)'))
allExpPercDF$sexGroup_reorder = factor(allExpPercDF$sexGroup,levels=c(
  'Unbiased (Compensated)',
  'Male and Female-bias (Overcompensated and Undercompensated)',
  'Female-bias (Undercompensated)',
  'Male-bias (Overcompensated)'))
allExpDFortho$sexGroup_reorder = factor(allExpDFortho$sexGroup,levels=c(
  'Unbiased (Compensated)',
  'Male and Female-bias (Overcompensated and Undercompensated)',
  'Female-bias (Undercompensated)',
  'Male-bias (Overcompensated)'))
allExpPercDFortho$sexGroup_reorder = factor(allExpPercDFortho$sexGroup,levels=c(
  'Unbiased (Compensated)',
  'Male and Female-bias (Overcompensated and Undercompensated)',
  'Female-bias (Undercompensated)',
  'Male-bias (Overcompensated)'))

## Plot proportions using frequencies
for(fileType in c("png","tiff","eps")){
  tempPlot = ggplot(allExpDF, aes(x = chromGroup, y = Freq, fill = sexGroup_reorder)) +
    geom_col(colour = "white", position = "fill", width = 0.5) +
    labs(x="", y="", fill="Expression Sex Bias") +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = c("gray","mediumpurple","salmon","lightblue")) +
    coord_flip() + theme_minimal()
  ggsave(plot=tempPlot,filename=paste(outFILE,"_all_Expression.",fileType,sep=""),device=fileType,
         height=2, width=15)
}

## Plot proportions using percentages (dropping unbiased)
## Set axis to be the max of each sex-bias group + the associated pronounced (rounded to nearest 5)
maxY <- max(c(
  aggregate(Percent ~ chromGroup, allExpPercDF[
    allExpPercDF$sexGroup %in% c(
      "Male and Female-bias (Overcompensated and Undercompensated)"
    ),], sum)$Percent,
  aggregate(Percent ~ chromGroup, allExpPercDF[
    allExpPercDF$sexGroup %in% c(
      "Male-bias (Overcompensated)"
    ),], sum)$Percent,
  aggregate(Percent ~ chromGroup, allExpPercDF[
    allExpPercDF$sexGroup %in% c(
      "Female-bias (Undercompensated)"
    ),], sum)$Percent))
maxYround <- ceiling(maxY / 5)*5
maxYortho <- max(c(
  aggregate(Percent ~ chromGroup, allExpPercDFortho[
    allExpPercDFortho$sexGroup %in% c(
      "Male and Female-bias (Overcompensated and Undercompensated)"
    ),], sum)$Percent,
  aggregate(Percent ~ chromGroup, allExpPercDFortho[
    allExpPercDFortho$sexGroup %in% c(
      "Male-bias (Overcompensated)"
    ),], sum)$Percent,
  aggregate(Percent ~ chromGroup, allExpPercDFortho[
    allExpPercDFortho$sexGroup %in% c(
      "Female-bias (Undercompensated)"
    ),], sum)$Percent))
maxYorthoround <- ceiling(maxYortho / 5)*5

# Plot split bar charts for all genes and orthologs (use max of the two Y limits so they are on the same scale)
for(fileType in c("png","tiff","eps")){
  tempPlot = ggplot(allExpPercDF[allExpPercDF$sexGroup_reorder!="Unbiased (Compensated)",], aes(x = chromGroup, y = Percent, fill = sexGroup_reorder)) +
    geom_col(colour = "white", position = "dodge", width=0.8) +
    geom_text(aes(label = Freq),  position = position_dodge(0.8), vjust = 0, size=2.5) +
    labs(x="", y="% Expressed Genes", fill="Expression Sex Bias") +
    scale_fill_manual(values = c("mediumpurple","salmon","lightblue")) +
    theme_minimal() + ylim(0, max(maxYround, maxYorthoround))
  ggsave(plot=tempPlot,filename=paste(outFILE,"_sexBias_Expression_split_bar.",fileType,sep=""),device=fileType,
         height=4, width=10)
  # Plot one-to-one orthologs only
  tempPlot = ggplot(allExpPercDFortho[allExpPercDFortho$sexGroup_reorder!="Unbiased (Compensated)",], aes(x = chromGroup, y = Percent, fill = sexGroup_reorder)) +
    geom_col(colour = "white", position = "dodge", width=0.8) +
    geom_text(aes(label = Freq),  position = position_dodge(0.8), vjust = 0, size=2.5) +
    labs(x="", y="% Expressed Genes", fill="Expression Sex Bias") +
    scale_fill_manual(values = c("mediumpurple","salmon","lightblue")) +
    theme_minimal() + ylim(0, max(maxYround, maxYorthoround))
  ggsave(plot=tempPlot,filename=paste(outFILE,"_sexBias_Expression_split_bar_ortho.",fileType,sep=""),device=fileType,
         height=4, width=10)
}


###### Plot proportion of sex-biased expression on X and A for male/female-biased expressed genes only
MFExpDF <- allExpDF[allExpDF$sexGroup %in% c(
  'Male and Female-bias (Overcompensated and Undercompensated)',
  'Male-bias (Overcompensated)',
  'Female-bias (Undercompensated)'),]

for(fileType in c("png","tiff","eps")){
  tempPlot = ggplot(MFExpDF, aes(x = chromGroup, y = Freq, fill = sexGroup_reorder)) +
    geom_col(colour = "white", position = "fill", width = 0.5) +
    labs(x="", y="", fill="Expression Sex Bias") +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = c("mediumpurple","salmon","lightblue")) +
    coord_flip() + theme_minimal()
  ggsave(plot=tempPlot,filename=paste(outFILE,"_MF_Expression.",fileType,sep=""),device=fileType,
         height=2, width=15)
}


## Output all counts
reorder_cols <- c("sexGroup")
cat(paste("Counts for ",outFILE, "_all_Expression:\n", sep=""))
print(allExpDF[order(allExpDF$chromGroup),!names(allExpDF) %in% reorder_cols], row.names = FALSE)
cat(paste("Counts for ",outFILE, "_all_Expression_ortho:\n", sep=""))
print(allExpDF[order(allExpDFortho$chromGroup),!names(allExpDFortho) %in% reorder_cols], row.names = FALSE)
cat(paste("\nCounts for ",outFILE, "_MF_Expression:\n", sep=""))
print(MFExpDF[order(MFExpDF$chromGroup),!names(MFExpDF) %in% reorder_cols], row.names = FALSE)

