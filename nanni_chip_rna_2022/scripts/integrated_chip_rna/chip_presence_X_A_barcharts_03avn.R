#!/usr/bin/env Rscript --vanilla

# Plot propotion of sex-biased histone mark on X and autosomes
#   (2 species counted separately but together in graph)
# Requires 6 arguments:
#   /path/to/species1_gene_annot.csv
#   /path/to/species2_gene_annot.csv
#   /path/to/ortho_gene_annot.csv
#   /path/to/output
#   species1
#   species2

# Ttest gene annotation files must contain at least "gene_sex_bias_ttest_foldchange" and "xsome" columns
# Gene annotation files must contain at least "gene_k4_sex", "gene_k27_sex" and "xsome" columns

library(readr)
library(ggplot2)
library(scales)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)

if(length(args)!=6) {
  stop("5 arguments required in this order: /path/to/species1_gene_annot.csv /path/to/species2_gene_annot.csv /path/to/ortho_gene_annot.csv /path/to/output species1 species2", call.=FALSE)
} else if (length(args)==6){
  if(endsWith(args[1],".csv")){
    inFile1 = args[1]
  } else {
    stop("First argument must be /path/to/species1_gene_annot.csv", call.=FALSE)
  }
  if(endsWith(args[2],".csv")){
    inFile2 = args[2]
  } else {
    stop("Second argument must be /path/to/species2_gene_annot.csv", call.=FALSE)
  }
  if(endsWith(args[3],".csv")){
    inFile3 = args[3]
  } else {
    stop("Third argument must be /path/to/ortho_gene_annot.csv", call.=FALSE)
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
one2one <- orthoDF[(orthoDF$flag_one2one_ortholog == 1),]
one2oneGene1 <- unlist(unique(one2one[,paste(name1,"_geneID", sep="")]), use.names = FALSE)
one2oneGene2 <- unlist(unique(one2one[,paste(name2,"_geneID", sep="")]), use.names = FALSE)

## Set species chromosome groups
geneDF1$chromGroup <- paste(fullName1, geneDF1$xsome)
geneDF2$chromGroup <- paste(fullName2, geneDF2$xsome)

## Set gene column names to be the same
names(geneDF2)[names(geneDF2) == geneCol2] <- geneCol1

## Combine both species with shared variables only (dropping GO and other literature lists)
## Drop rows with chromGroup that is not X or A
combGeneDF1 <- rbind(geneDF1[geneDF1$xsome %in% c("X","A"),
                            intersect(names(geneDF1), names(geneDF2))],
                    geneDF2[geneDF2$xsome %in% c("X","A"),
                            intersect(names(geneDF1), names(geneDF2))])

## Set categorical variables for k4 and k27 presence based on flags
combGeneDF2 <- mutate(combGeneDF1,
                      k4_group = case_when(flag_any_k4 == 0 ~ "No H3K4me3",
                                           flag_male_limited_k4 == 1 ~ "Male-limited H3K4me3",
                                           flag_female_limited_k4 == 1 ~ "Female-limited H3K4me3",
                                           flag_any_k4 == 1 ~ "H3K4me3 Present Both Sexes",
                                           TRUE ~ "oops"))
combGeneDF <- mutate(combGeneDF2,
                      k27_group = case_when(flag_any_k27 == 0 ~ "No H3K27me2me3",
                                           flag_male_limited_k27 == 1 ~ "Male-limited H3K27me2me3",
                                           flag_female_limited_k27 == 1 ~ "Female-limited H3K27me2me3",
                                           flag_any_k27 == 1 ~ "H3K27me2me3 Present Both Sexes",
                                           TRUE ~ "oops"))
## Subset on orthologs
combGeneOrtho <- combGeneDF[(unlist(combGeneDF[,paste(geneCol1)], use.names = FALSE) %in% one2oneGene1)
                           | (unlist(combGeneDF[,paste(geneCol1)], use.names = FALSE) %in% one2oneGene2),]

###### Plot proportion of H3K4me3 presence on X and A for all genes (including none)
## Get frequencies of chromGroup*k4_group
allopenDF <- as.data.frame(table(combGeneDF[,c("chromGroup","k4_group")]))
allopenOrtho <- as.data.frame(table(combGeneOrtho[,c("chromGroup","k4_group")]))

## Reorder gene_sex_bias_ttest_foldchange
allopenDF$k4_group_reorder = factor(allopenDF$k4_group,
                                     levels=c("No H3K4me3",
                                              "H3K4me3 Present Both Sexes",
                                              "Female-limited H3K4me3",
                                              "Male-limited H3K4me3"))
allopenOrtho$k4_group_reorder = factor(allopenOrtho$k4_group,
                                    levels=c("No H3K4me3",
                                             "H3K4me3 Present Both Sexes",
                                             "Female-limited H3K4me3",
                                             "Male-limited H3K4me3"))

for(fileType in c("png", "eps")){
  tempPlot = ggplot(allopenDF, aes(x = chromGroup, y = Freq, fill = k4_group_reorder)) +
    geom_col(colour = "white", position = "fill", width = 0.5) +
    labs(x="", y="", fill="H3K4me3") +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = c("black", "purple", "red", "blue")) +
    coord_flip() + theme_minimal()
  ggsave(plot=tempPlot,filename=paste(outFILE,"_all_H3K4me3.",fileType,sep=""),device=fileType,
         height=2, width=10)
}

for(fileType in c("png", "eps")){
  tempPlot = ggplot(allopenOrtho, aes(x = chromGroup, y = Freq, fill = k4_group_reorder)) +
    geom_col(colour = "white", position = "fill", width = 0.5) +
    labs(x="", y="", fill="H3K4me3") +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = c("black", "purple", "red", "blue")) +
    coord_flip() + theme_minimal()
  ggsave(plot=tempPlot,filename=paste(outFILE,"_all_H3K4me3_ortho.",fileType,sep=""),device=fileType,
         height=2, width=10)
}

###### Plot proportion of H3K4me3 on X and A only genes with a mark present
MFopenDF <- allopenDF[!(allopenDF$k4_group %in% c("No H3K4me3")),]

for(fileType in c("png","eps")){
  tempPlot = ggplot(MFopenDF, aes(x = chromGroup, y = Freq, fill = k4_group_reorder)) +
    geom_col(colour = "white", position = "fill", width = 0.5) +
    labs(x="", y="", fill="H3K4me3") +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = c("purple", "red", "blue")) +
    coord_flip() + theme_minimal()
  ggsave(plot=tempPlot,filename=paste(outFILE,"_MF_H3K4me3.",fileType,sep=""),device=fileType,
         height=2, width=10)
}


###### Plot proportion of H3K4me3 presence on X and A for all genes (including none)
## Get frequencies of chromGroup*gene_k27_sex
allcloseDF <- as.data.frame(table(combGeneDF[,c("chromGroup","k27_group")]))
allcloseOrtho <- as.data.frame(table(combGeneOrtho[,c("chromGroup","k27_group")]))

## Reorder gene_sex_bias_ttest_foldchange
allcloseDF$k27_group_reorder = factor(allcloseDF$k27_group,
                                      levels=c("No H3K27me2me3",
                                               "H3K27me2me3 Present Both Sexes",
                                               "Female-limited H3K27me2me3",
                                               "Male-limited H3K27me2me3"))
allcloseOrtho$k27_group_reorder = factor(allcloseOrtho$k27_group,
                                      levels=c("No H3K27me2me3",
                                               "H3K27me2me3 Present Both Sexes",
                                               "Female-limited H3K27me2me3",
                                               "Male-limited H3K27me2me3"))

for(fileType in c("png","eps")){
  tempPlot = ggplot(allcloseDF, aes(x = chromGroup, y = Freq, fill = k27_group_reorder)) +
    geom_col(colour = "white", position = "fill", width = 0.5) +
    labs(x="", y="", fill="H3K27me2me3") +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = c("black", "purple", "red", "blue")) +
    coord_flip() + theme_minimal()
  ggsave(plot=tempPlot,filename=paste(outFILE,"_all_H3K27me2me3.",fileType,sep=""),device=fileType,
         height=2, width=10)
}

for(fileType in c("png","eps")){
  tempPlot = ggplot(allcloseOrtho, aes(x = chromGroup, y = Freq, fill = k27_group_reorder)) +
    geom_col(colour = "white", position = "fill", width = 0.5) +
    labs(x="", y="", fill="H3K27me2me3") +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = c("black", "purple", "red", "blue")) +
    coord_flip() + theme_minimal()
  ggsave(plot=tempPlot,filename=paste(outFILE,"_all_H3K27me2me3_ortho.",fileType,sep=""),device=fileType,
         height=2, width=10)
}

###### Plot proportion of sex-biased H3K27me2me3 on X and A for male/female-biased genes only
MFcloseDF <- allcloseDF[!(allcloseDF$k27_group %in% c("No H3K27me2me3")),]

for(fileType in c("png","eps")){
  tempPlot = ggplot(MFcloseDF, aes(x = chromGroup, y = Freq, fill = k27_group_reorder)) +
    geom_col(colour = "white", position = "fill", width = 0.5) +
    labs(x="", y="", fill="H3K27me2me3") +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = c("purple", "red", "blue")) +
    coord_flip() + theme_minimal()
  ggsave(plot=tempPlot,filename=paste(outFILE,"_MF_H3K27me2me3.",fileType,sep=""),device=fileType,
         height=2, width=10)
}

## Output all counts
reorder_cols <- c("k4_group_reorder", "k27_group_reorder")
cat(paste("\nCounts for ",outFILE, "_all_H3K4me3:\n", sep=""))
print(allopenDF[order(allopenDF$chromGroup),!names(allopenDF) %in% reorder_cols], row.names = FALSE)
cat(paste("\nCounts for ",outFILE, "_all_H3K4me3_ortho:\n", sep=""))
print(allopenOrtho[order(allopenOrtho$chromGroup),!names(allopenOrtho) %in% reorder_cols], row.names = FALSE)
cat(paste("Counts for ",outFILE, "_MF_H3K4me3:\n", sep=""))
print(MFopenDF[order(MFopenDF$chromGroup),!names(MFopenDF) %in% reorder_cols], row.names = FALSE)
cat(paste("\nCounts for ",outFILE, "_all_H3K27me2me3:\n", sep=""))
print(allcloseDF[order(allcloseDF$chromGroup),!names(allcloseDF) %in% reorder_cols], row.names = FALSE)
cat(paste("\nCounts for ",outFILE, "_all_H3K27me2me3_ortho:\n", sep=""))
print(allcloseOrtho[order(allcloseOrtho$chromGroup),!names(allcloseOrtho) %in% reorder_cols], row.names = FALSE)
cat(paste("\nCounts for ",outFILE, "_MF_H3K27me2me3:\n", sep=""))
print(MFcloseDF[order(MFcloseDF$chromGroup),!names(MFcloseDF) %in% reorder_cols], row.names = FALSE)
