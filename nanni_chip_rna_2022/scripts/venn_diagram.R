#!/usr/bin/env Rscript --vanilla

# Plot venn diagram of 2, 3, or 4 groups

# require 2 comma separated lists and output file prefix:
# venn_diagram.R labels areas /path/to/output

# areas is the list of areas in the order of...
#     for 2: 10, 01, 11
#     for 3: 100, 010, 110, 001, 101, 011, 111
#     for 4: 1000, 0100, 1100, 0010, 1010, 0110, 1110, 0001,
#            1001, 0101, 1101, 0011, 1011, 0111, 1111


#install.packages("remotes")
#remotes::install_github("js229/Vennerable")

library(Vennerable)

getVenn <- function(n,v) {
  # Need to fix labels
  # Use something like https://rdrr.io/rforge/Vennerable/man/VennSetSetLabels.html
  V1 <- Venn(SetNames = n, Weight = c(0,v))
  gp <- VennThemes(compute.Venn(V1), colourAlgorithm = "sequential")
  if(length(n) <= 3){
    out <- plot(V1,show=list(Universe=FALSE), gp=gp)
  } else {
    out <- plot(V1,type='ellipses', show=list(Universe=FALSE), gp=gp)
  }
  return(out)
}

args = commandArgs(trailingOnly = TRUE)

names <- unlist(strsplit(args[1],","))

values <- as.integer(unlist(strsplit(args[2],",")))

outFile <- args[3]

png(paste(outFile,".png",sep=""))
getVenn(names,values)
dev.off()