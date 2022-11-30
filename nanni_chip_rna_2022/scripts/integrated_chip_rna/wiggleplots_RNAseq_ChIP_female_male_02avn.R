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

source('./wiggleplot_function_04avn.R')

# Get command line arguments
#args[1] = exon annotation file
#args[2] = counts file
#args[3] = gene name
#args[4] = output name
args <- commandArgs(TRUE)
ename <- args[1]
cname <- args[2]
gene_symbol <- args[3]
oname <- args[4]


# Import exon annotations

# Import exon annotation file (gene_symbol, chrom, exon_start, exon_end, trans_id)
exon_anno <- read.csv(ename)

# Import count file (chrom, pos, count1, count2, ...)
cov_count <- read.csv(cname)

## PNG
# Split dataframe into coords and counts
coords <- cov_count[,c(1,2)]
counts <- cov_count[,c(2,3,4)]

# Calc maximum counts and reorder
colorder <- apply(counts, 2, max)
mydatord <- counts[,order(-colorder)]

# Add back into dataframe
ordered_counts <- merge(coords, mydatord, by="pos")
ordered <- cov_count[,c(2,1,3,4)]

png(oname,width=1200,height=1200,units="px")
#par(mfrow=c(2,1));
mat <- matrix(c(1,2),nrow=2,ncol=1)
layout(mat, c(1), c(1,2))
plot.overlay.wiggle(cov_count,main.title=gene_symbol)
plot.gene.model(exon_anno,gene_symbol,model.ylim=300)

dev.off()
