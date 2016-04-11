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

source('/home/fnew/ufgi_share/SHARE/McIntyre_Lab/scripts/R/wiggleplot_function_fnn2.R')

# Get command line arguments
## args[1] = exon annotation file
## args[2] = counts file
## args[3] = gene name
## args[4] = output name
args <- commandArgs(TRUE)
ename <- args[1]
cname <- args[2]
gene_symbol <- args[3]
oname <- args[4]


# Import exon annotation file (gene_symbol, chrom, exon_start, exon_end, trans_id)
exon_anno <- read.csv(ename)

# Import count file (chrom, pos, count1, count2, ...)
cov_count <- read.csv(cname)

## PNG
png(oname,width=1200,height=800,units="px")
#pdf(oname,width=1200,height=800)

par(mfrow=c(2,1));
plot.overlay.wiggle(cov_count,main.title=gene_symbol)
plot.gene.model(exon_anno,gene_symbol,model.ylim=40)

dev.off()
