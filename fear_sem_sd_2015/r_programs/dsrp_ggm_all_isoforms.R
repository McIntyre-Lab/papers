# Based off of code found at:
# http://strimmerlab.org/software/genenet/download/ecoli-net.R

library(GeneNet)
mclab <- Sys.getenv("MCLAB")

# Import Side-by-side data where each column is an isoform
mydata <- read.csv('/tmp/data_for_ggm.csv',header=TRUE)

# Set row names to the RIL id (patRIL_matRIL)
row.names(mydata) <- mydata$sample

# Drop RIL id
mydata <- subset(mydata, select = -c(1))

# Calculate Partial Correlations using shrinkage estimator
pc <- ggm.estimate.pcor(mydata)

# Estimate p and q-values, along with posterior probabilites
mydata.edges <- network.test.edges(pc, direct=TRUE, fdr=TRUE)

# Only keep edges with a local fdr 0.2 (cutoff.ggm = 1-0.2)
mydata.net <- extract.network(mydata.edges, cutoff.ggm=0.8, cutoff.dir=0.8)

# Plot the Network
node.labels <- colnames(mydata)
network.make.dot(paste(mclab,"/cegs_sem_sd_paper/analysis_output/ggm/dsrp_ggm_isoforms_FDR2.dot",sep=""),mydata.net, node.labels, show.edge.labels=TRUE, main="Isoform GGM")
#system("dot -T png -o $MCLAB/cegs_sem_sd_paper/analysis_output/ggm/dsrp_ggm_isoforms_FDR2.png $MCLAB/cegs_sem_sd_paper/analysis_output/ggm/dsrp_ggm_isoforms_FDR2.dot")
