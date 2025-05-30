library(foreign)
library(lavaan)
library(simsem)

# Set ENVs and Folders
proj <- '/scratch/lfs/mcintyre/cegs_sem_sd_paper'

# Read in FILES
## Original data for running baseline SEM to simulate True data
origDat <- read.csv(paste0(proj,'/cegs_adding_genes_simulation/cegsv_by_gene_sbs.csv'), header=TRUE)

## Null mean and variance data
nullDat <- read.csv(paste0(proj,'/cegs_adding_genes_simulation/null_mean_and_variances.csv'), colClasses=c('NULL', 'numeric', 'numeric'), header=TRUE)

# Read in command line arguments
args <- commandArgs(TRUE)
if (length(args) > 0){
    numGenes <- as.numeric(args[1])
} else {
    numGenes <- 7419
}

# Create True model (no covariance model)
model <- '
    Sxl ~ snf + Spf45 + vir + fl_2_d
    tra ~ vir + fl_2_d + Sxl
    fru ~ tra + tra2
    dsx ~ tra + tra2
    Yp2 ~ her + ix + dsx

    # Exogenous Cov set to 0
    fl_2_d ~~ 0*her + 0*ix + 0*snf + 0*Spf45 + 0*tra2 + 0*vir
    her    ~~ 0*ix + 0*snf + 0*Spf45 + 0*tra2 + 0*vir 
    ix     ~~ 0*snf + 0*Spf45 + 0*tra2 + 0*vir
    snf    ~~ 0*Spf45 + 0*tra2 + 0*vir 
    Spf45  ~~ 0*tra2 + 0*vir
    tra2   ~~ 0*vir
    '

fit <- lavaan(model, origDat, fixed.x=FALSE, int.ov.free=TRUE, int.lv.free=FALSE,
              auto.fix.first=TRUE, auto.fix.single=TRUE, auto.var=TRUE,
              auto.cov.lv.x=TRUE, auto.th=TRUE, auto.delta=TRUE,
              auto.cov.y=FALSE)

# Function to simulate NULL data based on real 
simNull <- function(mydat, numGenes){
    # Randomly sample mean-variance from null data (with replacement)
    simDist <- mydat[sample(nrow(mydat), numGenes, replace=TRUE),]
    row.names(simDist) <- paste0('gene',seq(1,nrow(simDist)))

    # Simulate 75 samples from normal distribution using mean-variance from null data
    simDat <- apply(simDist, 1, function(x) rnorm(75, x[1], x[2]))
    return(simDat)
}

# Create simulated datasets
# Create True simulation by simulating from reald SEM data (creates dataset modelDat)
modelDat <- generate(fit, 75)

# Create NULL simulation (creates dataset simDat)
simDat <- simNull(nullDat, numGenes)

# Merge data sets together and output
newDat <- cbind(simDat, modelDat)
write.foreign(newDat, args[2], args[3], package="SAS")
