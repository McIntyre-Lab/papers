library(lavaan)
library(simsem)

# Set ENVs and Folders
proj <- '/scratch/lfs/mcintyre/cegs_sem_sd_paper'

# Read in FILES
## Original data for running baseline SEM to simulate True data
origDat <- read.csv(paste0(proj,'/dspr_adding_links_simulation/dspr_sbs_gene_level_sym.csv'), header=TRUE)

numSamples <- 75

# Create True model (no covariance model)
model <- '
    Sxl ~ snf + Spf45 + vir + fl_2_d
    tra ~ vir + fl_2_d + Sxl
    fru ~ tra + tra2
    Yp2 ~ her + ix + tra + tra2

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

# Create simulated datasets
# Create True simulation by simulating from reald SEM data (creates dataset modelDat)
modelDat <- generate(fit, numSamples)

write.csv(modelDat,row.names=FALSE, quote=FALSE)
