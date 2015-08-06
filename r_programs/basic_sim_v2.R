library(lavaan)
library(simsem)

# read in the null distribution of means and variances
mclab <- Sys.getenv('MCLAB')
mydat <- read.csv(paste0(mclab,'/cegs_sem_sd_paper/exported_data/cegsv_by_gene_sbs.csv'), header=TRUE)

# no cov model
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

fit <- lavaan(model, mydat, fixed.x=FALSE, int.ov.free=TRUE, int.lv.free=FALSE,
              auto.fix.first=TRUE, auto.fix.single=TRUE, auto.var=TRUE,
              auto.cov.lv.x=TRUE, auto.th=TRUE, auto.delta=TRUE,
              auto.cov.y=FALSE,meanstructure=TRUE)

#output <- sim(1000, n=nrow(mydat), model=fit, generate=fit)
modelDat <- generate(fit, 75)
