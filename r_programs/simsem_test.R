# This is a test of the simSEM package.

library(simsem)

# Create factor loading matrix
loading <- matrix(0, 6, 2) #create a matrix of all 0s
loading[1:3, 1] <- c("a1", "a2", "a3") #specify free parameters with texts
loading[4:6, 2] <- c("a4", "a5", "a6")

loadingValues <- matrix(0, 6, 2)  
loadingValues[1:3, 1] <- 0.7  
loadingValues[4:6, 2] <- 0.7
LY <- bind(loading, loadingValues)

# Create psi, error corr matrix
error.cor <- matrix(0, 6, 6)
diag(error.cor) <- 1
RTE <- binds(error.cor)

# create latent corr matrix
latent.cor <- matrix(NA, 2, 2) #specify a 2x2 matrix of NAs
diag(latent.cor) <- 1 #set the diagonal of the matrix to 1
RPS <- binds(latent.cor, 0.5) #Defaults to making all NA values in the matrix .5

# Put them all together to make a model object
CFA.Model <- model(LY = LY, RPS = RPS, RTE = RTE, modelType="CFA")

# simulate 200 data points and save
dat <- generate(CFA.Model, 200)
write.csv(dat, file='/home/jfear/mclab/cegs_sem_sd_paper/simulation/test_cfa.csv')

# Analyze dataset using the lavaan package and get summary.
out <- analyze(CFA.Model, dat)
parameterEstimates(out)
fitMeasures(out)
resid(out)
fitted(out)




# OpenMX specification

Avalues <- matrix(0, 8, 8)
Avalues[1:3, 7] <- c(0.7, 0.7, 0.7) # Factor loadings from f1 to x1-x3
Avalues[4:6, 8] <- c(0.7, 0.7, 0.7) # Factor loadings from f1 to x1-x3

Afree <- matrix(FALSE, 8, 8)
Afree[1:3, 7] <- c(TRUE, TRUE, TRUE) # The marker-variable approach is used
Afree[4:6, 8] <- c(TRUE, TRUE, TRUE)

A <- mxMatrix(type="Full", nrow=8, ncol=8, values=Avalues, free=Afree, byrow=TRUE, name="A")


# Error variances of nine indicators and then the variances of three factors
Svalues <- diag(c(0.51, 0.51, 0.51, 0.51, 0.51, 0.51, 0.51, 0.51)) 
# Factor covariances
Svalues[7, 8] <- Svalues[7, 8] <- 0.5 # Covariance between f1 and f2
Svalues[8, 7] <- Svalues[8, 7] <- 0.5 # Covariance between f1 and f2


Sfree <- matrix(FALSE, 8, 8)
diag(Sfree) <- TRUE
Sfree[7, 8] <- Sfree[7, 8] <- TRUE
Sfree[8, 7] <- Sfree[8, 7] <- TRUE
S <- mxMatrix(type="Symm", nrow=8, ncol=8, values=Svalues, free=Sfree, byrow=TRUE, name="S")

Fvalues <- cbind(diag(6), matrix(0, 6, 3))
F <- mxMatrix(type="Full", nrow=6, ncol=8, free=FALSE, values=Fvalues, byrow=TRUE, name="F")

Mvalues <- rep(0, 8)

Mfree <- c(rep(TRUE, 6), rep(FALSE, 2))

M <- mxMatrix(type="Full", nrow=1, ncol=8, values=Mvalues, free=Mfree, name="M")


popModel <- mxModel("Three Factor Model",
    type="RAM", 
    A, S, F, M,
    mxRAMObjective("A","S","F","M", dimnames=c(paste0("y", 1:6), "f1", "f2"))
)



out2 <- analyze(popModel, dat)
summary(out)
