###################################################################################
#
# ESTIMATION OF EXPRESSION RESIDUALS USING PEER
# I am using an expression dataset consisting of all fusions (exonic regions) that
# are flagged as "on" in all three cell types. This consists of 163713 fusions
#
# PEER reference: Stegle, O. et al, 2012, Nat Protoc, 7(3), 500-507
# PEER can be downloaded here: https://github.com/PMBio/peer
#
####################################################################################

# Import libraries
import peer
import scipy as SP
import numpy

# Load expression data. This has been formatted as samples as columns, fusions as rows
expr = SP.loadtxt('expr_data_for_peer.csv', delimiter=',')
expr.shape

# Need to transpose the data, so that samples are rows, fusions are columns
expr_t = numpy.transpose(expr)
expr_t.shape

# Set the model object
model = peer.PEER()

# Add means to model
model.setPhenoMean(expr_t)
model.getPhenoMean().shape

# Set the number of factors to output (n=10)
model.setNk(10)
model.getNk()

# Update the model
model.update()

# Get output: factors, weights, precision, residuals
factors = model.getX()
factors.shape
weights = model.getW()
weights.shape
precision = model.getAlpha()
precision.shape
residuals = model.getResiduals()
residuals.shape

# Export data
numpy.savetxt("peer_factors.csv", factors, delimiter=",", fmt="%10.18f")
numpy.savetxt("peer_weights.csv", weights, delimiter=",", fmt="%10.18f")
numpy.savetxt("peer_precision.csv", precision, delimiter=",", fmt="%10.18f")
numpy.savetxt("peer_residuals.csv", residuals, delimiter=",", fmt="%10.18f")

