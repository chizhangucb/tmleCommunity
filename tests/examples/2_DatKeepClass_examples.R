#***************************************************************************************
# Example 1: storing, managing, subsetting and manipulating a data with continuous A
data(indSample.iid.cA.cY_list)
indSample.iid.cA.bY <- indSample.iid.cA.bY_list$indSample.iid.cA.bY
psi0.Y <- indSample.iid.cA.bY_list$psi0.Y  # 0.333676
# Assume that W2 has no effect on neither A nor Y, so no need to put into nodes
nodes <- list(Ynode = "Y", Anodes = "A", WEnodes = c("W1", "W3", "W4"))  
tmleCom_Options(nbins = 10, maxNperBin = nrow(indSample.iid.cA.bY))
#***************************************************************************************

#***************************************************************************************
# 1.1 Specifying the stochastic intervention of interest gstar
#***************************************************************************************
# Interested in the effect of a shift of delta(W1, W3, W4) of the current treatment
define_f.gstar <- function(data, ...) {
  shift.mu <- 0.3 * data[,"W1"] + 0.6 * data[,"W3"] - 0.14 * data[,"W4"]
  shift.val <- rnorm(n = NROW(data), mean = shift.mu, sd = 0.5)
  shifted.new.A <- data[, "A"] - shift.val
  return(shifted.new.A)
}

#***************************************************************************************
# 1.2 Creating an R6 object of DatKeepClass (to store the input data)
#***************************************************************************************
# Don't normalize continous covariates by setting norm.c.sVars = FALSE
OData_R6 <- DatKeepClass$new(Odata = indSample.iid.cA.bY, 
                             nodes = nodes, norm.c.sVars = FALSE)  
# names of all variables that are in input data and specified in nodes
OData_R6$names.sVar  # "Y"  "A"  "W1" "W3" "W4" 
# names of all continuous variables that are in input data and specified in nodes
OData_R6$names.c.sVar  # "A"  "W3" "W4" 
# a sub dataframe of the input data, including all variables in nodes
head(OData_R6$dat.sVar) 
# the number of observations of the input data
OData_R6$nobs  # 10000
OData_R6$get.sVar.type("A")  # "contin"
OData_R6$get.sVar.type()  # Provide a list of types of all variables

#***************************************************************************************
# 1.3 Manipulating the input data by adding observed outcomes and observation weights
#***************************************************************************************
# Add YnodeVals (a vector of outcomes) to both public and private field 
OData_R6$addYnode(YnodeVals = indSample.iid.cA.bY[, nodes$Ynode], det.Y = FALSE)  
head(OData_R6$YnodeVals)  # set YnodeVals[det.Y=TRUE] to NA in public field
head(OData_R6$noNA.Ynodevals)  # YnodeVals in private field is protected (without NAs)  
OData_R6$addObsWeights(obs.wts = rep(c(1,2), 5000))  # Add a vector of sampling weights
ODaOData_R6ta$addObsWeights(obs.wts = 1)  # Assumed to be all 1 (i.e., equally weighted)

#***************************************************************************************
# 1.4 Creating an new R6 object of DatKeepClass under stochastic intervention g.star
# Generate new exposures under user-specific intervention f.g_fun
#***************************************************************************************
OData.gstar_R6 <- DatKeepClass$new(Odata = indSample.iid.cA.bY, nodes = nodes)
# Create 1 new Odata and replace A under g0 in Odata with A* under g.star
set.seed(12345)
OData.gstar_R6$make.dat.sVar(p = 1, f.g_fun = define_f.gstar) 
dim(OData.gstar_R6$dat.sVar)  # 10000     4
# Create 3 new Odatas and repalce A with A*
OData.gstar_R6$make.dat.sVar(p = 3, f.g_fun = define_f.gstar) 
dim(OData.gstar_R6$dat.sVar)  # 30000     4
# Since A* is stochastically generated, each p may produce different values of A*
head(OData.gstar_R6$dat.sVar[1:10000, ])
head(OData.gstar_R6$dat.sVar[10001:20000, ])
