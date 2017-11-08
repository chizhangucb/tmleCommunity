#***************************************************************************************
# Example: storing, managing, subsetting and manipulating a data with continuous A
data(indSample.iid.cA.cY_list)
indSample.iid.cA.bY <- indSample.iid.cA.bY_list$indSample.iid.cA.bY
psi0.Y <- indSample.iid.cA.bY_list$psi0.Y  # 0.333676
# Assume that W2 has no effect on neither A nor Y, so no need to put into nodes
nodes <- list(Ynode = "Y", Anodes = "A", WEnodes = c("W1", "W3", "W4"))  
tmleCom_Options(nbins = 10, maxNperBin = nrow(indSample.iid.cA.bY))
#***************************************************************************************

#***************************************************************************************
# 1.1 Specifying the stochastic intervention of interest
#***************************************************************************************
# Interested in the effect of a shift of delta(W1, W3, W4) of the current treatment
define_f.gstar <- function(data, ...) {
  shift.val <- 0.3 * data[,"W1"] + 0.6 * data[,"W3"] - 0.14 * data[,"W4"]
  shifted.new.A <- data[, "A"] - shift.val
  return(shifted.new.A)
}

#***************************************************************************************
# 1.2 Creating an R6 object of DatKeepClass (to store the input data)
#***************************************************************************************
# Don't normalize continous covariates by setting norm.c.sVars = FALSE
OData <- DatKeepClass$new(Odata = indSample.iid.cA.bY, nodes = nodes, norm.c.sVars = FALSE)  
# names of all variables that are in input data and specified in nodes
OData$names.sVar  # "Y"  "A"  "W1" "W3" "W4" 
# names of all continuous variables that are in input data and specified in nodes
OData$names.c.sVar  # "A"  "W3" "W4" 
# a sub dataframe of the input data, including all variables in nodes
head(OData$dat.sVar) 
# the number of observations of the input data
OData$nobs  # 10000
OData$get.sVar.type("A")  # "contin"
OData$get.sVar.type()  # Provide a list of types of all variables

#***************************************************************************************
# 1.3 Manipulating the input data by adding observed outcome and observation weights
#***************************************************************************************
OData$addYnode(YnodeVals = indSample.iid.cA.bY[, nodes$Ynode], det.Y = FALSE)  
head(OData$YnodeVals)  # Adding public YnodeVals & setting det.Y values to NA
head(OData$noNA.Ynodevals)  # Adding actual observed Y as protected (without NAs)
OData$addObsWeights(obs.wts = rep(c(1,2), 5000))  # Add a vectopr of observation (sampling) weights
OData$addObsWeights(obs.wts = 1)  # If not specified, assumed to be all 1 (i.e., equally weighted)

#***************************************************************************************
# 1.4 
#***************************************************************************************
OData.gstar <- DatKeepClass$new(Odata = indSample.iid.cA.bY, nodes = nodes, norm.c.sVars = FALSE)
# Create (p=) 1 new Odata (nobs obs at a time) and replace A under g0 in Odata with A^* under g.star
# Generate new exposures under user-specific arbitrary intervention f.g_fun.
OData.gstar$make.dat.sVar(p = 1, f.g_fun = define_f.gstar) 
dim(OData.gstar$dat.sVar)  # 10000     4
OData.gstar$make.dat.sVar(p = 3, f.g_fun = define_f.gstar) 
dim(OData.gstar$dat.sVar)  # 30000     4
