### Load data
# A is normal with mu for each observation being a function of (W1, W2, W3, W4), sd = 1;
data(sampleDat_iidcontAContY)
dat_iidcontAContY <- sampleDat_iidcontAContY$dat_iidcontAContY
psi0.Y <- sampleDat_iidcontAContY$psi0.Y  # 3.309084
psi0.Ygstar <- sampleDat_iidcontAContY$psi0.Ygstar  # 3.50856
nodes <- list(Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), Enodes = NULL, Crossnodes = NULL)

define_f.gstar <- function(shift.val, truncBD, rndseed = NULL) {
  shift.const <- shift.val
  trunc.const <- truncBD
  f.gstar <- function(data, ...) {
    print("shift.const: " %+% shift.const)
    set.seed(rndseed)
    A.mu <- 0.86 * data[,"W1"] + 0.41 * data[,"W2"] - 0.34 * data[,"W3"] + 0.93 * data[,"W4"]
    untrunc.A <- rnorm(n = nrow(data), mean = A.mu + shift.const, sd = 1)
    r.new.A <- exp(0.5 * shift.const * (untrunc.A - A.mu - shift.const / 2))
    trunc.A <- ifelse(r.new.A > trunc.const, untrunc.A - shift.const, untrunc.A)
    return(trunc.A)
  }
  return(f.gstar)
}
f.gstar <- define_f.gstar(shift = sampleDat_iidcontAContY$shift.val, truncBD = sampleDat_iidcontAContY$truncBD, 
                          rndseed = sampleDat_iidcontAContY$rndseed)

tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidcontAContY))
OData <- DatKeepClass$new(Odata = dat_iidcontAContY, nodes = nodes, norm.c.sVars = FALSE)  # don't normalize continous covariates here
OData$names.sVar  # "Y"  "A"  "W1" "W2" "W3" "W4" (names of all variables in input data )
OData$names.c.sVar  # "A"  "W3" "W4" (names of all continuous variables in input data)
head(OData$dat.sVar)  # a subset data.frame of Odata that includes all variables in nodes
nobs <- OData$nobs; nobs # 10000

obsYvals <- dat_iidcontAContY[, nodes$Ynode]
OData$addYnode(YnodeVals = obsYvals)  # The same as OData$addYnode(YnodeVals = obsYvals, det.Y = F)
head(OData$YnodeVals)  # Adding public YnodeVals & setting det.Y values to NA
head(OData$noNA.Ynodevals)  # Adding actual observed Y as protected (without NAs)
OData$addObsWeights(obs.wts = rep(c(1,2), 5000))  # Add a vectopr of observation (sampling) weights
OData$addObsWeights(obs.wts = 1)  # If not specified, assumed to be all 1 (i.e., equally weighted)

OData$get.sVar.type("A")  # "contin"
OData$get.sVar.type()  # Provide a list of types of all variables

OData.gstar <- DatKeepClass$new(Odata = dat_iidcontAContY, nodes = nodes, norm.c.sVars = FALSE)
# Create (p=) 1 new Odata (nobs obs at a time) and replace A under g0 in Odata with A^* under g.star
OData.gstar$make.dat.sVar(p = 1, f.g_fun = f.gstar)  # Generate new exposures under user-specific arbitrary intervention f.g_fun.
dim(OData.gstar$dat.sVar)  # 10000     5
OData.gstar$make.dat.sVar(p = 3, f.g_fun = f.gstar) 
dim(OData.gstar$dat.sVar)  # 30000     5
