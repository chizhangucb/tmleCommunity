data(sampleDat_iidcontAContY)
dat_iidcontAContY <- sampleDat_iidcontAContY$dat_iidcontAContY
nodes <- list(Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), Enodes = NULL, Crossnodes = NULL)

# Parse the formulae for actual covariate names in A, W & E
# If the formula is not specified (i.e., NULL), variable(s) input into Anodes.lst will be treated as outcome variable(s),
# and all variables input into Wnodes.lst will be used as explanatory variable(s).
Q.sVars <- define_regform(NULL, Anodes.lst = nodes$Ynode, Wnodes.lst = nodes[c("Anodes", "Wnodes", "Enodes")])
tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidcontAContY))
Qreg <- RegressionClass$new(outvar = Q.sVars$outvars, predvars = Q.sVars$predvars, subset_vars = (!rep_len(FALSE, nobs)))
Qreg$estimator  # "speedglm__glm"
