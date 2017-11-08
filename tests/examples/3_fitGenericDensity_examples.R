data(indSample.iid.cA.bY_list)
indSample.iid.cA.bY <- indSample.iid.cA.bY_list$indSample.iid.cA.bY
tmleCom_Options(gestimator = "speedglm__glm", maxNperBin = nrow(indSample.iid.cA.bY))
gvars$verbose <- TRUE  # Print status messages (global setting)

# Define a stochastic intervention
define_f.gstar <- function(data, shift.rate, ...) {
  print("rate of shift: " %+% shift.rate)
  shifted.new.A <- data[, "A"] * shift.rate
  return(shifted.new.A)
}
f.gstar <- define_f.gstar(shift.rate = 0.5)

# Run
h_gN <- fitGenericDensity(data = indSample.iid.cA.bY, Anodes = "A", 
                          Wnodes = c("W1", "W2", "W3", "W4"), 
                          f_gstar = NULL, lbound = 0)$h_gstar
h_gstar <- fitGenericDensity(data = indSample.iid.cA.bY, Anodes = "A",
                             Wnodes = c("W1", "W2", "W3", "W4"), 
                             f_gstar = f.gstar, lbound = 0)$h_gstar
