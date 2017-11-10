data(indSample.iid.cA.cY_list)
indSample.iid.cA.cY <- indSample.iid.cA.cY_list$indSample.iid.cA.cY
tmleCom_Options(gestimator = "speedglm__glm", maxNperBin = nrow(indSample.iid.cA.cY),
                bin.method = "dhist", nbins = 8)
gvars$verbose <- TRUE  # Print status messages (global setting)

# Define a stochastic intervention
define_f.gstar <- function(shift.rate, ...) {
  eval(shift.rate)
  f.gstar <- function(data, ...) {
    print("rate of shift: " %+% shift.rate)
    shifted.new.A <- data[, "A"] - mean(data[, "A"]) * shift.rate
    return(shifted.new.A)
  }
  return(f.gstar)
}
f.gstar <- define_f.gstar(shift.rate = 0.5)

# Under current treatment mechanism g0
h_gN <- fitGenericDensity(data = indSample.iid.cA.cY, Anodes = "A", 
                          Wnodes = c("W1", "W2", "W3", "W4"), 
                          f_gstar = NULL, lbound = 0)$h_gstar
# Under stochastic intervention gstar
h_gstar <- fitGenericDensity(data = indSample.iid.cA.cY, Anodes = "A",
                             Wnodes = c("W1", "W2", "W3", "W4"), 
                             f_gstar = f.gstar, lbound = 0)$h_gstar
