########### General utilities / Global Vars
`%+%` <- function(a, b) paste0(a, b)
is.VecMatDf <- function(x) { is.vector(x) || is.matrix(x) || is.data.frame(x) }

#------------------ function CreateInputs ------------------
# Purpose: Bound Y by mapping to Ystar if applicable, and
# set bounds on on Q and enforce on user-specified values
# returns:
#    Ystar - outcome values (between [0,1] if maptoYstar=TRUE)
#    Q - matrix of user-specified values
#    Qbounds - bounds on predicted values for Q (10% wider at each end than observed range of Y)
#       	 		 (-Inf,+Inf) is default for linear regression
#    ab - bounding levels used to transform Y to Ystar
#----------------------------------------------------------
CreateInputs <- function(Y, Qbounds, alpha, maptoYstar){
  if (all(Y >= 0 & Y <= 1)) Qbounds <- c(0, 1)
  if (is.null(Qbounds)) {
    if (maptoYstar) { 
      # 10% wider at each end than observed range of Y
      Qbounds <- range(Y) 
      Qbounds <- Qbounds + 0.1*c(-abs(Qbounds[1]), abs(Qbounds[2]))
    } else {
      Qbounds <- c(-Inf, Inf)  # for linear regression
    }
  }
  ab <- c(0, 1)  # Default
  Ystar <- Y
  if (maptoYstar) { 
    Ystar <- bound(Y, Qbounds)  # Bound Y by Qbounds
    if (0 >= alpha | 1 <= alpha) {
      alpha <- .995
      warning(paste("\n\talpha must be between 0 and 1, alpha reset to",alpha,"\n"), immediate. = TRUE)
    }
    ab <- range(Ystar, na.rm=TRUE)
    Ystar[is.na(Ystar)] <- 0  # For missing outcomes, treat as 0
    Ystar <- (Ystar-ab[1]) / diff(ab)  # Bound Ystar into [0, 1]. Transformed by (Y-a)/(b-a)
    Qbounds <- c(1 - alpha, alpha)
  }	
  return(list(Ystar=Ystar, Qbounds=Qbounds, ab=ab, maptoYstar = maptoYstar))
} 


#------------------ function CheckInputs ------------------
# Purpose: initial checks on data passed in
#----------------------------------------------------------
CheckInputs <- function(data, nodes, Qform, hform.g0, hform.gstar, fluctuation, Qbounds, obs.wts, community.wts, f_gstar1, f_gstar2) {
  datfactor <- CheckExistFactors(data)
  formulas <- list(Qform, hform.g0, hform.gstar) 
  validFormula <- sapply(formulas, function(x) {identical(class(try(as.formula(x))), "formula")})
  validNames <- c(as.vector(unlist(nodes)), ".")
  validTerms <- rep(TRUE, length(formulas))
  validTerms[validFormula] <- sapply(formulas[which(validFormula)], function(x) {
    is.null(x) || all(all.names(as.formula(x), functions=FALSE) %in% validNames)})
  validFluct <- fluctuation %in% c("logistic", "linear")
  formwarns <- c("\tInvalid regression formula for 'Qform' " %+% deparse(formulas[[1]]),
                 "\tInvalid regression formula for 'hform.g0' " %+% deparse(formulas[[2]]),
                 "\tInvalid regression formula for 'hform.gstar' " %+% deparse(formulas[[3]]))
  termwarns <- c("\tInvalid term name in regression formula for 'Qform' " %+% deparse(formulas[[1]]),
                 "\tInvalid term name in regression formula for 'hform.g0' " %+% deparse(formulas[[2]]),
                 "\tInvalid term name in regression formula for 'hform.gstar' " %+% deparse(formulas[[3]]))
  f_gstar1.part1 <- is.null(f_gstar1) || is.function(f_gstar1) || is.VecMatDf(f_gstar1)
  if (!is.null(f_gstar1) && !is.function(f_gstar1)) {
    f_gstar1 <- as.data.frame(f_gstar1)
    f_gstar1.part2 <- NCOL(f_gstar1)==length(nodes$Anodes) && (NROW(f_gstar1)==NROW(data) || NROW(f_gstar1)==1)
  } else {
    f_gstar1.part2 <- TRUE
  }
  f_gstar2.part1 <- is.null(f_gstar2) || is.function(f_gstar2) || is.VecMatDf(f_gstar2)
  if (!is.null(f_gstar2) && !is.function(f_gstar2)) {
    f_gstar2 <- as.data.frame(f_gstar2)
    f_gstar2.part2 <- NCOL(f_gstar2)==length(nodes$Anodes) && (NROW(f_gstar2)==NROW(data) || NROW(f_gstar2)==1)
  } else {
    f_gstar2.part2 <- TRUE
  }
  pass <- c(is.data.frame(data), is.null(datfactor), is.null(obs.wts) || (length(obs.wts)==NROW(data) && all(obs.wts >= 0)),
            is.null(community.wts) || (all(dim(community.wts)==c(length(unique(data[, nodes$communityID])), 2) && community.wts[, 2] >= 0)),
            is.null(community.wts) || (setequal(sort(community.wts[,1]), sort(unique(data[, nodes$communityID]))) && 
                                       !anyDuplicated(community.wts[,1])),
            f_gstar1.part1, f_gstar1.part2, f_gstar2.part1, f_gstar2.part2,
            validFormula, validTerms, validFluct, is.null(Qbounds) || length(Qbounds)==2)
  warning_messages <- c("\tThe input data must be a data frame",
                        "\tNo factor column(s) allowed in the input data, consider removing or recoding such column(s) as strings: " 
                        %+% paste0(datfactor, collapse=' , ') %+% "\n", 
                        "\t'obs.wts', must contain the same number of non-negative observations as 'data' does\n",
                        "\t'community.wts', must contain the same number of non-negative communities as 'data' does, and has 2 columns\n",
                        "\t'community.wts', must contain the same community IDs as data doee (diff order is ok but no duplicate allows\n",
                        "\t'f_gstar1', should be either a function, or a vector, or a matrix/ data frame\n",
                        "\tIf 'f_gstar1' is a vector/matrix/data.frame, it should have a length (number of rows) NROW(data)\n",
                        "\t'f_gstar2', should be either a function, or a vector, or a matrix/ data frame\n",
                        "\tIf 'f_gstar2' is a vector/matrix/data.frame, it should have a length (number of rows) NROW(data)\n",
                        formwarns, termwarns, "\tfluctuation should be logistic or linear\n",
                        "\tQbounds should have length 2\n")
  if(!all(pass)) warning("\n", warning_messages[!pass], immediate.=TRUE)
  return(all(pass))
}


#------------------------------------ CheckExistFactors -------------------------------------
# Purpose: Returns NULL if no factors exist, otherwise return the name of the factor variable(s)
#--------------------------------------------------------------------------------------------
CheckExistFactors <- function(data) {
  testvec <- unlist(lapply(data, is.factor))
  if (any(testvec)) {
    return(names(data)[which(testvec)])
  } else {
    return(NULL)
  }
}


#------------------------------- function SuppressGivenWarnings -------------------------------
# Purpose: if warning is in ignoreWarningList, ignore it; otherwise post it as usual
#---------------------------------c------------------------------------------------------------
SuppressGivenWarnings <- function(expr, warningsToIgnore) {
  h <- function (w) {
    if (w$message %in% warningsToIgnore) invokeRestart( "muffleWarning" )
  }
  withCallingHandlers(expr, warning = h )
}


#------------------------------------- GetWarningsToSuppress ----------------------------------
# Purpose: create an ignoreWarningList
#---------------------------------c------------------------------------------------------------
GetWarningsToSuppress <- function(update.step=FALSE) {
  warnings.to.suppress <- c("glm.fit: fitted probabilities numerically 0 or 1 occurred",
                            "prediction from a rank-deficient fit may be misleading")
  if (update.step) {
    warnings.to.suppress <- c(warnings.to.suppress, "glm.fit: algorithm did not converge")
  }
  return(warnings.to.suppress)
}


#------------------------------------------ bound -------------------------------------------
# Purpose: set outliers to min/max allowable values, assuming x contains only numerical data
#--------------------------------------------------------------------------------------------
bound <- function(x, bounds){
  x[x<min(bounds)] <- min(bounds)
  x[x>max(bounds)] <- max(bounds)
  return(x)
}

#---------------------------------------- LhsVars -------------------------------------------
# Purpose: Extract the variables on the left side of a formula 
#--------------------------------------------------------------------------------------------
LhsVars <- function(f) {
  f <- as.formula(f)
  return(as.character(f[[2]]))
}

#------------------------------------ CheckVarNameExists ------------------------------------
# Purpose: throw exception if varname doesn't exist
#--------------------------------------------------------------------------------------------
CheckVarNameExists <- function(data, varname) {
  idvar <- names(data) %in% varname
  if (sum(idvar) < 1) stop("variable name " %+% varname %+% " not found in data input")
  if (sum(idvar) > 1) stop("more than one column in the input data has been matched to name " 
                            %+% varname %+% ". Consider renaming some of the columns: " %+% 
                            paste0(names(data)[idvar], collapse=","))
  return(invisible(NULL))
}
