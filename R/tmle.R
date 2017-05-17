tmle <- function(Y,A,W,Z=NULL, Delta=rep(1,length(Y)),  
                 Q=NULL, Q.Z1=NULL, Qform=NULL, Qbounds=NULL, 
                 Q.SL.library=c("SL.glm", "SL.step", "SL.glm.interaction"), cvQinit=FALSE,
                 g1W=NULL, gform=NULL, gbound=0.025, 
                 pZ1=NULL, g.Zform=NULL,
                 pDelta1=NULL, g.Deltaform=NULL, 
                 g.SL.library=c("SL.glm", "SL.step", "SL.glm.interaction"),
                 family="gaussian", fluctuation="logistic", 
                 alpha  = 0.995, id=1:length(Y), V = 5, verbose=FALSE) {
  # Initializations
  psi.tmle <- varIC <- CI <- pvalue <- NA
  # W <- as.matrix(W)  # if W is a dataframe with factors, changing to matrix will blow it
  colnames(W) <- .setColnames(colnames(W), NCOL(W), "W")
  if (identical(family, binomial)) {
    family == "binomial"
  }else if (identical(family, gaussian)) {
    family == "gaussian"
  }else if (identical(family, poisson)) {  
    family == "poisson"
    if(is.null(Qform)){ 		# only glm for Q when family=poisson
      Qform <- paste("Y~A", paste(colnames(W), collapse="+"), sep="+")
    }
  }
  if(is.null(A) | all(A==0)){
    A <- rep(1, length(Y))
  }
  
  if(!.verifyArgs(Y,Z,A,W,Delta, Qform, gform, g.Zform, g.Deltaform)){
    stop()
  }
  
  maptoYstar <- fluctuation=="logistic"
  
  if(!is.null(Z) & !is.null(pZ1)){
    if(NCOL(pZ1)==2){
      colnames(pZ1) <- c("A0", "A1")
    } else {
      stop("pZ1 must be an nx2 matrix: [P(Z=1|A=0,W),P(Z=1|A=1,W)]\n")
    }
  }
  
  if(any(Delta!=1) & !is.null(pDelta1)){
    if(NCOL(pDelta1)==2 & is.null(Z)){
      pDelta1 <- cbind(pDelta1, pDelta1)
      colnames(pDelta1) <- c("Z0A0", "Z0A1", "Z1A0", "Z1A1")  # won't use the second set
    } else if(NCOL(pDelta1)==4){
      colnames(pDelta1) <- c("Z0A0", "Z0A1", "Z1A0", "Z1A1") 
    }else {
      if(is.null(Z)){
        stop("pDelta1 must be an nx2 matrix: [P(Delta=1|A=0,W), P(Delta=1|A=1,W)]\n")
      } else {
        stop("pDelta1 must be an nx4 matrix:\n  [P(Delta=1|Z=0,A=0,W), P(Delta=1|Z=0,A=1,W), P(Delta=1|Z=1,A=0,W), P(Delta=1|Z=1,A=1,W)]\n")
      }
    }
  }
  
  if(is.null(Z)){ 
    Z <- rep(1, length(Y))	
  }
  CDE <- length(unique(Z))>1
  
  # Stage 1	
  stage1 <- .initStage1(Y, A, Q, Q.Z1, Delta, Qbounds, alpha, maptoYstar, family)		
  Q <- suppressWarnings(estimateQ(Y=stage1$Ystar,Z,A,W, Delta, Q=stage1$Q, Qbounds=stage1$Qbounds, Qform, 
                                  maptoYstar=maptoYstar, SL.library=Q.SL.library, 
                                  cvQinit=cvQinit, family=family, id=id, V = V, verbose=verbose))
  
  # Stage 2
  if(length(gbound)==1){
    if(length(unique(A))==1 & length(unique(Z))==1){  # EY1 only, no controlled direct effect
      gbound <- c(gbound,1)
    } else {
      gbound <- c(gbound, 1-gbound)
    }
  }
  g <- suppressWarnings(estimateG(d=data.frame(A,W), g1W, gform, g.SL.library, id=id, V = V, verbose, "treatment mechanism", outcome="A")) 
  g$bound <- gbound
  if(g$type=="try-error"){
    stop("Error estimating treatment mechanism (hint: only numeric variables are allowed)") 
  }
  
  if(!CDE){
    g.z <- NULL
    g.z$type="No intermediate variable"
    g.z$coef=NA
    g.Delta <- suppressWarnings(estimateG(d=data.frame(Delta, Z=1, A, W), pDelta1, g.Deltaform, 
                                          g.SL.library,id=id, V = V, verbose = verbose, "missingness mechanism", outcome="D")) 
    g1W.total <- .bound(g$g1W*g.Delta$g1W[,"Z0A1"], gbound)
    g0W.total <- .bound((1-g$g1W)*g.Delta$g1W[,"Z0A0"], gbound)  
    if(all(g1W.total==0)){g1W.total <- rep(10^-9, length(g1W.total))}
    if(all(g0W.total==0)){g0W.total <- rep(10^-9, length(g0W.total))}
    H1W <- A/g1W.total
    H0W <- (1-A)/g0W.total
    
    suppressWarnings(
      epsilon <- coef(glm(stage1$Ystar~-1 + offset(Q$Q[,"QAW"]) + H0W + H1W, family=Q$family, subset=Delta==1))
    )
    epsilon[is.na(epsilon)] <- 0  # needed for EY1 calculation
    Qstar <- Q$Q + c((epsilon[1]*H0W + epsilon[2]*H1W), epsilon[1]/g0W.total, epsilon[2]/g1W.total)
    colnames(Qstar) <- c("QAW", "Q0W", "Q1W")
    Ystar <- stage1$Ystar
    if (maptoYstar) {
      Qstar <- plogis(Qstar)*diff(stage1$ab)+stage1$ab[1] 
      Q$Q <- plogis(Q$Q)*diff(stage1$ab)+stage1$ab[1]
      Ystar <- Ystar*diff(stage1$ab)+stage1$ab[1]
    } else if (family == "poisson"){  	
      Q$Q <- exp(Q$Q)				  
      Qstar <- exp(Qstar)
    }
    colnames(Q$Q) <- c("QAW", "Q0W", "Q1W")
    Q$Q <- Q$Q[,-1]
    res <- calcParameters(Ystar, A, I.Z=rep(1, length(Ystar)), Delta, g1W.total, g0W.total, Qstar, 
                          mu1=mean(Qstar[,"Q1W"]), mu0=mean(Qstar[,"Q0W"]), id, family)
    returnVal <- list(estimates=res, Qinit=Q, g=g, g.Z=g.z, g.Delta=g.Delta, Qstar=Qstar[,-1], epsilon=epsilon) 
    class(returnVal) <- "tmle"
  } else {
    returnVal <- vector(mode="list", length=2)
    g.z <- suppressWarnings(estimateG(d=data.frame(Z,A,W), pZ1, g.Zform, g.SL.library, id=id, V = V, 
                                      verbose, "intermediate variable", outcome="Z"))
    g.Delta <- suppressWarnings(estimateG(d=data.frame(Delta,Z, A, W), pDelta1, g.Deltaform, 
                                          g.SL.library,id=id, V=V, verbose, "missingness mechanism", outcome="D")) 
    ZAD <- cbind(D1Z0A0 = .bound((1-g$g1W)*(1-g.z$g1W[,"A0"])*g.Delta$g1W[,"Z0A0"], gbound),
                 D1Z0A1 = .bound(g$g1W*(1-g.z$g1W[,"A1"])*g.Delta$g1W[,"Z0A1"], gbound),
                 D1Z1A0 = .bound((1-g$g1W)*g.z$g1W[,"A0"]*g.Delta$g1W[,"Z1A0"], gbound),
                 D1Z1A1 = .bound(g$g1W*g.z$g1W[,"A1"]*g.Delta$g1W[,"Z1A1"], gbound))
    adjustZero <- colSums(ZAD)==0
    ZAD[,adjustZero] <- 10^-9
    for (z in 0:1){
      H0W <- (1-A)*(Z==z)/ZAD[,z*2+1]
      H1W <- A*(Z==z) /ZAD[,z*2+2]
      suppressWarnings(
        epsilon <- coef(glm(stage1$Ystar~-1 + offset(Q$Q[,"QAW"]) + H0W + H1W, family=Q$family,
                            subset=(Delta==1 & Z==z)))
      )  			
      
      hCounter <- cbind(1/ZAD[,z*2+1], 1/ZAD[,z*2+2])
      Qstar <- Q$Q[,c(1, z*2+2, z*2+3)] + c((epsilon[1]*H0W + epsilon[2]*H1W), 
                                            epsilon[1]*hCounter[,1], epsilon[2]*hCounter[,2])
      colnames(Qstar) <- c("QAW", "Q0W", "Q1W") 
      newYstar <- stage1$Ystar     
      if (maptoYstar) {
        Qstar <- plogis(Qstar)*diff(stage1$ab)+stage1$ab[1] 
        Qinit.return <- plogis(Q$Q[,c(1, z*2+2, z*2+3)])*diff(stage1$ab)+stage1$ab[1]
        newYstar <- stage1$Ystar*diff(stage1$ab)+stage1$ab[1]
      } else if (family == "poisson"){ 
        Qinit.return <- exp(Q$Q[,c(1, z*2+2, z*2+3)]) 
        Qstar <- exp(Qstar)
      }
      colnames(Qinit.return) <- c("QAW", "Q0W", "Q1W")
      res <- calcParameters(newYstar, A,I.Z=as.integer(Z==z), Delta, g1W=ZAD[,z*2+2], g0W=ZAD[,z*2+1],  
                            Qstar, mu1=mean(Qstar[,"Q1W"]), mu0=mean(Qstar[,"Q0W"]), id, family)
      Qreturn <- Q
      Qreturn$Q <- Qinit.return[,-1]
      returnVal[[z+1]] <- list(estimates=res, Qinit=Qreturn, g=g, g.Z=g.z, g.Delta=g.Delta, 
                               Qstar=Qstar[,-1], epsilon=epsilon)
    }
    class(returnVal[[1]]) <- class(returnVal[[2]]) <- "tmle"
    class(returnVal) <- "tmle.list"
  }
  return(returnVal)
}

