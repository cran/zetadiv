###########################################################################################################################################
##    THESE FUNCTIONS ARE COPIED FROM THE PACKAGE SCAM V.1.2.0, BY NATALYA PYA , BECAUSE THEY ARE NOT EXPORTED FROM THE SCAM PACKAGE     ##
##    Natalya Pya (2016). scam: Shape Constrained Additive Models. R package version 1.2-0. https://CRAN.R-project.org/package=scam      ##
###########################################################################################################################################

#################
## fix for scam##
##This function is a copy of scam() from the package scam V.1.2.0 by Natalya Pya, with a slight modification to avoid errors when the initial smoothing terms are 0 (see the comments in the function.
#################

Scam.zeta <- function (formula, family = stats::gaussian(), data = list(), gamma = 1, 
                       sp = NULL, weights = NULL, offset = NULL, optimizer = "bfgs", 
                       optim.method = c("Nelder-Mead", "fd"), scale = 0, devtol = 1e-08, 
                       steptol = 1e-08, check.analytical = FALSE, del = 1e-04, start = NULL, 
                       etastart, mustart, keepData = FALSE, not.exp = FALSE) 
{
  G <- mgcv::gam(formula, family, data, fit = FALSE)
  n.terms <- length(G$smooth)
  n <- nrow(G$X)
  intercept <- G$intercept
  gp <- mgcv::interpret.gam(formula)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$formula <- gp$fake.formula
  mf$family <- mf$control <- mf$scale <- mf$knots <- mf$sp <- mf$min.sp <- mf$H <- mf$select <- mf$gamma <- mf$method <- mf$fit <- mf$paraPen <- mf$G <- mf$optimizer <- mf$optim.method <- mf$not.exp <- mf$in.out <- mf$... <- NULL
  mf[[1]] <- as.name("model.frame")
  pmf <- mf
  mf <- eval(mf, parent.frame())
  G$offset <- as.vector(stats::model.offset(mf))
  if (is.null(G$offset)) 
    G$offset <- rep.int(0, n)
  if (is.null(weights)) 
    weights <- rep.int(1, n)
  fam.name <- G$family[1]
  if (scale == 0) {
    if (fam.name == "binomial" || fam.name == "poisson") 
      sig2 <- 1
    else sig2 <- -1
  }
  else {
    sig2 <- scale
  }
  if (sig2 > 0) 
    scale.known <- TRUE
  else scale.known <- FALSE
  Q <- penalty_pident(G)
  if (!is.null(sp)) {
    neg <- FALSE
    if (length(sp) != length(G$off)) {
      warning("Supplied smoothing parameter vector is too short - ignored.")
      sp <- NULL
    }
    else if (sum(is.na(sp))) {
      warning("NA's in supplied smoothing parameter vector - ignoring.")
      sp <- NULL
    }
    else {
      good <- sp < 0
      if (sum(good) > 0) {
        warning("Supplied smoothing parameter vector has negative values - ignored.")
        neg <- TRUE
      }
    }
    if (neg) 
      sp <- NULL
  }
  env <- new.env()
  assign("start", rep(0, 0), envir = env)
  assign("dbeta.start", rep(0, 0), envir = env)
  assign("sp.last", rep(0, 0), envir = env)
  q.f <- rep(0, n.terms)
  for (i in 1:n.terms) {
    q.f[i] <- ncol(G$smooth[[i]]$S[[1]]) + 1
  }
  G$S <- Q$S
  G$q.f <- q.f
  G$q0 <- G$off[1] - 1
  G$p.ident <- Q$p.ident
  G$n.terms <- n.terms
  G$weights <- weights
  G$sig2 <- sig2
  G$scale.known <- scale.known
  G$not.exp <- not.exp
  if (!keepData) 
    rm(data)
  object <- list()
  if (is.null(sp)) {
    start <- etastart <- mustart <- NULL
    y <- G$y
    family <- G$family
    nobs <- NROW(y)
    eval(family$initialize)
    G$y <- y
    def.sp <- initial.sp.scam(G, Q, q.f = q.f, n.terms = n.terms, 
                              family = family, intercept = intercept, offset = G$offset, 
                              env = env, weights = weights, devtol = 1e-04, steptol = 1e-04)
    rho <- log(def.sp+0.001)                               ##THIS WAS MODIFIED FROM THE ORIGINAL FUNCTION SCAM FROM THE SCAM PACKAGE V.1.2.0 TO AVOID LOG(0)
    ptm <- proc.time()
    re <- estimate.scam(G = G, optimizer = optimizer, optim.method = optim.method, 
                        rho = rho, gamma = gamma, env = env, check.analytical = check.analytical, 
                        del = del, devtol = devtol, steptol = steptol)
    CPU.time <- proc.time() - ptm
    best <- re
    object$gcv.ubre <- re$gcv.ubre
    object$dgcv.ubre <- re$dgcv.ubre
    object$aic <- re$aic
    best$p.ident <- Q$p.ident
    best$S <- Q$S
    object$optimizer <- optimizer
    object$edf1 <- re$edf1
    object$termcode <- re$termcode
    if (optimizer == "bfgs") {
      object$check.grad <- re$check.grad
      object$dgcv.ubre.check <- re$dgcv.ubre.check
    }
  }
  else {
    best <- scam::scam.fit(G = G, sp = sp, gamma = gamma, devtol = devtol, 
                           steptol = steptol, env = env)
    object$aic <- best$aic
    object$optimizer <- "NA"
  }
  best$n.smooth <- object$n.smooth <- n.terms
  best$formula <- object$formula <- formula
  best$family <- object$family <- G$family
  best$smooth <- object$smooth <- G$smooth
  best$model <- object$model <- G$mf
  object$R <- best$R
  if (is.null(object$R)) {
    rr <- scam::scam.fit(G = G, sp = best$sp, gamma = gamma, devtol = devtol, 
                         steptol = steptol, env = env)
    object$R <- rr$R
  }
  object$df.residual <- nrow(best$X) - sum(best$edf)
  object$sp <- best$sp
  names(object$sp) <- names(G$sp)
  if (sum(is.na(names(object$sp))) != 0) {
    for (i in 1:n.terms) names(object$sp)[i] <- object$smooth[[i]]$label
  }
  object$deviance <- best$deviance
  object$residuals <- best$residuals
  object$conv <- best$conv
  post <- scam.fit.post(y = G$y, X = G$X, object = best, sig2 = sig2, 
                        offset = G$offset, intercept = G$intercept, weights = weights, 
                        scale.known = scale.known)
  object$edf <- post$edf
  object$edf1 <- post$edf1
  object$trA <- post$trA
  names(object$edf) <- G$term.names
  names(object$edf1) <- G$term.names
  object$null.deviance <- post$nulldev
  object$var.summary <- G$var.summary
  object$cmX <- G$cmX
  object$model <- G$mf
  object$full.sp <- G$full.sp
  if (!is.null(object$full.sp)) 
    names(object$full.sp) <- names(G$full.sp)
  object$na.action <- attr(G$mf, "na.action")
  object$df.null <- post$df.null
  object$Ve <- post$Ve
  object$Vp <- post$Vb
  object$Ve.t <- post$Ve.t
  object$Vp.t <- post$Vb.t
  object$sig2 <- post$sig2
  object$coefficients <- best$beta
  object$coefficients.t <- best$beta.t
  object$beta <- best$beta
  object$beta.t <- best$beta.t
  object$pterms <- G$pterms
  object$terms <- G$terms
  object$assign <- G$assign
  object$nsdf <- G$nsdf
  object$y <- G$y
  if (keepData) 
    object$data <- data
  object$offset <- G$offset
  object$not.exp <- G$not.exp
  object$scale.estimated <- !scale.known
  object$prior.weights <- weights
  object$weights <- best$w
  object$fitted.values <- best$mu
  object$linear.predictors <- best$eta
  object$call <- cl
  object$p.ident <- Q$p.ident
  object$intercept <- G$intercept
  object$min.edf <- G$min.edf
  object$gamma <- gamma
  object$iter <- best$iter
  if (is.null(sp)) 
    object$CPU.time <- CPU.time
  else object$CPU.time <- NULL
  if (is.null(sp)) {
    if (optimizer == "bfgs") {
      object$bfgs.info <- list()
      object$bfgs.info$conv <- re$conv.bfgs
      object$bfgs.info$iter <- re$iterations
      object$bfgs.info$grad <- re$dgcv.ubre
    }
    else if (optimizer == "nlm.fd" || optimizer == "nlm") {
      object$nlm.info <- list()
      object$nlm.info$conv <- re$conv
      object$nlm.info$iter <- re$iterations
      object$nlm.info$grad <- re$dgcv.ubre
    }
    else if (optimizer == "optim") {
      object$optim.info <- list()
      object$optim.info$conv <- re$conv
      object$optim.info$iter <- re$iterations
      object$optim.method <- re$optim.method
    }
  }
  if (scale.known) 
    object$method <- "UBRE"
  else object$method <- "GCV"
  if (G$nsdf > 0) 
    term.names <- colnames(G$X)[1:G$nsdf]
  else term.names <- array("", 0)
  if (n.terms) 
    for (i in 1:n.terms) {
      k <- 1
      for (j in G$smooth[[i]]$first.para:G$smooth[[i]]$last.para) {
        term.names[j] <- paste(G$smooth[[i]]$label, ".", 
                               as.character(k), sep = "")
        k <- k + 1
      }
    }
  names(object$coefficients) <- term.names
  names(object$coefficients.t) <- term.names
  ynames <- if (is.matrix(G$y)) 
    rownames(G$y)
  else names(G$y)
  names(object$residuals) <- ynames
  class(object) <- c("scam", "glm", "lm")
  object
}

#################################################################
## function to get initial estimates of smoothing parameters...##
#################################################################

initial.sp.scam <- function(G,Q,q.f,n.terms,family,intercept,offset, env= env,
                            weights,devtol=1e-4,steptol=1e-4) 
{  ## function to get initial estimates of smoothing parameters
  ## step 1: set sp=rep(0.5,p) and estimate hessian...
  b <- scam::scam.fit(G=G,sp=rep(0.5,length(G$off)), devtol, steptol, env=env) 
  H <- crossprod(b$wX1) - b$E
  ## step 2:...
  n.p <- length(Q$S) ## number of penalty matrices
  def.sp <- array(0,n.p) ## initialize the initial sp values
  j <- 1
  for (i in 1:n.terms)
  {   for (kk in 1:length(G$smooth[[i]]$S))
  {   start <- G$off[j]
  finish <- start + ncol(G$smooth[[i]]$S[[kk]])-1
  # matrix norm of the Hessian elements penalized by S[[kk]]...
  Hi.norm <- sum(H[start:finish,start:finish]*H[start:finish,start:finish]) 
  Si.norm <- sum(G$smooth[[i]]$S[[kk]]*G$smooth[[i]]$S[[kk]])
  def.sp[j] <- (Hi.norm/Si.norm)^0.5
  j <- j+1
  }
  }
  ## Create again new environments with `start' initially empty...
  env <- new.env()
  assign("start",rep(0,0),envir=env)
  assign("dbeta.start",rep(0,0),envir=env)
  assign("sp.last",rep(0,0),envir=env)
  def.sp
}


#########################################################
## function to get list of penalty matrices and        ## 
## vector of parameter identifications .....           ##
#########################################################


penalty_pident <- function(object)
{  ## function to get the list of penalties and vector of model parameters 
  ## identifications from the gam() setting...
  n.terms <- length(object$smooth)  # number of terms in the model
  q <- ncol(object$X)          # total number of parameters
  cons.terms <- rep(0,n.terms) # define whether each term is constrained or not
  for (i in 1:n.terms)
  {   if (!is.null(object$smooth[[i]]$p.ident))
    cons.terms[i] <- 1  
  }
  p.ident <- rep(0,q) # initialize vector of parameter identifications
  # with `1' - for a parameter to be exponentiated, `0' - otherwise
  off.terms <- rep(0,n.terms) # starting points for each term
  off <- object$off
  if (n.terms ==length(off))
    off.terms <- off
  else 
  {   off.terms[1] <- off[1]
  k <- 1
  l <- 1
  while (l<length(off))
  {   if (off[l]!=off[l+1])
  {   off.terms[k+1] <- off[l+1] 
  k <- k+1; l <- l+1 
  } 
    else l <- l+1
  }
  
  }
  for (i in 1:n.terms)
  {   if (cons.terms[i]==1) 
    p.ident[off.terms[i]:(off.terms[i]+ncol(object$smooth[[i]]$S[[1]])-1)] <- 
      object$smooth[[i]]$p.ident
  }
  ## getting the list of penalty matrices in terms of the full model vector of coefficients...
  S <- list()
  j <- 1
  for(i in 1:n.terms)
  { for (kk in 1:length(object$smooth[[i]]$S))
  {    S[[j]] <- matrix(0,q,q) # initialize penalty matrix
  S[[j]][off.terms[i]:(off.terms[i]+ncol(object$smooth[[i]]$S[[kk]])-1),
         off.terms[i]:(off.terms[i]+ncol(object$smooth[[i]]$S[[kk]])-1)] <- object$smooth[[i]]$S[[kk]]
  j <- j+1       
  }
  }
  object$S <- S 
  object$p.ident <- p.ident
  object
} ## end penalty_pident



#######################################################################
## function to get null deviance and covariance matrices after fit   ##
#######################################################################


scam.fit.post <- function(y,X,object,sig2,offset,intercept,
                          weights,scale.known)
{  ## Function to compute null deviance and covariance matrices after a scam fit.
  ## covariance matrix should use expected Hessian, so re-computation of factors 
  ## is required.  
  ## object - object from estimate.scam()
  n <- nobs <- NROW(y) # number of observations
  linkinv <- object$family$linkinv
  dev.resids <- object$family$dev.resids
  
  wtdmu <- if (intercept) sum(weights * y)/sum(weights) 
  else linkinv(offset)
  
  nulldev <- sum(dev.resids(y, wtdmu, weights))
  
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - as.integer(intercept)
  
  # calculating the approximate covariance matrices 
  # (dealing with the expected Hessian of the log likelihood) ...
  
  if (!scale.known) sig2 <- object$scale.est
  ## get the inverse of the expected Hessian...
  ## wX1 <- sqrt(object$w1)[1:n]*object$X1
  wX11 <- rbind(object$wX1,object$rS)
  q <- ncol(object$wX1)
  Q <- qr(wX11,LAPACK=TRUE) 
  R <- qr.R(Q)
  rp <- 1:ncol(R)
  rp[Q$pivot] <- rp ## reverse pivot X=Q%*%R[,rp]
  if (mgcv::Rrank(R)==ncol(R))  ## no need to truncate, can just use QR
  {  P <- backsolve(R,diag(q))[rp,]
  K <- qr.Q(Q)[1:n,]
  } else {  ## need SVD step
    R <- R[,rp] ## unpivoted R
    s1 <- svd(R)
    d.inv1 <- rep(0,q)
    good1 <- s1$d >= max(s1$d)*.Machine$double.eps^.5
    d.inv1[good1] <- 1/s1$d[good1]
    P <- t(d.inv1[good1]*t(s1$v[,good1]))
    K <- qr.qy(Q,rbind(s1$u,matrix(0,n,q))[,good1])[1:n,]         
  }
  Vb <- tcrossprod(P) * sig2 
  ## P%*%t(P)*sig2 # Bayesian posterior covariance matrix for the parameters 
  Ve <- crossprod(K%*%t(P)) *sig2
  #PKt%*%t(PKt)*sig2 # covariance matrix of the parameter estimators 
  ## Delta method to get covariance matrix for the reparametrized parameters...
  df.p <- rep(1,q)
  df.p[object$iv] <- object$beta.t[object$iv]
  Vb.t <- t(df.p*t(df.p*Vb))
  Ve.t <- t(df.p*t(df.p*Ve))
  
  ## calculating edf and trA...
  KtILQ1R <- crossprod(object$L*object$I.plus*K,object$wX1) ## t(object$L*object$I.plus*K)%*%object$wX1
  F <- P%*%(KtILQ1R)
  edf <- diag(F) ## effective degrees of freedom
  edf1 <- 2*edf - rowSums(t(F)*F) ## alternative
  trA <- sum(edf)
  list (nulldev=nulldev, df.null=nulldf,Vb=Vb,Vb.t=Vb.t,Ve=Ve,Ve.t=Ve.t,
        sig2=sig2,edf=edf,edf1=edf1,trA=trA)
} ## end of scam.fit.post



### the following three functions are for use in place of exp(beta)
### notExp() is similar to that in R package mgcv() of Simon N Wood
### in positivity ensuring beta parameters re-parameterization.... they have `better' 
### over/underflow characteristics, but is still continuous to second
### derivative. 
### DnotExp() calculates the first derivative
### D2notExp() gets the second derivative 

notExp <- function(x){
  f <- x
  ind <- x > 1
  f[ind] <- exp(1)*(x[ind]^2+1)/2
  ind <- (x <= 1)&(x > -1)
  f[ind] <- exp(x[ind])
  ind <- (x <= -1)
  f[ind] <-  exp(1)*(x[ind]^2+1)/2; f[ind]<-1/f[ind]
  f
}

DnotExp <- function(x) {
  ## first derivative of notExp()...
  f <- x
  ind <- x > 1
  f[ind] <- exp(1)*x[ind]
  ind <- (x <= 1)&(x > -1)
  f[ind] <- exp(x[ind])
  ind <- (x <= -1)
  f[ind] <-  -4*x[ind]/exp(1)/(x[ind]^2+1)^2
  f
}

D2notExp <- function(x) {
  ## second derivative of notExp()...
  f <- x
  ind <- x > 1
  f[ind] <- exp(1)
  ind <- (x <= 1)&(x > -1)
  f[ind] <- exp(x[ind])
  ind <- (x <= -1)
  f[ind] <-  (12*x[ind]^2-4)/exp(1)/(x[ind]^2+1)^3
  f
}


D3notExp <- function(x) {
  ## third derivative of notExp()...
  f <- x
  ind <- x > 1
  f[ind] <- 0
  ind <- (x <= 1)&(x > -1)
  f[ind] <- exp(x[ind])
  ind <- (x <= -1)
  f[ind] <-  48*x[ind]*(1-x[ind]^2)/exp(1)/(x[ind]^2+1)^4
  f
}





###############################################################
## loading functions, copied from mgcv() package of Simon Wood
#################################################################

# print.scam.version <- function()
# { library(help=scam)$info[[1]] -> version
#   version <- version[pmatch("Version",version)]
#   um <- strsplit(version," ")[[1]]
#   version <- um[nchar(um)>0][2]
#   hello <- paste("This is scam ",version,".",sep="")
#   packageStartupMessage(hello)
# }


# .onAttach <- function(...) { 
#   print.scam.version()
# }

##.onUnload <- function(libpath) library.dynam.unload("scam", libpath)









#########################################################
# Function to return gcv/ubre ...                      ##
#########################################################

gcv.ubre <- function(rho,G,gamma,env) 
{  ## function to get GCV.UBRE value for optim()...
  if (length(rho)!= length(G$off)) stop (paste("length of rho and n.terms has to be the same"))
  sp <- exp(rho)
  b <- scam::scam.fit(G=G, sp=sp, env=env) 
  if (G$scale.known) #  value of Mallow's Cp/UBRE/AIC ....
  {  n <- nrow(G$X)
  gcv.ubre <- b$dev/n - G$sig2 +2*gamma*b$trA*G$sig2/n
  }  else   # value of GCV ...
    gcv.ubre <- b$gcv
  return(gcv.ubre)
}

#########################################################
## function to get the gradient of the gcv/ubre.....   ##
#########################################################

gcv.ubre.derivative <- function(rho,G, gamma,env, check.analytical=FALSE, del)  
{  ## function to return derivative of GCV or UBRE for optim...
  scam::gcv.ubre_grad(rho, G, gamma,env,check.analytical, del)$gcv.ubre.rho
}


#############################################################################
## for nlm() function to get the gcv/ubre and gradient of the gcv/ubre.....##
#############################################################################

dgcv.ubre.nlm <- function(rho,G, gamma,env, check.analytical=FALSE, del) 
{  ## GCV UBRE objective function for nlm
  gg <- scam::gcv.ubre_grad(rho, G, gamma,env,check.analytical, del) 
  attr(gg$gcv.ubre,"gradient") <- gg$gcv.ubre.rho
  gg$gcv.ubre
}



#######################################################
#### estimate.scam()....                             ##
#######################################################


estimate.scam <- function(G,optimizer,optim.method,rho, gamma,env,
                          check.analytical, del, devtol, steptol)
{  ## function to select smoothing parameter...
  if (!(optimizer %in% c("bfgs", "nlm", "optim","nlm.fd")) )
    stop("unknown outer optimization method")
  if (optimizer == "bfgs") ## minimize GCV/UBRE by BFGS...
  {  b <- bfgs_gcv.ubre(scam::gcv.ubre_grad,rho=rho, G=G,gamma=gamma,env=env,
                        check.analytical=check.analytical, del=del) 
  sp <- exp(b$rho)
  object <- b$object
  object$gcv.ubre <- b$gcv.ubre
  object$dgcv.ubre <- b$dgcv.ubre
  object$termcode <- b$termcode
  object$check.grad <- b$check.grad
  object$dgcv.ubre.check <- b$dgcv.ubre.check
  object$conv.bfgs <- b$conv.bfgs
  object$iterations <- b$iterations
  }
  else if (optimizer=="optim")  ## gr=gcv.ubre.derivative
  {  if (!(optim.method[1] %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN")))
  {   warning("unknown optim() method, `L-BFGS-B' were used")
    optim.method[1] <- "L-BFGS-B"
  }
    if (is.na(optim.method[2])) 
    {   warning("the second parameter of optim.method argument is not supplied, 
                finite-difference approximation of the gradient were used")
      grr <- NULL
    }  else if (!(optim.method[2] %in% c("fd","grad")))
    {   warning("only `fd' and `grad' options are possible, finite-difference 
                approximation of the gradient were used")
      grr <- NULL
    }  else if (optim.method[2] == "grad")
      grr <- gcv.ubre.derivative
    else 
      grr <- NULL
    b <- stats::optim(par=rho,fn=gcv.ubre, gr=grr, method=optim.method[1],G=G, gamma=gamma,env=env) 
    sp <- exp(b$par)
    gcv.ubre <- b$value
    dgcv.ubre <- NULL
    iterations <- b$counts
    termcode <- b$convergence
    if (termcode == 0)
      conv <- "Successful completion"
    else if (termcode == 1)  
      conv <- "The iteration limit `maxit' had been reached"
    else if (termcode == 10)  
      conv <- "Degeneracy of the Nelder-Mead simplex"
    else if (termcode == 51)  
      conv <- "A warning from the `L-BFGS-B' method; see help for `optim' for further details"
    else if (termcode == 52)  
      conv <- "An error from the `L-BFGS-B' method; see help for `optim' for further details"
  }
  else if (optimizer=="nlm.fd") ## nlm() with finite difference derivatives...
  {  b <- stats::nlm(f=gcv.ubre, p=rho,iterlim=100, G=G, gamma=gamma,env=env) 
  }
  else if (optimizer=="nlm")  ## nlm() with analytical derivatives...
  { b <- stats::nlm(f=dgcv.ubre.nlm, p=rho,iterlim=100,G=G,gamma=1,env=env,
                    check.analytical=check.analytical, del=del) 
  }
  if (optimizer== "nlm.fd" || optimizer== "nlm") 
  {   sp <- exp(b$estimate)
  gcv.ubre <- b$minimum
  dgcv.ubre <- b$gradient
  iterations <- b$iterations
  termcode <- b$code
  if (termcode == 1)
    conv <- "Relative gradient is close to zero, current iterate is probably solution"
  else if (termcode == 2)  
    conv <- "Successive iterates within tolerance, current iterate is probably solution"
  else if (termcode == 3)  
    conv <- "Last global step failed to locate a point lower than `estimate'. Either 
  `estimate' is an approximate local minimum of the function or 
  `steptol' is too small"
  else if (termcode == 4)  
    conv <- "Iteration limit exceeded"
  else if (termcode == 5)  
    conv <- "Maximum step size `stepmax' exceeded five consecutive 
  times. Either the function is unbounded below, becomes asymptotic 
  to a finite value from above in some direction or stepmax is too small"   
  }
  ## fit the model using the optimal sp from "optim" or "nlm"...
  if (optimizer== "nlm.fd" || optimizer== "nlm" || optimizer== "optim")
  {  object <- scam::scam.fit(G=G, sp=sp,env=env, devtol=devtol, steptol=steptol) 
  object$gcv.ubre <- gcv.ubre
  object$dgcv.ubre <- dgcv.ubre 
  object$termcode <- termcode
  object$conv <- conv
  object$iterations <- iterations 
  }
  if (optimizer=="optim")
  {  object$optim.method <- rep(NA,2)
  object$optim.method[1] <- optim.method[1]
  if (!is.null(grr))
    object$optim.method[2] <- "grad"
  }
  object$sp <- sp
  object$q.f <- G$q.f
  object$p.ident <- G$p.ident
  object$S <- G$S
  object
  }





## BFGS for GCV/UBRE minimization 
## with modifications on convergence control, scaling, etc...

#####################################################
## BFGS for gcv/ubre miminization..                ##
#####################################################

bfgs_gcv.ubre <- function(fn=scam::gcv.ubre_grad, rho, ini.fd=TRUE, G, gamma=1, env,
                          n.pen=length(rho), typx=rep(1,n.pen), typf=1, steptol= 1e-7, 
                          gradtol = 6.0554*1e-06, maxNstep = 5, maxHalf = 30, 
                          check.analytical, del)  
{  ## fn - GCV/UBRE Function which returs the GCV/UBRE value and its derivative wrt log(sp)
  ## rho - log of the initial values of the smoothing parameters
  ## ini.fd - if TRUE, a finite difference to the Hessian is used to find the initial inverse Hessian
  ## typx - vector whose component is a positive scalar specifying the typical magnitude of sp
  ## typf - a positive scalar estimating the magnitude of the gcv near the minimum
  ## gradtol - a scalar giving a tolerance at which the gradient is considered
  ##            to be close enougth to 0 to terminate the algorithm
  ## steptol - a positive scalar giving the tolerance at which the scaled distance between
  ##          two successive iterates is considered close enough to zero to terminate the algorithm 
  ## maxNstep - a positive scalar which gives the maximum allowable step length
  ## maxHalf - a positive scalar which gives the maximum number of step halving in "backtracking"
  
  Sp <- 1/typx     # diagonal of the scaling matrix 
  ## storage for solution track
  rho1 <- rho
  old.rho <- rho
  not.exp <- G$not.exp
  b <- fn(rho,G,gamma=1, env, check.analytical=FALSE, del) 
  old.score <- score <- b$gcv.ubre
  score.plot <- rep(0,200) ## for plotting the gcv
  score.plot[1] <- score
  grad <- b$dgcv.ubre
  ## The initial inverse Hessian ...
  if (ini.fd) 
  {  B <- matrix(0,n.pen,n.pen)
  for (j in 1:n.pen) 
  {  rho2 <- rho; rho2[j] <- rho[j] + 1e-6
  b2 <- fn(rho2,G,gamma=1,env,check.analytical=FALSE, del) 
  B[,j] <- (b2$dgcv.ubre - grad)/1e-6   
  }
  B <- B + t(B)
  eh <- eigen(B)
  ind <- eh$values < 0
  eh$values[ind] <- -eh$values[ind] 
  B <- eh$vectors%*%(t(eh$vectors)/eh$values)
  } 
  else  
    B <- diag(n.pen)*100
  score.scale <- b$scale.est + score  ## ??? for UBRE
  unconv.ind <- abs(grad) > score.scale*gradtol # checking the gradient is within the tolerance
  if (!sum(unconv.ind))  ## if at least one is false
    unconv.ind <- unconv.ind | TRUE
  consecmax <- 0
  ## Quasi-Newton algorithm to minimize GCV...
  for (i in 1:200) 
  {  ## compute a BFGS search direction ...
    Nstep <- 0*grad   ## initialize the quasi-Newton step
    Nstep[unconv.ind] <- -drop(B[unconv.ind, unconv.ind]%*%grad[unconv.ind])
    Dp <- Sp*Nstep
    Newtlen <- (sum(Dp^2))^.5  ## Euclidean norm to get the Newton step
    if (Newtlen > maxNstep)  ## reduce if ms is greater than the max allowable
    {  Nstep <- maxNstep*Nstep/Newtlen
    Newtlen <- maxNstep
    }
    maxtaken <- FALSE
    retcode <- 2
    initslope <- sum(grad * Nstep) ## initial slope
    rellength <- max(abs(Nstep)/max(abs(rho),1/Sp)) ## relative length of rho for the stopping criteria
    minalpha <- steptol/rellength
    c1 <- 1e-4   ## constant for the sufficient decrease condition 
    alpha <- 1   ## initialize step length
    ii <- 0 ## initialize the number of "step halving"
    step <- alpha*Nstep 
    ## step length selection ...
    curv.condition <- TRUE
    repeat 
    {  rho1 <- rho + alpha*Nstep 
    b <- fn(rho=rho1,G,gamma=1,env,check.analytical=FALSE, del) 
    score1 <- b$gcv.ubre
    if (score1 <= score+c1*alpha*initslope) 
    {   grad1 <- b$dgcv.ubre
    newslope <- sum(grad1 * Nstep)
    curv.condition <- TRUE
    if (newslope < 0.9*initslope) # the curvature condition is not satisfied
    {   if (alpha == 1 && Newtlen < maxNstep)
    {   maxalpha <- maxNstep/Newtlen
    repeat 
    {  old.alpha <- alpha
    old.score1 <- score1
    alpha <- min(2*alpha, maxalpha)
    rho1 <- rho + alpha*Nstep
    b <- fn(rho=rho1,G,gamma=1, env,
            check.analytical=FALSE, del) 
    score1 <- b$gcv.ubre
    if (score1 <= score+c1*alpha*initslope)
    {   grad1 <- b$dgcv.ubre
    newslope <- sum(grad1*Nstep)
    }
    if (score1 > score+c1*alpha*initslope)
      break
    if (newslope >= 0.9*initslope)
      break
    if (alpha >= maxalpha)
      break
    } 
    }
      if ((alpha < 1) || (alpha>1 && (score1>score+c1*alpha*initslope)))
      {   alpha.lo <- min(alpha, old.alpha)
      alpha.diff <- abs(old.alpha - alpha)
      if (alpha < old.alpha)
      {  sc.lo <- score1
      sc.hi <- old.score1
      }    
      else
      {  sc.lo <- old.score1
      sc.hi <- score1 
      }
      repeat
      {  alpha.incr <- -newslope*alpha.diff^2/(2*(sc.hi-(sc.lo+newslope*alpha.diff)))
      if (alpha.incr < 0.2*alpha.diff) 
        alpha.incr <- 0.2*alpha.diff
      alpha <- alpha.lo+alpha.incr
      rho1 <- rho + alpha*Nstep
      b <- fn(rho=rho1,G,gamma=1, env,
              check.analytical=FALSE, del) 
      score1 <- b$gcv.ubre
      if (score1 > score+c1*alpha*initslope)
      {  alpha.diff <- alpha.incr
      sc.hi <- score1
      }
      else
      {   grad1 <- b$dgcv.ubre
      newslope <- sum(grad1*Nstep)
      if (newslope < 0.9*initslope)
      {  alpha.lo <- alpha
      alpha.diff <- alpha.diff-alpha.incr
      sc.lo <- score1
      }
      }
      if (newslope >= 0.9*initslope)
        break
      if (alpha.diff < minalpha)
        break 
      }
      if (newslope < 0.9*initslope)   ## couldn't satisfy curvature condition
      {   curv.condition <- FALSE
      score1 <- sc.lo
      rho1 <- rho + alpha.lo*Nstep
      b <- fn(rho=rho1,G,gamma=1,env,check.analytical=FALSE, del) 
      } 
      } ## end of "if ((alpha < 1) || (alpha>1 && ..."
    }  ## end of "if (newslope < 0.9*initslope) ..."
    retcode <- 0
    if (newslope < 0.9*initslope) ## couldn't satisfy curvature condition
      curv.condition <- FALSE
    if (alpha*Newtlen > 0.99*maxNstep)  
      maxtaken <- TRUE
    } ## end of "if (score1 <= ...) ..."
    else if (alpha < minalpha) ## no satisfactory rho+ can be found suff-ly distinct from previous rho
    {   retcode <- 1
    rho1 <- rho
    b <- fn(rho=rho1,G,gamma=1,env,check.analytical=FALSE,del)  
    }
    else   ## backtracking to satisfy the sufficient decrease condition...
    {   ii <- ii+1
    if (alpha == 1) ## first backtrack, quadratic fit
    {   alpha.temp <- -initslope/(score1-score-initslope)/2
    }
    else {  ## all subsequent backtracts, cubic fit 
      A1 <- matrix(0,2,2)
      bb1 <-rep(0,2)
      ab <- rep(0,2)
      A1[1,1] <- 1/alpha^2
      A1[1,2] <- -1/old.alpha^2
      A1[2,1] <- -old.alpha/alpha^2
      A1[2,2] <- alpha/old.alpha^2
      bb1[1] <- score1-score-alpha*initslope
      bb1[2] <- old.score1 -score-old.alpha*initslope
      ab <- 1/(alpha-old.alpha)*A1%*%bb1
      disc <- ab[2]^2-3*ab[1]*initslope
      if (ab[1] == 0) ## cubic is a quadratic
        alpha.temp <- -initslope/ab[2]/2
      else     ## legitimate cubic
        alpha.temp <- (-ab[2]+disc^0.5)/(3*ab[1])
      if (alpha.temp > 0.5*alpha)
        alpha.temp <- 0.5*alpha   
    }
    old.alpha <- alpha
    old.score1 <- score1
    if (alpha.temp <= 0.1*alpha) alpha <- 0.1*alpha
    else alpha <- alpha.temp 
    }
    if (ii == maxHalf)
      break
    if (retcode < 2)
      break 
    } ## end of REPEAT for the step length selection
    ## rho1 is now new point. 
    step <- alpha*Nstep
    old.score <-score
    old.rho <- rho
    rho <- rho1
    old.grad <- grad
    score <- score1
    grad <- b$dgcv.ubre 
    score.plot[i+1] <- score
    ## update B...
    yg <- grad - old.grad
    rr <- 1/sum(yg * step)
    skipupdate <- TRUE 
    ## skip update if `step' is sufficiently close to B%*%yg ...
    for (i in 1:n.pen)
    {  closeness <- step[i]-B[i,]%*%yg  
    if (abs(closeness) >= gradtol*max(abs(grad[i]),abs(old.grad[i])))   
      skipupdate <- FALSE
    }
    ## skip update if curvature condition is not satisfied...
    if (!curv.condition)
      skipupdate <- TRUE
    if (!skipupdate) 
    {  B <- B - rr * step %*% crossprod(yg,B) # (t(yg)%*%B) 
    B <- B - rr*tcrossprod((B %*% yg),step) + rr *tcrossprod(step)  # B - rr*(B %*% yg) %*% t(step) + rr * step %*% t(step)
    }
    
    ## check the termination condition ...
    termcode <- 0
    if (retcode ==1) 
      termcode <- 3
    else if (max(abs(grad)*max(abs(rho),1/Sp)/max(abs(score),typf))<= gradtol)
      termcode <- 1
    else if (max(abs(rho-old.rho)/max(abs(rho),1/Sp))<= steptol)
      termcode <- 2
    else if (i==200) 
      termcode <- 4 
    else if (maxtaken) ## step of length maxNstep was taken
    {  consecmax <- consecmax +1
    if (consecmax ==5)
      termcode <- 5 # limit of 5 maxNsteps was reached
    }
    else consecmax <- 0 
    ##---------------------
    if (termcode > 0)
      break
    else  ## if not converged...
    {   converged <- TRUE
    score.scale <- b$scale.est + score
    unconv.ind <- abs(grad) > score.scale * gradtol
    if (sum(unconv.ind))
      converged <- FALSE
    if (abs(old.score - score) > score.scale*gradtol) 
    {  if (converged)
      unconv.ind <- unconv.ind | TRUE
    converged <- FALSE
    }
    } # end of ELSE
  } ## end of the Quasi-Newton algorithm 
  
  ## printing why the algorithm terminated... 
  if (termcode == 1)
    ct <- "Full convergence"
  else if (termcode == 2)
    ct <- "Successive iterates within tolerance, current iterate is probably solution"
  else if (termcode == 3)
    ct <- "Last step failed to locate a lower point than old.rho"
  else if (termcode == 4)
    ct <- "Iteration limit reached"
  else if (termcode ==5)
    ct <- "Five conseqcutive steps of length maxNstep have been taken" 
  list (gcv.ubre=score, rho=rho, dgcv.ubre=grad, iterations=i, B=B, conv.bfgs = ct, object=b$object, score.plot=score.plot[1:(i+1)], termcode = termcode, check.grad= b$check.grad,
        dgcv.ubre.check = b$dgcv.ubre.check) 
} ## end bfgs_gcv.ubre


