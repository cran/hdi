lasso.proj <- function(x, y, family = "gaussian",
                       standardize = TRUE,
                       multiplecorr.method = "holm",
                       N = 10000,
                       parallel = FALSE, ncores = 4,
                       sigma = NULL, ## sigma estimate provided by the user
                       Z = NULL)     ## Z or Thetahat provided by the user
{
  ## Purpose:
  ## An implementation of the LDPE method http://arxiv.org/abs/1110.2563
  ## which is identical to
  ## http://arxiv.org/abs/1303.0518 
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Return values:
  ## pval: p-values for every parameter (individual tests)
  ## pval.corr:  multiple testing corrected p-values for every parameter
  ## betahat:    initial estimate by the scaled lasso of \beta^0
  ## bhat:       de-sparsified \beta^0 estimate used for p-value calculation
  ## sigmahat:   \sigma estimate coming from the scaled lasso
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 18 Oct 2013 (initial version),
  ## in part based on an implementation of the ridge projection method
  ## ridge-proj.R by P. Buehlmann + adaptations by L. Meier.

  ####################
  ## Get data ready ##
  ####################
  
  n <- nrow(x)
  p <- ncol(x)

  if(standardize)
    sds <- apply(x, 2, sd)
  else
    sds <- rep(1, p)

  ## *center* (scale) the columns
  x <- scale(x, center = TRUE, scale = standardize)

  dataset <- switch(family,
                    "gaussian" = {
                      list(x = x, y = y)
                    },
                    {
                      switch.family(x = x, y = y,
                                    family = family)
                    })

  ## center the columns the response to get rid of the intercept
  x <- scale(dataset$x, center = TRUE, scale = FALSE)
  y <- scale(dataset$y, scale = FALSE)
  y <- as.numeric(y)

  ## Warning: should we allow user to still specify their
  ## own Z here?
  ## Z <- NULL##will have to recalculate Z
  
  ## force sigmahat to 1 when doing glm!
  if(family == "binomial")
    sigma <- 1

  ######################################
  ## Calculate Z using nodewise lasso ##
  ######################################
  
  if(is.null(Z)){
    nodewiselasso.out <- score.nodewiselasso(x = x,
                                             wantTheta = FALSE,
                                             parallel = parallel,
                                             ncores = ncores)
    Z <- nodewiselasso.out$out
  }else{
    ## Check if normalization is fulfilled
    if(!all.equal(rep(1, p), colSums(Z * x) / n, tolerance = 10^-8)){
      print("Z not properly normalized...")
      print("...we assume all(colSums(Z * x) / n == rep(1, p))")
      print("...rescaling Z ourselves.")
      Z <- score.rescale(Z = Z, x = x)
    }
  }
  

  ###################################
  ## Projection estimator and bias ##
  ###################################
  
  bproj <- t(Z) %*% y / n
  
  ## Bias estimate based on lasso estimate
  scaledlassofit <- scalreg(X = x, y = y, lam0 = "univ")
  betalasso      <- scaledlassofit$coefficients

  ## Get estimated standard deviation
  if(is.null(sigma))
    sigmahat <- scaledlassofit$hsigma
  else
    sigmahat <- sigma

  ## Subtract bias
  bias <- numeric(p)
  for(j in 1:p){ ## replace loop?
    bias[j] <- (t(Z[,j]) %*% x[,-j]) %*% betalasso[-j] / n
  }
  
  bproj <- bproj - bias

  #########################
  ## p-Value calculation ##
  #########################

  ## Determine normalization factor
  scaleb        <- n / (sigmahat * sqrt(colSums(Z^2))) ##sqrt(diag(crossprod(Z)))
  bprojrescaled <- bproj * scaleb
  
  ## Calculate p-value
  pval <- 2 * pnorm(abs(bprojrescaled), lower.tail = FALSE)

  cov2 <- crossprod(Z)

  #################################
  ## Multiple testing correction ##
  #################################
  
  if(multiplecorr.method == "WY"){
    ## Westfall-Young like procedure as in ridge projection method,
    ## P.Buhlmann & L.Meier
    ## method related to the Westfall - Young procedure
    ## constants left out since we'll rescale anyway
    ## otherwise cov2 <- crossprod(Z)/n
    pcorr <- p.adjust.wy(cov = cov2, pval = pval, N = N)
  }else{
    if(multiplecorr.method %in% p.adjust.methods){
      pcorr <- p.adjust(pval,method = multiplecorr.method)
    }else{
      stop("Unknown multiple correction method specified")
    }
  }## end multiple testing correction
  
  ## Also return the confidence intervals
  se   <- 1 / scaleb

  
  ##############################################
  ## Function to calculate p-value for groups ##
  ##############################################
  
  pre <- preprocess.group.testing(N = N, cov = cov2, alt = FALSE)
  
  group.testing.function <- function(group, alt = TRUE){
    calculate.pvalue.for.group(brescaled  = bprojrescaled,
                               group      = group,
                               individual = pval,
                               ##correct    = TRUE,
                               alt        = alt,
                               zz2        = pre)
  }

  cluster.group.testing.function <-
    get.clusterGroupTest.function(group.testing.function =
                                  group.testing.function, x = x)
  
  ############################
  ## Return all information ##
  ############################
  
  out <- list(pval        = as.vector(pval),
              pval.corr   = pcorr,
              groupTest   = group.testing.function,
              clusterGroupTest = cluster.group.testing.function,
              sigmahat    = sigmahat,
              standardize = standardize,
              sds         = sds,
              bhat        = bproj / sds,
              se          = se / sds,
              betahat     = betalasso / sds,
              family      = family,
              method      = "lasso.proj",
              call        = match.call())

  names(out$pval) <- names(out$pval.corr) <- names(out$bhat) <-
    names(out$sds) <- names(out$se) <- names(out$betahat) <-
      colnames(x)

  class(out) <- "hdi"
  return(out)
}

