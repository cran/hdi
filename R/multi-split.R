multi.split <- function(x, y, B = 50, fraction = 0.5,
                        model.selector = lasso.cv,
                        classical.fit = lm.pval,
                        gamma = seq(0.05, 0.99, by = 0.01),
                        args.model.selector = NULL,
                        args.classical.fit = NULL,
                        return.nonaggr = FALSE,
                        return.selmodels = FALSE,
                        trace = FALSE)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  2 Apr 2013, 11:52

  n <- nrow(x)
  p <- ncol(x)
  
  ## Matrix of bootstrap p-values:
  ## rows = sample-splits
  ## cols = predictors
  pvals <- matrix(1, nrow = B, ncol = p)
  colnames(pvals) <- colnames(x)

  if(return.selmodels){
    sel.model.all <- matrix(FALSE, nrow = B, ncol = p)
    colnames(sel.model.all) <- colnames(x)
  }else{ ## safe memory space in case no output is wanted regarding sel. models
    sel.model.all <- NULL
  }
  
  n.left <- floor(n * fraction)
  n.right <- n - n.left
      
  for(b in 1:B){
    try.again <- TRUE
    repeat.count <- 0
    while(try.again){
      if(trace)
        cat("...Split", b, "\n")
        
      ## Perform sample-splitting; sample *without* replacement
      split <- sample(1:n, n.left, replace = FALSE)
    
      ## x.left is used to perform model selection
      x.left <- x[split,]
      y.left <- y[split]
      
      ## x.right is used to calculate p-values
      x.right <- x[-split,]
      y.right <- y[-split]
      
      sel.model <- do.call(model.selector,
                           args = c(list(x=x.left, y=y.left),
                               args.model.selector))

      p.sel <- length(sel.model)

      ## Classical situation:
      ## A model with intercept is used, hence p.sel + 1 < nrow(x.right),
      ## otherwise, p-values can not be calculated
      if(p.sel > 0 & p.sel < nrow(x.right) - 1){ 
        sel.pval <- do.call(classical.fit,
                            args = c(list(x = x.right[,sel.model],
                                y = y.right), args.classical.fit))
        ## Bonferroni on small model
        pvals[b, sel.model] <- pmin(sel.pval * p.sel, 1)

        if(return.selmodels)
          sel.model.all[b, sel.model] <- TRUE
        
        try.again <- FALSE ## break the loop, continue with next sample-split
      }
      ## Empty model selected:
      ## Do nothing in that case. Matrix already filled with 1's.
      ## Print out information for the sake of completeness
      if(p.sel == 0){
        if(trace)
          cat("......Empty model selected. That's ok...\n")

        try.again <- FALSE ## break the loop, continue with next sample-split
      }
      ## Too large model selected for classical fitter
      if(p.sel >= n.right - 1){ ## p.sel + 1 < n.right for p-val calculation
        try.again <- TRUE ## re-do sample splitting
        repeat.count <- repeat.count + 1
        warning("Too large model selected in a sample-split")
      }
      if(repeat.count > 20){ ## too prevent never-ending loops
        stop("More than 20 sample splits resulted in too large models...giving up")
        try.again <- FALSE
      }
    } ## end while(try.again)
  }
  
  ## For loop is not very innovative, but it does it's job...
  pvals.current <- which.gamma <- numeric(p)
  for(j in 1:p){ ## loop through all predictors
    quant.gamma <- quantile(pvals[,j], gamma) / gamma

    if(length(gamma) > 1)
      penalty <- (1 - log(min(gamma)))
    else
      penalty <- 1

    pvals.pre <- min(quant.gamma) * penalty
             
    pvals.current[j] <- pmin(pvals.pre, 1)
    
    which.gamma[j] <- which.min(quant.gamma)
  }
  
  names(pvals.current) <- names(which.gamma) <- colnames(x)

  if(!return.nonaggr) ## Overwrite pvals with NULL if no output is wanted
    pvals <- NULL
      
  out <- list(pval          = pvals.current,
              gamma.min     = gamma[which.gamma],
              pvals.nonaggr = pvals,
              sel.models    = sel.model.all,
              method        = "multi.split",
              call          = match.call())
  
  class(out) <- "hdi"
  
  return(out)
}





