multi.split <- function(x, y, B = 100, fraction = 0.5,
                        ci = TRUE, ci.level = 0.95,
                        model.selector = lasso.cv,
                        classical.fit = lm.pval,
                        classical.ci  = lm.ci,
                        parallel = FALSE, ncores = 4,
                        gamma = seq(ceiling(0.05*B)/B,1-1/B,by=1/B),                        
                        args.model.selector = NULL,
                        args.classical.fit = NULL,
                        args.classical.ci = NULL,
                        return.nonaggr = FALSE,
                        return.selmodels = FALSE,
                        trace = FALSE)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  2 Apr 2013, 11:52
  ## Updated with confidence interval calculation, Ruben Dezeure (5 Feb 2014)

  n <- nrow(x)
  p <- ncol(x)
  
  n.left  <- floor(n * fraction)
  n.right <- n - n.left

  ## Helper function for one sample split
  ## all variables are taken from the function environment
  
  oneSplit <- function(...){
    pvals.v         <- rep(1, p)
    sel.model.all.v <- logical(p) ## FALSE by default
    lci.v           <- rep(-Inf, p)
    uci.v           <- rep(Inf, p)
    centers.v       <- rep(NA, p)
    ses.v           <- rep(Inf, p)
    df.res          <- NA

    try.again    <- TRUE
    repeat.count <- 0
    
    while(try.again){
      ## Perform sample-splitting; sample *without* replacement
      split <- sample(1:n, n.left, replace = FALSE)
    
      ## x.left is used to perform model selection
      x.left <- x[split,]
      y.left <- y[split]
      
      ## x.right is used to calculate p-values
      x.right <- x[-split,]
      y.right <- y[-split]
      
      sel.model <- do.call(model.selector,
                           args = c(list(x = x.left, y = y.left),
                             args.model.selector))

      p.sel <- length(sel.model)

      ## Classical situation:
      ## A model with intercept is used, hence p.sel + 1 < nrow(x.right),
      ## otherwise, p-values can *not* be calculated
      if(p.sel > 0 & p.sel < nrow(x.right) - 1){
        sel.pval <- do.call(classical.fit,
                            args = c(list(x = x.right[,sel.model],
                              y = y.right), args.classical.fit))
        ## Bonferroni on small model
        pvals.v[sel.model] <- pmin(sel.pval * p.sel, 1) ## new
        
        if(ci){ ## Calculations of confidence intervals
          if(identical(classical.fit, lm.pval)){
            ## Calculate ci's and save all necessary information
            tmp.fit.lm <- lm(y.right ~ x.right[,sel.model],
                             args.classical.fit)
            ## Based on code from getAnywhere(confint.lm)
            a   <- (1 - ci.level) / 2
            a   <- c(a, 1 - a)
            fac <- qt(a, tmp.fit.lm$df.residual)
            sel.ses     <- sqrt(diag(vcov(tmp.fit.lm)))[-1] ## without intercept
            sel.centers <- coef(tmp.fit.lm)[-1]
            sel.ci      <- sel.centers + sel.ses %o% fac
            centers.v[sel.model] <- sel.centers
            lci.v[sel.model]     <- sel.ci[,1] ## new
            uci.v[sel.model]     <- sel.ci[,2] ## new
            ses.v[sel.model]     <- sel.ses    ## change?
            df.res               <- tmp.fit.lm$df.residual ## change?          
          }else{
            ## do the primitive ci interval aggregation
            sel.ci <- do.call(classical.ci,
                              args =  c(list(x = x.right[, sel.model],
                                y = y.right, level = 1 - (1 - ci.level) / 2),
                                args.classical.ci))
            lci.v[sel.model] <- sel.ci[, 1] ## new
            uci.v[sel.model] <- sel.ci[, 2] ## new
          }
        } ## end if(ci)

        if(return.selmodels)
          sel.model.all.v[sel.model] <- TRUE
        
        try.again <- FALSE ## break the loop, continue with next sample-split
      }
      ## Empty model selected:
      ## Do nothing in that case. Matrix already filled with 0's.
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
      if(repeat.count > 20){ ## to prevent never-ending loops
        stop("More than 20 sample splits resulted in too large models...giving up")
        try.again <- FALSE
      }
    } ## end while(try.again)
    return(list(pvals         = pvals.v,
                sel.model.all = sel.model.all.v,
                centers       = centers.v,
                ses           = ses.v,
                df.res        = df.res,
                lci           = lci.v,
                uci           = uci.v))
  }

  ######################
  ## Sample-splitting ##
  ######################
  
  if(parallel){
    if(trace)
      cat("...starting parallelization of sample-splits\n")
 
    split.out <- mclapply(1:B, oneSplit, mc.cores = ncores)
  }
  else{
    split.out <- list(); length(split.out) <- B
    for(b in 1:B){
      if(trace)
        cat("...Split", b, "\n")
      
      split.out[[b]] <- oneSplit()
    }
  }

  #####################
  ## Deparse objects ##
  #####################
  
  myDeparse <- function(name){
    matrix(unlist(lapply(split.out, "[[", name)), nrow = B, byrow = TRUE)
  }

  ## Matrix of bootstrap p-values (etc.)
  ## rows = sample-splits
  ## cols = predictors

  pvals <- myDeparse("pvals"); colnames(pvals) <- colnames(x)
  
  if(return.selmodels){
    sel.model.all <- myDeparse("sel.model.all")
    colnames(sel.model.all) <- colnames(x)
  }else{ ## safe memory space in case no output is wanted regarding sel. models
    sel.model.all <- NA
  }
  
  lci           <- myDeparse("lci")
  uci           <- myDeparse("uci")
  centers       <- myDeparse("centers")
  ses           <- myDeparse("ses")
  df.res        <- unlist(lapply(split.out, "[[", "df.res"))

  ##############################
  ## Calculate final p-values ##
  ##############################
  
  ## For loop is not very innovative, but it does it's job...
  pvals.current <- which.gamma <- numeric(p)
  if(!(0.05 %in% gamma))
    warning("0.05 is not in the gamma range due to the choice of B, the results might be incorrect. Pick a B such that some integer multiple of 1/B equals 0.05 to avoid this.")
  
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

  if(ci && identical(classical.fit,lm.pval)){
    vars <- ncol(lci)
    ##for every separate variable, aggregate the ci
    if(any(is.na(sel.model.all)))
      {
        s0 <- NA
      }else{
        s0 <- apply(sel.model.all,1,sum)
      }
    ##calculate single-testing confidence intervals
    new.ci <- mapply(aggregate.ci,
                     lci = split(lci, rep(1:vars, each = B)),
                     rci = split(uci, rep(1:vars, each = B)),
                     centers = split(centers,rep(1:vars, each = B)),
                     ses = split(ses, rep(1:ncol(ses), each = B)),
                     df.res = list(df.res = df.res),
                     gamma.min = min(gamma),
                     ##multi.corr = TRUE,##temporarily trying multiple testing corrected ci
                     multi.corr = FALSE,##single testing confidence intervals
                     verbose = FALSE,
                     s0=list(s0=s0),
                     var=1:vars)#for debug information to have the var of this aggregate call
    lci.current <- t(new.ci)[,1]
    uci.current <- t(new.ci)[,2]
  }else{
    lci.current <- apply(lci, 2, median)
    uci.current <- apply(uci, 2, median)
  }

  if(!return.nonaggr) ## Overwrite pvals with NULL if no output is wanted
    pvals <- NA

  names(lci.current) <- names(uci.current) <- names(pvals.current)
  
  out <- list(pval          = NA,
              pval.corr     = pvals.current,
              lci           = if (ci) lci.current else NA,
              uci           = if (ci) uci.current else NA,
              gamma.min     = gamma[which.gamma],
              pvals.nonaggr = pvals,
              sel.models    = sel.model.all,
              ci.level      = ifelse(ci, ci.level, NA),
              method        = "multi.split",
              call          = match.call())
  
  class(out) <- "hdi"
  
  return(out)
}

##aggregate the ci over multiple splits for one single beta_j
aggregate.ci <- function(lci,rci,centers,
                         ses,
                         df.res,
                         gamma.min=0.05,
                         multi.corr=FALSE,
                         verbose=FALSE,
                         var,
                         s0){
  ##detect the -Inf Inf intervals
  inf.ci <- is.infinite(lci)|is.infinite(rci)
  no.inf.ci <- sum(inf.ci)##this we will use later on
  if(verbose)
    {
      print("number of Inf ci")
      print(no.inf.ci)
    }
  if((no.inf.ci == length(lci)) || (no.inf.ci >= (1-gamma.min)*length(lci)))
    {##we only have infinite ci or more than 1-gamma.min of the splits have an inf ci
      return(c(-Inf,Inf))
    }
  ##remove the infinite ci from the input
  lci <- lci[!inf.ci]
  rci <- rci[!inf.ci]
  centers <- centers[!inf.ci]
  df.res <- df.res[!inf.ci]
  ses <- ses[!inf.ci]
  s0 <- s0[!inf.ci]
  
  ci.lengths <- rci-lci
  ci.info <- list(lci=lci,
                  rci=rci,
                  centers=centers,
                  ci.lengths=ci.lengths,
                  no.inf.ci=no.inf.ci,
                  ses=ses,
                  s0=s0,
                  df.res=df.res,
                  gamma.min=gamma.min,
                  multi.corr=multi.corr)
  ##find an inside point: we need to find a point that is definitely in the confidence intervals
  inner <- find.inside.point.gammamin(low=min(centers),
                                      high=max(centers),
                                      ci.info=ci.info,
                                      verbose=verbose)

  ##inner <- max(lci)
  outer <- min(lci)
  ##finding good bounds for our bisection method
  new.bounds <- find.bisection.bounds.gammamin(shouldcover=inner,
                                               shouldnotcover=outer,
                                               ci.info=ci.info,
                                               verbose=verbose)
  inner <- new.bounds$shouldcover
  outer <- new.bounds$shouldnotcover
  
  l.bound <- bisection.gammamin.coverage(outer=outer,
                                         inner=inner,
                                         ci.info=ci.info,
                                         verbose=verbose)
  if(verbose)
    {
      print("lower bound ci aggregated is")
      print(l.bound)
    }

  ##inner <- min(rci)
  outer <- max(rci)
  new.bounds <- find.bisection.bounds.gammamin(shouldcover=inner,
                                               shouldnotcover=outer,
                                               ci.info=ci.info,
                                               verbose=verbose)
  inner <- new.bounds$shouldcover
  outer <- new.bounds$shouldnotcover
  
  u.bound <- bisection.gammamin.coverage(inner=inner,
                                         outer=outer,
                                         ci.info=ci.info,
                                         verbose=verbose)
  if(verbose)
    {
      print("upper bound ci aggregated is")
      print(u.bound)
    }
  return(c(l.bound,u.bound))
}

find.inside.point.gammamin <- function(low,
                                       high,
                                       ci.info,
                                       verbose)
{
  ##searching over the range low and high for a point that is inside
  ##does it cover gammamin
  ##we increase the granularity up until a certain level and then give up
  range.length <- 10
  test.range <- seq(low,high,length.out=range.length)
  cover <- mapply(does.it.cover.gammamin,
                  beta.j=test.range,
                  ci.info=list(ci.info=ci.info))
  while(!any(cover))
    {
      range.length <- 10*range.length
      test.range <- seq(low,high,length.out=range.length)
      cover <- mapply(does.it.cover.gammamin,
                      beta.j=test.range,
                      ci.info=list(ci.info=ci.info))
      if(range.length>10^3)
        {
          print("FOUND NO INSIDE POINT")
          print("number of splits")
          print(length(ci.info$centers))
          print("centers")
          print(ci.info$centers)
          print("ses")
          print(ci.info$ses)
          print("df.res")
          print(ci.info$df.res)
          
          stop("couldn't find an inside point between low and high. The confidence interval doesn't exist!")
        }
    }
  if(verbose)
    {
      print("Found an inside point at granularity of ")
      print(range.length)
    }
  inside.point <- min(test.range[cover])##take the smalles value that was covered, which value we take is not of importance actually
  return(inside.point)
}

find.bisection.bounds.gammamin <- function(shouldcover,
                                           shouldnotcover,
                                           ci.info,
                                           verbose)
{
  reset.shouldnotcover <- FALSE
  if(does.it.cover.gammamin(beta.j=shouldnotcover,ci.info=ci.info))
    {
      reset.shouldnotcover <- TRUE
      if(verbose)
        print("finding a new shouldnotcover bound")
      ##need to find a shouldnotcover further away from this point
      ##the direction we move in is shouldnotcover-shouldcover
      while(does.it.cover.gammamin(beta.j=shouldnotcover,ci.info=ci.info))
        {
          old <- shouldnotcover
          shouldnotcover <- shouldnotcover + (shouldnotcover - shouldcover)
          shouldcover <- old##update the should cover bound too!
        }
      if(verbose)
        {
          print("new")
          print("shouldnotcover")
          print(shouldnotcover)
          print("shouldcover")
          print(shouldcover)
        }
    }
  ##Is it possible that these get triggered consecutively?, no :0
  if(!does.it.cover.gammamin(beta.j=shouldcover,ci.info=ci.info))
    {
      if(reset.shouldnotcover)
        stop("Problem: we first reset shouldnotcover and are now resetting shouldcover, this is not supposed to happen")
      if(verbose)
        print("finding a new shouldcover bound")
      while(!does.it.cover.gammamin(beta.j=shouldcover,ci.info=ci.info))
        {##Problem, it is possible that there is no coverage!?!?!?
          ##This could be if we jump over the CI!!
          ##TODO: fix this, NEED A SMARTER WAY TO FIND AN INSIDE POINT
          old <- shouldcover
          shouldcover <- shouldcover + (shouldcover - shouldnotcover)
          shouldnotcover <- old
        }
      if(verbose)
        {
          print("new")
          print("shouldnotcover")
          print(shouldnotcover)
          print("shouldcover")
          print(shouldcover)
        }
    }
  return(list(shouldcover=shouldcover,
              shouldnotcover=shouldnotcover))
}

check.bisection.bounds.gammamin <- function(shouldcover,
                                            shouldnotcover,
                                            ci.info,
                                            verbose)
{
  if(does.it.cover.gammamin(beta.j=shouldnotcover,
                            ci.info=ci.info))
    {
      stop("shouldnotcover bound is covered! we need to decrease it even more! (PLZ implement)")
    }else{
      if(verbose)
        print("shouldnotcover bound is not covered, this is good")
    }
  
  if(does.it.cover.gammamin(beta.j=shouldcover,
                            ci.info=ci.info))
    {
      if(verbose)
        print("shouldcover is covered!, It is a good covered bound")
    }else{
      stop("shouldcover is a bad covered bound, it is not covered!")
    }
}

bisection.gammamin.coverage <- function(outer,
                                        inner,
                                        ci.info,
                                        verbose,
                                        eps.bound=10^(-7))
{
  check.bisection.bounds.gammamin(shouldcover=inner,
                                  shouldnotcover=outer,
                                  ci.info=ci.info,
                                  verbose=verbose)
  ##do.bisection
  eps <- 1
  
  while(eps > eps.bound)
    {
      ##calc on the outer + inner /2 and see if the thing covers in this
      middle <- (outer+inner)/2
      if(does.it.cover.gammamin(beta.j=middle,
                                ci.info=ci.info))
        {
          inner <- middle
        }else{
          outer <- middle
        }
      eps <- abs(inner-outer)
    }
  solution <- (inner+outer)/2
  if(verbose)
    {
      print("finished bisection")
      print("eps is")
      print(eps)
    }
  return(solution)
}

does.it.cover.gammamin <- function(beta.j,
                                   ci.info)#,
                                        #non.central.t=FALSE)
{
  if(missing(ci.info))
    stop("ci.info is missing to the function does.it.cover.gammamin")
  ##extract ci.info
  centers <- ci.info$centers
  ci.lengths <- ci.info$ci.lengths
  no.inf.ci <- ci.info$no.inf.ci
  ses <- ci.info$ses
  df.res <- ci.info$df.res
  gamma.min <- ci.info$gamma.min
  multi.corr <- ci.info$multi.corr
  s0 <- ci.info$s0
  
  ##Warning!: this is also affected by noncentral vs central t dist
  pval.rank <- rank(-abs(beta.j-centers)/(ci.lengths/2))##the rank of the pvalue in increasing order, - sign to reverse rank
  nsplit <- length(pval.rank) + no.inf.ci##the number of ci + the inf ci we left out

  gamma.b <- pval.rank/nsplit
  if(multi.corr)
    {
      if(any(is.na(s0)))
        stop("need s0 information to be able to create multiple testing corrected pvalues")
      level <- (1-0.05*gamma.b/(1-log(gamma.min)*s0))
    }else{
      level <- (1-0.05*gamma.b/(1-log(gamma.min)))
    }
  ##from the getAnywhere(confint.lm) code
  a <- (1 - level)/2
  a <- 1-a

  ##Warning!: this is also affected by noncentral vs central t dist
  ##if(non.central.t)
  ##  {
  ##    ##CHECK: is this the correct way of doing things?
  ##    ##browser()
  ##    ##non-central t-distribution! --> not symmetric anymore!!!!
  ##    lfac <- qt(1-a,df.res,ncp=(centers-beta.j)/ses)
  ##    ufac <- qt(a,df.res,ncp=(centers-beta.j)/ses)
  ##    ##Are these the correct non centers to pick? not the point where we want to test!?!?
  ##    ##TODO DEBUG
  ##    nlci <- beta.j+lfac*ses
  ##    nrci <- beta.j+ufac*ses
  ##    if(all(gamma.b <= gamma.min))
  ##      {##the fraction of non Inf ci of all splits is smaller than gamma.min
  ##        beta.is.in <- TRUE
  ##      }else{
  ##        ##should not check gamma.b < 1 because we already removed the -Inf Inf ci!
  ##        beta.is.in <- all(nlci[gamma.b > gamma.min ] <= beta.j) && all(nrci[gamma.b > gamma.min ] >= beta.j)##the centers are in the ci around beta.j/c
  ##      }
  ##  }else{
  
  ##central t-distribution, approximations are made! if df big enough this shouldn't be a problem
  fac <- qt(a,df.res)
  
  nlci <- centers - fac*ses
  nrci <- centers + fac*ses
  if(all(gamma.b <= gamma.min))
    {##the fraction of non Inf ci of all splits is smaller than gamma.min
      beta.is.in <- TRUE
    }else{
      ##should not check gamma.b < 1 because we already removed the -Inf Inf ci!
      beta.is.in <- all(nlci[gamma.b > gamma.min] <= beta.j) && all(nrci[gamma.b > gamma.min] >= beta.j)##beta_j is in the required intervals
    }
  
  return(beta.is.in)
}



