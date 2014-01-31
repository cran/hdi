print.hdi <- function(x, ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  3 Apr 2013, 15:12

  if(is.list(x))
    method <- x$method
  else
    method <- ""

  #####################
  ## Multi-Splitting ##
  #####################
  
  if(method == "multi.split"){
    cat("alpha = 0.01:")
    cat(" Selected predictors:", which(x$pval <= 0.01), "\n")
    cat("alpha = 0.05:")
    cat(" Selected predictors:", which(x$pval <= 0.05), "\n")
    cat("------\n")
    cat("Familywise error rate controlled at level alpha.\n")
  }
  
  ###############
  ## Stability ##
  ###############
  else if(method == "stability"){
    cat("Selected predictors:\n")
    cat("--------------------\n")
    if(is.null(x$select)){
      cat("none\n")
    }else{
      print(x$select)
    }
   
    ##for(i in 1:length(x$select)){
    ##  cat("EV = ", x$EV[[i]], ":", sep = "")
    ##  cat(" Selected predictors:", x$select[[i]], "\n")
    #}
    cat("--------------------\n")
    cat("Expected number of false positives controlled at level", x$EV, "\n")
  }

  #########################
  ## Cluster Lower-Bound ##
  #########################
  else if(method == "clusterLowerBound"){
    cat("\n lower l1-norm group bounds in a hierarchical clustering ")
    cat("\n lower bound on l1-norm of all regression coefficients:",
        signif(max(x$lowerBound),4))
    
    if(sum(x$isLeaf & x$lowerBound>0)==1){
      tmp <- sum(x$noMembers[which(x$isLeaf & x$lowerBound > 0)])
      cat("\n only 1 significant non-overlapping cluster with", tmp ,
          if(tmp == 1) "member" else "members")
    }else{
      cat("\n number of non-overlapping significant clusters       :",
          sum(x$isLeaf & x$lowerBound > 0))
      tmp <- range(x$noMembers[which(x$isLeaf& x$lowerBound > 0)])
      if(diff(tmp) == 0){
        cat("\n with", tmp[1],
          if(tmp[1] == 1) "member in each non-overlapping cluster" else
            "members in each non-overlapping cluster")
      }else{
        cat("\n with", paste(tmp, collapse = " up to "), "members each")
      }
    }
    cat("\n ")
  }

  ######################
  ## Other situations ##
  ######################
  else{
    print.default(x, ...)
  }
}

plot.clusterLowerBound <- function(x, cexfactor = 1, yaxis = "members",
                                   col = NULL, ...){
  hh <- x$noMembers
  hh2 <- sqrt(x$lowerBound)
  
  if(yaxis!="members"){
    hh2 <- sqrt(x$noMembers)
    hh  <- x$lowerBound
  }
  
  hh2 <- hh2 / max(hh2)
  xvec <- x$position
  
  col <- if(yaxis == "members") rgb(0.8,0.2,0.2,0.6) else rgb(0.2,0.2,0.8,0.6)

  plot(xvec, hh, cex = 1, axes = FALSE, xlab = "",
       ylab = if(yaxis == "members") "cluster size" else "lower l1-norm bound",
       pch = 20, col = "white", log = if(yaxis == "members") "y" else "")
    
  axis(2)
  coll <- rgb(0.1, 0.1, 0.1, 0.7)
  
  for(k in length(x$position):1){
    if((nm <- x$leftChild[k]) > 0){
      lines(c(xvec[k],  xvec[nm]), c(hh[k], hh[k]),  col=coll)
      lines(c(xvec[nm], xvec[nm]), c(hh[k], hh[nm]), col=coll)
    }
    if((nm <- x$rightChild[k]) > 0){
      lines(c(xvec[k],  xvec[nm]), c(hh[k],hh[k]),  col = coll)
      lines(c(xvec[nm], xvec[nm]), c(hh[k],hh[nm]), col= coll)
    }
  }
  points(xvec, hh, cex = sqrt(hh2) * 4 * cexfactor, xlab = "", col = "white",
         pch = 20)
  points(xvec, hh, cex = sqrt(hh2) * 4 * cexfactor, xlab = "", col = col,
         pch = 20)
}




