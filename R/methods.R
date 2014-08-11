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
    
    if(sum(x$isLeaf & x$lowerBound > 0)==1){
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

plot.clusterGroupBound <- function(x, cexfactor = 1, yaxis = "members",
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

plot.clusterGroupTest <- function(x, alpha = 0.05, ...)
{
  ## Purpose:
  ## facilitate the creation of clusterGrouptest based on the
  ## group.testing.function
  ## Author: Ruben Dezeure, original version 26th of May 2014
  
  p <- length(x$hh$order)
  
  ## need to change the ordering for the plot function plot.clusterGroupBound
  clusters   <- x$clusters
  lowerbound <- ifelse(x$pval <= 0.05, 1, 0)
  leftChild  <- x$leftChild
  rightChild <- x$rightChild
  hh         <- x$hh
  
  ## order clustering according to increasing sizes
  orig.indexing <- 1:length(clusters)
  sizes         <- unlist(lapply(clusters,FUN=length))
  new.indexing  <- sort(sizes,index.return=TRUE)$ix
  clusters      <- clusters[new.indexing]
  lowerbound    <- lowerbound[new.indexing]
  
  ## need to reorganise this reordering of leftchild rightchild better!
  w.child               <- leftChild[new.indexing] > 0
  n.leftChild           <- rep(-1, length(leftChild))
  n.leftChild[w.child]  <- match(leftChild[new.indexing][w.child], new.indexing)
  w.child               <- rightChild[new.indexing] > 0
  n.rightChild          <- rep(-1, length(rightChild))
  n.rightChild[w.child] <- match(rightChild[new.indexing][w.child], new.indexing)
  
  ## get the positioning
  ord      <- hh$order##ordering from the cluster
  position <- numeric(length(clusters))
  for(i in 1:length(clusters)){
    position[i] <- mean(((1:length(ord)))[ord %in% clusters[[i]]] / p)
  }
  ## Get left child and right child of each signif node from the clustering
  
  input <- list("noMembers"  = unlist(lapply(clusters,FUN=length)),
                ## A vector containing the number of members in each group
                "lowerBound" = lowerbound,
                ## The lower bound on the l1-norm in each group
                "position"   = position,
                ## The position on the x-axis of each group (used for plotting)


                "leftChild"  = n.leftChild,
                ## Gives the index of the group that corresponds to the left
                ## child node in the tested tree (negative values correspond
                ## to leaf nodes)
                "rightChild" = n.rightChild)
                ## Same as leftChild for the right child of each node

  plot.clusterGroupBound(input,...)
}

confint.hdi <- function(object, parm, level = 0.95, ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 27 Jun 2014, 11:30

  if(!object$method %in% c("lasso.proj", "ridge.proj", "multi.split"))
    stop("Not supported object type")
  
  obj <- switch(object$method,
                "lasso.proj"  = object$bhat,
                "ridge.proj"  = object$bhat,
                "multi.split" = object$pval)
  
  if(is.null(names(obj)))
    pnames <- 1:length(obj)
  else
    pnames <- names(obj)
      
  if(missing(parm))
    parm <- pnames
  else if (is.numeric(parm))
    parm <- pnames[parm]

  if(object$method %in% c("lasso.proj", "ridge.proj")){
    quant <- qnorm(1 - (1 - level) / 2)

    if(object$method == "ridge.proj")
      delta <- object$delta[parm]
    else
      delta <- rep(0, length(parm))
        
    add   <- object$se[parm] * quant + object$se[parm] * delta
    m     <- cbind(object$bhat[parm] - add, object$bhat[parm] + add)
  }
  else if(object$method == "multi.split"){
    if(is.na(object$ci.level))
      stop("No confidence interval information available in fitted object")
    
    m <- cbind(object$lci[parm], object$uci[parm])
    
    if(level != object$ci.level)
      warning("Ignoring argument level, using ci.level of fitted object")
  }

  colnames(m) <- c("lower", "upper")
  rownames(m) <- parm

  return(m)
}

