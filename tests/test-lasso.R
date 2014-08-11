#####################################
## Load stuff for testing purposes ##
#####################################

##- library(glmnet)
##- library(scalreg)
##- 
##- setwd("/u/meierluk/R/Pkgs/hdi/pkg/R")
##- 
##- source("hdi.R")
##- source("methods.R")
##- source("lasso-proj.R")
##- source("helpers.R")
##- source("helpers.nodewise.R")
##- 
##- ##########################
##- ## Load riboflavin data ##
##- ##########################
##- 
##- load("/u/meierluk/research/annualReview/Rcode/dsmN71.rda")

library(hdi)

data(riboflavin)

x <- riboflavin[,-1]
y <- riboflavin[,1]

dim(x)
##- [1]   71 4088
length(y)
##- [1] 71

x.use <- x[,1:50]

######################
## Lasso projection ##
######################

set.seed(3) ## because of cv
fit.lasso <- lasso.proj(x = x.use, y = y)

## Check standardization
set.seed(3)
fit.lasso2 <- lasso.proj(x = 2 + 4 * x.use, y = y)
stopifnot(all.equal(fit.lasso$pval, fit.lasso2$pval))

stopifnot(all.equal(max(abs(range(fit.lasso$bhat / fit.lasso2$bhat - 4))), 0))

## confidence intervals
ci.lasso  <- confint(fit.lasso, level = 0.95)
ci.lasso2 <- confint(fit.lasso2, level = 0.95)

stopifnot(all.equal(ci.lasso, ci.lasso2 * 4))


#########
## GLM ##
#########

