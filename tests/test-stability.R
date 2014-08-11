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
##- source("stability.R")
##- source("helpers.R")
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

x.use <- x[,1:200]

#########################
## Stability selection ##
#########################

x <- matrix(rnorm(100*1000), nrow = 100, ncol = 1000)
y <- x[,1] + x[,2] + rnorm(100)

fit.stab <- stability(x, y, EV = 1)
fit.stab
fit.stab$freq[1:10] ## selection frequency of the first 10 predictors

##- fit.stab2 <- stability(x, y, EV = 1, parallel = TRUE, trace = TRUE)
##- fit.stab2
##- 
##- stopifnot(all.equal(fit.stab$select, fit.stab2$select))
##- 
##- fit.stab3 <- stability(x, y, EV = 500, parallel = TRUE)
##- fit.stab3






