library(MASS)
library(ggplot2)
library(gridExtra)
library(viridis)
library(tidyverse)

targetSpecificity <- .95 # or 0.995

## Function to generate data ----

genDat <- function(n_neg, n_pos, test=c("spike","rbd")) {
  if (test=="spike") {
    # Positive
    mixes <- sample(2, n_pos, prob=c(0.6511623, 0.3488377), replace=T)
    pos <- numeric(n_pos)
    pos[mixes==1] <- rgamma(sum(mixes==1), 8.4160019, scale=0.5324562)
    pos[mixes==2] <- rlnorm(sum(mixes==2), 1.0566340, 0.1197868)
    
    # Negative
    mixes <- sample(2, n_neg, prob=c(0.1274695, 0.8725305), replace=T)
    neg <- numeric(n_neg)
    neg[mixes==1] <- rlnorm(sum(mixes==1), 0.7531151, 0.4223818)
    neg[mixes==2] <- rlnorm(sum(mixes==2), -0.1331322, 0.3230064)
  } else if (test=="rbd") {
    # Positive
    mixes <- sample(3, n_pos, prob=c(0.36340028, 0.59119533, 0.04540438), replace=T)
    
    pos <- numeric(n_pos)
    pos[mixes==1] <- rgamma(sum(mixes==1), 5.867495, scale=2.315552)
    pos[mixes==2] <- rgamma(sum(mixes==2), 13.334534, scale=1.483064)
    pos[mixes==3] <- rgamma(sum(mixes==3), 1.671092, scale=4.706356)
    
    # Negative
    mixes <- sample(3, n_neg, prob=c(0.27266311, 0.69342221, 0.03391467), replace=T)
    neg <- numeric(n_neg)
    neg[mixes==1] <- rgamma(sum(mixes==1), 6.6175403, scale=0.2899964)
    neg[mixes==2] <- rlnorm(sum(mixes==2), -0.04849108, 0.22474440)
    neg[mixes==3] <- rlnorm(sum(mixes==3), 1.4875934, 0.8820271 )
  } else {
    stop("wrong dist \n")
  }
  
  list(pos=pos, neg=neg)
}

## Functions to fit cutoffs ----

fitPareto <- function(dat, ct, Q) {
  mu <- quantile(dat,ct)
  quant <- (Q-ct)/(1-ct)
  qexp(quant, 1/(mean(dat[dat>=mu])-mu))+mu
}

estCutoffs <- function(dat, Q) {
  # Empirical
  emperical <- quantile(dat, Q)
  
  # Normal and lognormal
  norm <- qnorm(Q, mean=mean(dat), sd=sd(dat))
  ldat <- log(dat)
  lnorm <- qlnorm(Q, meanlog=mean(ldat), sdlog=sd(ldat))
  
  # MAD and log MAD
  madnorm <- qnorm(Q, mean=median(dat), sd=mad(dat))
  lmadnorm <- qlnorm(Q, meanlog=median(ldat), sdlog=mad(ldat))
  
  # Pareto 0.9 and Pareto 0.95
  expFit.9 <- fitPareto(dat,.9, Q)
  expFit.95 <- fitPareto(dat,.95, Q)
  
  c(emp=as.numeric(emperical),
    norm=as.numeric(norm), lnorm=lnorm, 
    madnorm=madnorm, lmadnorm=lmadnorm, 
    pareto.9=as.numeric(expFit.9), pareto.95=as.numeric(expFit.95))
}

## Functions to estimate prevalence ----

RGadj <- function(pos, sens, Q) {
  (pos + Q-1)/(Q+sens-1)
}

estPrev <- function(posTraining, testing, cuts, Q) {
  sens <- sapply(cuts, function(x) mean(posTraining>=x))
  positivity <- sapply(cuts, function(x) mean(testing>=x))
  RGadj(pos=positivity, sens=sens, Q)
}

## Functions to estimate sens/spec/pos/neg ----

getSens <- function(tab) {
  tab["1","1"]/(tab["1","1"]+tab["0","1"])
}

getSpec <- function(tab) {
  tab["0","0"]/(tab["0","0"]+tab["1","0"])
}

getPosPred <- function(tab) {
  tab["1","1"]/(tab["1","0"]+tab["1","1"])
}

getNegPred <- function(tab) {
  tab["0","0"]/(tab["0","0"]+tab["0","1"])
}

getTotalError <- function(tab) {
  (tab["0","0"]+tab["1","1"])/sum(tab)
}

## The simulation function ----

sim <- function(truePrev, test, n_neg, n_pos, n_test, Q) {
  # Generate data
  training <- genDat(n_neg=n_neg, n_pos=n_pos, test=test)
  testingdatsep <- genDat(n_neg=round(n_test*(1-truePrev)), 
                          n_pos=round(n_test*truePrev), 
                          test=test)
  testingdat <- c(testingdatsep$pos, testingdatsep$neg)
  
  # Estimate cutoffs
  cuts <- estCutoffs(training$neg, Q)
  
  # Specificity, sensitivity, positive predictive, negative predictive, error
  testChar <- sapply(cuts, function(x) {
    results <- ifelse(testingdat>=x, 1, 0)
    tab <- table(factor(results, levels=c(0,1)), 
                 factor(c(rep(1, length(testingdatsep$pos)), rep(0, length(testingdatsep$neg))),
                        levels=c(0,1)),
                 useNA="ifany")
    c(sens=getSens(tab),
      spec=getSpec(tab),
      pospred=getPosPred(tab),
      negpred=getNegPred(tab),
      totalerror=getTotalError(tab))
  })
  
  # Estimate Prevalence
  prevalence <- estPrev(posTraining=training$pos, 
                        testing=testingdat, 
                        cuts=cuts, Q=Q)
  
  # Test normality for hybrid methods
  hemp <- 1 # hybrid empirical
  if (shapiro.test(training$neg)$p.val > 0.05) {
    hemp <- 2
  } else if (shapiro.test(log(training$neg))$p.val > 0.05) {
    hemp <- 3
  }
  
  hp95 <- hp9 <- hemp
  ind <- which(hemp==1)
  hp9[ind] <- 6 # hybrid Pareto 0.9
  hp95[ind] <- 7 # hybrid Pareto 0.95
  
  cuts <- c(cuts, hemp=cuts[hemp], hp9=cuts[hp9], hp95=cuts[hp95])
  testChar <- cbind(testChar, testChar[,c(hemp, hp9,hp95)])
  prevalence <- c(prevalence, hemp=prevalence[hemp], hp9=prevalence[hp9], hp95=prevalence[hp95])
  
  
  
  return(list(cuts=cuts,
              testChar=testChar,
              prevalence=prevalence,
              hybrid=hemp
  ))
}

## Run it ----

allSettings <- expand.grid(c(.05, .3), c("spike","rbd"), c(50, 200), c(200, 500))
allSettings <- cbind(allSettings[,c(1:3,3,4)])

set.seed(87245)

system.time({
  simResults <- lapply(1:nrow(allSettings), function(i) {
    replicate(10000, sim(truePrev=allSettings[i,1], 
                         test=allSettings[i,2], 
                         n_neg=allSettings[i,3], 
                         n_pos=allSettings[i,4], 
                         n_test=allSettings[i,5],
                         Q=targetSpecificity),
              simplify=FALSE)
  })
})

save(simResults, allSettings, file=paste0("cutoff-simulation-",targetSpecificity,".RData"))
