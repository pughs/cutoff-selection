library(knitr)
library(gridExtra)
library(MASS)
library(ggplot2)
library(xtable)
library(viridis)

targetSpecificity <- .95 # or 0.995

## Read in data ----

load("training-data.RData") 
load("testing-data.RData")

negs <- dat[dat$neut=="Negative" & !is.na(dat$neut),]
pos <- dat[dat$neut=="Positive" & !is.na(dat$neut),]

## Beeswarm plot of the data ----

spike <- ggplot(dat[!is.na(dat$neut),]) + 
  geom_jitter(data=dat[!is.na(dat$neut),], aes(x=neut, y=Spike.PN, col=neut), alpha=.5, 
              position=position_jitter(width=.2)) +
  labs(y="Spike P/N", x="Neutralization assay result", color="Neut result") +
  theme_bw()

rbd <- ggplot(dat[!is.na(dat$neut),]) + 
  geom_jitter(data=dat[!is.na(dat$neut),], aes(x=neut, y=RBD.PN, col=neut), alpha=.5, 
              position=position_jitter(width=.2)) +
  labs(y="RBD P/N", x="Neutralization assay result", color="Neut result") +
  theme_bw()

grid.arrange(spike, rbd, ncol=2)

## Estimate cutoffs functions ----

fitPareto <- function(dat, ct, Q) {
  mu <- quantile(dat,ct)
  quant <- (Q-ct)/(1-ct)
  qexp(quant, 1/(mean(dat[dat>=mu])-mu))+mu
}

estCutoffs <- function(dat) {
  Q <- targetSpecificity
  emperical <- quantile(dat, Q) # Empirical
  norm <- qnorm(Q, mean(dat), sd(dat)) # Normal 
  lnorm <- qlnorm(Q, mean(log(dat)), sd(log(dat))) # Lognormal
  madnorm <- qnorm(Q, median(dat), mad(dat)) # MAD
  lmadnorm <- qlnorm(Q, median(log(dat)), mad(log(dat))) # Log MAD
  expFit.9 <- fitPareto(dat,.9, Q) # Pareto 0.9
  expFit.95 <- fitPareto(dat,.95, Q) # Pareto 0.95
  c(emp=as.numeric(emperical),
    norm=as.numeric(norm), lnorm=lnorm,  
    madnorm=madnorm, lmadnorm=lmadnorm, 
    pareto.9=as.numeric(expFit.9), pareto.95=as.numeric(expFit.95))
}

# Get cutoffs ----

cutSpike <- estCutoffs(negs$Spike.PN)
cutRBD <- estCutoffs(negs$RBD.PN)

a <- rbind(cutSpike, cutRBD)
colnames(a) <- c("Emp", "Normal", "Log Norm", "MAD", "Log MAD", "Par 10%", "Par 5%")
rownames(a) <- c("Spike", "RBD")

xtable(a, digits=1)

## Estimate sensitivity ----

## Spike ##
posOfTotal <- sapply(cutSpike, function(x) mean(dat$Spike.PN>x))
sensEst <- sapply(cutSpike, function(x) mean(pos$Spike.PN>x))
sensSpike <- sensEst

## RBD ##
posOfTotal <- sapply(cutRBD, function(x) mean(dat$RBD.PN>x))
sensEst <- sapply(cutRBD, function(x) mean(pos$RBD.PN>x))
sensRBD <- sensEst

a <- rbind(sensSpike, sensRBD)
colnames(a) <- c("Emp","Normal","Log Norm", "MAD", "Log MAD", "Par 10%", "Par 5%")
rownames(a) <- c("Spike", "RBD")
round(a,2)

## Estimate specificity ----

specSpike <- sapply(cutSpike, function(x) mean(negs$Spike.PN<=x))
specRBD <- sapply(cutRBD, function(x) mean(negs$RBD.PN<=x))

a <- rbind(specSpike, specRBD)
colnames(a) <- c("Emp","Normal","Log Norm", "MAD", "Log MAD", "Par 10%", "Par 5%")
rownames(a) <- c("Spike", "RBD")
round(a,2)

## Prevalence of testing data ----

RGadj <- function(pos, sens) {
  (pos + .995-1)/(.995+sens-1)
}

## Spike ##
allSpike <- testingSpike
posOfTotal <- sapply(cutSpike, function(x) mean(allSpike>x))
rgSpike <- RGadj(posOfTotal, sensSpike)

## RBD ##
allRBD <- testingRBD
posOfTotal <- sapply(cutRBD, function(x) mean(allRBD>x))
rgRBD <- RGadj(posOfTotal, sensRBD)

## Table
a <- rbind(rgSpike, rgRBD)
colnames(a) <- c("Empirical", "Normal", "Log Normal", "MAD Normal", "Log MAD Normal", "Pareto 10%", "Pareto 5%")
rownames(a) <- c("Spike", "RBD")

round(a,2)

## Plot cutoffs with testing data too ----

thecols <- magma(20)

spike_all_data <- data.frame(dat[,c("Spike.PN","RBD.PN","neut")])
spike_all_data <- rbind(spike_all_data,
                        data.frame(Spike.PN=testingSpike, RBD.PN=testingRBD, neut="Testing"))
spike_all_data$neut <- factor(spike_all_data$neut)

plotRes2 <- function(x, test=c("RBD.PN", "Spike.PN")) {
  pltDat <- data.frame(val=x, 
                       method=c("Empirical",rep("Normal",2),rep("MAD normal",2), rep("Pareto",2)),
                       spec=as.factor(c(1,1:2,1:2,1,1)))
  levels(pltDat$spec) <- c("No transformation","Log-transformation")
  
  lab <- ifelse(test=="RBD.PN", "RBD P/N ratio", "Spike P/N ratio")
  
  test <- spike_all_data[!is.na(dat$neut),test]
  
  pltDat[6, "method"] <- "Pareto 0.9"
  pltDat[7, "method"] <- "Pareto 0.95"
  
  if (targetSpecificity==0.95) {
    pltDat <- pltDat[-7,]
  }
  
  set.seed(1) # For the jitter
  
  ggplot() + 
    geom_jitter(data=spike_all_data[!is.na(spike_all_data$neut),], aes(x=neut, y=test), alpha=.5, 
                position=position_jitter(width=.2)) +
    labs(y=lab, x="Data type", color="Method", linetype=" ") +
    geom_hline(aes(yintercept=val, color=method, linetype=spec), data=pltDat, size=.7) + 
    scale_colour_manual(values = c("#505BAC", "#E9616E","#FBBA72", "#94C9A9", "#4B816F"),
                        limits = c("Empirical","Normal","MAD normal", "Pareto 0.9", "Pareto 0.95")) +
    ggtitle(paste0("Target specificity=",targetSpecificity)) +
    theme_bw() + 
    theme(plot.title=element_text(size=12, margin=margin(t=15,b=-13)))
}

plotRes2(cutSpike, "Spike.PN")
# ggsave("spike-plot.pdf", width=4.8, height=4)
plotRes2(cutRBD, "RBD.PN")
# ggsave("rbd-plot.pdf", width=4.8, height=4)

## Test normality of training data ----

shapiro.test(negs$Spike.PN)
shapiro.test(log(negs$Spike.PN))
shapiro.test(negs$RBD.PN)
shapiro.test(log(negs$RBD.PN))
