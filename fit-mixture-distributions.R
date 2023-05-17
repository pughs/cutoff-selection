library(ltmix)

load("./mixture-model-data.RData")

# Spike ----
## Spike positives ----

ot <- ltmmCombo(spikePos, G=3) 
ot <- ltmmCombo(spikePos, G=2)

set.seed(29253252)
ot <- ltmm(spikePos, G=2, distributions=c('gamma','lognormal'), trunc=0)
hist(spikePos, prob=T, breaks=30, main="2 comp spike pos")
curve(ot$fitted_pdf(x),add=T,col="red")

## Spike negatives ----

ot <- ltmmCombo(spikeNeg, G=3) 
ot <- ltmmCombo(spikeNeg, G=2)

set.seed(29253252)
ot <- ltmm(spikeNeg, G=2, distributions=c('lognormal','lognormal'), trunc=0)
hist(spikeNeg, prob=T, breaks=30, main="2 comp spike neg")
curve(ot$fitted_pdf(x),add=T,col="red")

# RBD ----
## RBD positives ----

ot <- ltmmCombo(rbdPos, G=3)
ot <- ltmmCombo(rbdPos, G=2)

set.seed(23893) 
ot <- ltmm(rbdPos, G=3, distributions=c('gamma','gamma','gamma'), trunc=0)
hist(rbdPos, prob=T, breaks=20, main="RBD 3 comp pos")
curve(ot$fitted_pdf(x),add=T,col="red")

## RBD negatives ----

ot <- ltmmCombo(rbdNeg, G=3)
ot <- ltmmCombo(rbdNeg, G=2)

set.seed(23893) 
ot <- ltmm(rbdNeg, G=3, distributions=c('gamma','lognormal','lognormal'), trunc=0)
hist(rbdNeg, prob=T, breaks=100, main="RBD 3 comp neg", xlim=c(0,6))
curve(ot$fitted_pdf(x),add=T,col="red")
