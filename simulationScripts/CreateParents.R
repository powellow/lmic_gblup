# Cleanup
rm(list = ls())

### SetUp Environment ------------------------

library("default"); #install.packages("default")
default(hist.default) = list(col="limegreen",border="white",ylab="",xlab="")
#library(devtools)
#install_bitbucket("hickeyjohnteam/AlphaSimR@master")
library(AlphaSimR)

### Population & Genetic Parameters --------------------

nParents = 2000
nChr = 10
nQtl = 1000
nSnp = 5000

initMeanG = 0
initVarG = 0.2
initVarEnv = 1.8
VarE = 1.8

### ----------------------------------------------------

macsCommand = function(size){
  paste(1.00E+08,"-t",4.00E-05,"-r",4.00E-05,"-I 1",size,
        "-eN 0.0045 1.25 -eN 0.006 2 -eN 0.0385 3 -eN 0.1135 4 -eN 0.1635 5 -eN 0.4385 10 -eN 0.5885 15 -eN 0.8385 25 -eN 8.2885 100 -s",sample.int(1e8,1))
}


FOUNDERPOP = runMacs(nInd=nParents,
                     nChr=nChr,
                     segSites=nQtl+nSnp,
                     manualCommand = macsCommand(nParents*2),manualGenLen = 1)

SP = SimParam$new(FOUNDERPOP)
SP$setTrackPed(TRUE) #creates a tracked pedigree. Usefull for retrospective problem solving
SP$setGender("yes_sys")

if(nSnp>0){
  SP$addSnpChip(nSnpPerChr=nSnp)
}

SP$addTraitAD(nQtlPerChr=nQtl,mean=initMeanG,var=initVarG, useVarA=F)

default(hist.default) = list(col="limegreen",border="white",ylab="",xlab="")

### FOUNDER POPULATION 

initPop = newPop(FOUNDERPOP)
initPop = setPheno(initPop,varE=initVarEnv)
save.image(file='Macs.RData')
