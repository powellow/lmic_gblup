# Cleanup
rm(list = ls())
load("Macs.RData")

### Input arguments ------------------------
  
args = commandArgs(trailingOnly=TRUE)
#nSiresPerGen = 250 ; 
nSiresPerGen = as.numeric(args[1])

### SetUp Environment ------------------------

library("default"); #install.packages("default")
default(hist.default) = list(col="limegreen",border="white",ylab="",xlab="")

#library(devtools)
#install_bitbucket("hickeyjohnteam/AlphaSimR@master")
library(AlphaSimR)

#setwd("/exports/cmvm/eddie/eb/groups/tier2_hickey_group/personal/OwenPhD/AlphaBayesWithFixedAndRandomEffects/ASR/")  
Sires = paste0(nSiresPerGen,"Sires")
dir.create(paste0("./",Sires)); setwd(paste0("./",Sires))
dir.create("./Sires") ; dir.create("./Dams") ; dir.create("./Sons")


### Simulation Parameters ------------------------------

BI_sires = 225
BI_dams = 1000
nInd_Sel = 80000
Sel_sires = nSiresPerGen
Sel_dams = 40000

nBIyrs = 5
nSelyrs = nBIyrs+6


# Run breeding programme --------------------------------------------------


for (year in 1:nBIyrs){
  if (year ==1){
    # Extract parents
    females = initPop[initPop@gender=="F"] 
    males = initPop[initPop@gender=="M"] 
    rm(FOUNDERPOP) #No longer needed
  }
  else{
    females = offspring[offspring@gender=="F"] 
    males = offspring[offspring@gender=="M"]  
  }
  
  maleParents = selectInd(males,BI_sires,use="pheno")
  femaleParents = selectInd(females,BI_dams,use="pheno")
  
  if (year == nBIyrs){
    sire_id = sample(maleParents@id,nInd_Sel,replace=T)
    dam_id = sample(rep(femaleParents@id,(nInd_Sel/length(femaleParents@id))),replace=F)
    mate_plan = cbind(dam_id,sire_id) ; mate_plan = as.matrix(mate_plan)
  }
  else{
    sire_id = sample(maleParents@id,nParents,replace=T)
    dam_id = sample(rep(femaleParents@id,(nParents/length(femaleParents@id))),replace=F)
    mate_plan = cbind(dam_id,sire_id) ; mate_plan = as.matrix(mate_plan)
  }
  
  n_offspring = unname(table(sire_id))
  mean(n_offspring); var(n_offspring); hist(n_offspring,main=paste("Number of Offspring Per Sire - Generation",year),xlab = "No. of Offspring", ylab = "No. of Sires")
  
  offspring = hybridCross(femaleParents,maleParents,crossPlan = mate_plan,returnHybridPop = F)
  offspring = setPheno(pop=offspring,varE=VarE,reps=1)
}

for (year in ((nBIyrs+1):nSelyrs)){
  females = offspring[offspring@gender=="F"]
  females = setPheno(pop=females,varE=VarE,reps=1) 
  males = offspring[offspring@gender=="M"]  
  males = setPheno(pop=males,varE=VarE,reps=1)

  maleParents = selectInd(males,Sel_sires,use="pheno")
  maleParents = setPheno(pop=maleParents,varE=VarE,reps=1)
  femaleParents = selectInd(females,Sel_dams,use="pheno")
  femaleParents = setPheno(pop=femaleParents,varE=VarE,reps=1)

  sire_id = sample(maleParents@id,nInd_Sel,replace=T)
  dam_id = sample(rep(femaleParents@id,(nInd_Sel/length(femaleParents@id))),replace=F)
  mate_plan = cbind(dam_id,sire_id) ; mate_plan = as.matrix(mate_plan)
  
  n_offspring = unname(table(sire_id))
  mean(n_offspring); var(n_offspring); hist(n_offspring,main=paste("Number of Offspring Per Sire - TBV Sel - Generation",year),xlab = "No. of Offspring", ylab = "No. of Sires")
  
  offspring = hybridCross(femaleParents,maleParents,crossPlan = mate_plan,returnHybridPop = F)
  offspring = setPheno(pop=offspring,varE=VarE,reps=1)
  
  if (year == nSelyrs){
    cat("Writing Records...\n")
    writeRecords(femaleParents,"Dams",1,includeHaplo=FALSE,
                 reps=1,append=F)
    writeRecords(maleParents,"Sires",1,includeHaplo=FALSE,
                 reps=1,append=F)
    writeRecords(offspring,"Sons",1,includeHaplo=FALSE,
                 reps=1,append=F)
  }
}

save.image(file='BreedingPopulation.RData')


