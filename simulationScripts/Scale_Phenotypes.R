#!/exports/cmvm/eddie/eb/groups/hickey_group/microsoft-r-open-MRO-3.4.0/source/bin/Rscript

#install.packages("default",repos='http://cran.rstudio.com/' ) ;
### Load Environment --------------------------------------------------------
library("default");
default(hist.default) = list(col="limegreen",border="white",ylab="",xlab="")
#install.packages('data.table',repos='http://cran.rstudio.com/'); 
library(data.table)
#install.packages("readr",repos='http://cran.rstudio.com/')
#install.packages("tibble",repos='http://cran.rstudio.com/')
options(scipen=999)

### Parse Arguments ---------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
Desired_Final_Var = args[1]; Desired_Final_Var = strsplit(Desired_Final_Var, "_")
Desired_Final_Avar = as.numeric(unlist(Desired_Final_Var)[1])
Desired_Final_Hvar = as.numeric(unlist(Desired_Final_Var)[2])

nSiresPerGen=as.numeric(args[2])
nInd=as.numeric(args[3])
nIndPerHerd = c(1,2,4,8,16,32)


Sires = paste0(c(nSiresPerGen,"Sires"))

Genotype = fread("../../Sons/genotype.txt") 

########################################## LOAD IN PEDIGREE & IDENTIFY GRANDPARENTAL, PARENTAL & PROGENY GENERATIONS #########################################
print(getwd())

for (HS in nIndPerHerd){

  DamPed = read.table("../info.txt_subsetted",header=F); names(DamPed) = c("id","mother","father","reps","fixEff") 
  counts = length(DamPed$id) #counts will be used later on in script to properly assign Herd Effects.
  ### Trace Offspring of Dams
  Ped = read.table("../../Sons/info.txt",header=T); names(DamPed) = names(Ped)

  off_index = which(!is.na(match(DamPed$id,Ped$mother))); d_ind = DamPed$id[off_index]; off_index = match(d_ind,Ped$mother)
  min(off_index); max(off_index)
  Sons = Ped[off_index,]
  
  Sires = unique(Sons$father)
  Dams = unique(Sons$mother)

  #### Set Initial Parameters ####
  nInd = nrow(DamPed)
  print(paste("nInd =",nInd,sep=" ")) 
  nSires = length(unique(Sires))
  print(paste("nSires =",nSires,sep=" "))
  nDams = length(unique(Dams))
  print(paste("nDams =",nDams,sep=" "))
  nHerds = floor(nInd/HS)
  print(paste("nHerds =",nHerds,sep=" "))
  nSons = nrow(Sons)
  print(paste("nSons=",nSons,sep=" "))

  ########################################## ESTIMATE ADDITIVE GENETIC & PHENOTYPIC VARIANCES IN THE PROGENY POPULATION #########################################

  ### Read in phenotype and genetic value data
  SonsPheno = read.table("../../Sons/pheno.txt",header=F) ; colnames(SonsPheno) = c("Pheno") ; SonsPheno = SonsPheno$Pheno[off_index] #reorder phenotypes to match the sampled dataframe
  Pheno = read.table("../pheno.txt_subsetted",header=F) ; colnames(Pheno) = c("Pheno") 
  SonsGV = read.table("../../Sons/gv.txt",header=F) ; colnames(SonsGV) = c("GV") ; SonsGV = SonsGV$GV[off_index] #reorder genetic values to match the sampled dataframe
  GV = read.table("../gv.txt_subsetted",header=F) ; colnames(GV) = c("GV") 

  Dam = data.frame(cbind(DamPed,GV,Pheno)) ; rm(Pheno) ; rm(GV) #combine all values into a single dataframe
  Sons = data.frame(cbind(Sons,SonsGV,SonsPheno)) ; rm(SonsPheno) ; rm(SonsGV) #combine all values into a single dataframe
  ## Calculate Sample Variances
  Dam_BV = as.numeric(as.character(Dam$GV)) # list of Dam genetic values
  Dam_Avar = var(as.numeric(as.character(Dam$GV))) # variance of Dam genetic values
  Dam_Pvar = var(as.numeric(as.character(Dam$Pheno))) # variance of Dam phenotype values
  Dam_Herit = Dam_Avar/Dam_Pvar # Dam sample heritability

  print(paste("Dam Population additive genetic is",Dam_Avar,sep=" "))
  print(paste("Dam Population phenotypic variance is",Dam_Pvar,sep=" "))
  print(paste("Dam Population heritability is",Dam_Herit,sep=" "))

  ## Calculate Sample Variances
  Sons_BV = Sons$SonsGV # list of Offspring genetic values
  Sons_Avar = var(Sons$SonsGV) # variance of Offspring genetic values
  Sons_Pvar = var(as.numeric(as.character(Sons$Pheno))) # variance of Offspring phenotype values
  Sons_Herit = Sons_Avar/Sons_Pvar # Offspring sample heritability

  print(paste("Sons Population additive genetic is",Sons_Avar,sep=" "))
  print(paste("Sons Population phenotypic variance is",Sons_Pvar,sep=" "))
  print(paste("Sons Population heritability is",Sons_Herit,sep=" "))

  ########################################## VARIANCE ESTIMATES - FINISHED ########################################################################################

  ########################################## SCALE & STANDARDIZE GENETIC VARIANCE  ########################################################################################

  mean = mean(as.numeric(as.character(Dam$GV)))
  new_mean = mean
  Scaled_Dam_BV<-(sqrt(Desired_Final_Avar)*((Dam_BV-mean(Dam_BV))/sd(Dam_BV)))+new_mean
  New_Dam_Avar = var(Scaled_Dam_BV)

  print(paste("ReScaled Additive Genetic variance of Dams is",New_Dam_Avar,sep=" "))

  mean = mean(as.numeric(Sons$SonsGV))
  new_mean = mean
  Scaled_Sons_BV<-(sqrt(Desired_Final_Avar)*((Sons_BV-mean(Sons_BV))/sd(Sons_BV)))+new_mean
  New_Sons_Avar = var(Scaled_Sons_BV)
  Scaled_Sons_BV = Scaled_Sons_BV

  print(paste("ReScaled Additive Genetic variance of Sons is",New_Sons_Avar,sep=" "))


  ########################################### SCALE & STANDARDIZE GENETIC VARIANCE - FINISHED ########################################################################################

  ########################################### SAMPLE HERD EFFECT VARIANCE ########################################################################################

  mean=0 # set mean herd effect to zero.
  sd_hVar=sqrt(Desired_Final_Hvar) # calculate the standard deviation of the desired herd effect variance
  Herd_Effect = rnorm(nHerds,mean,sd_hVar) # sample herd effects from a normal distribution (mean = 0, var = desired[supplied on command line])

  ########################################## ASSIGN IDS TO HERDS ########################################################################################

  #Script to assign Dams to Herds, when there are variable Herd sizes 

  HerdAssign <- function(nInd,nHerds,HS){
  if (nHerds<=nInd){ 
   IndPerHerd = rpois(nHerds,HS); print(paste("Individuals sampled is",sum(IndPerHerd),sep=" ")) #sample a poisson distribution of HerdSizes, report final population size
   discrep = sum(IndPerHerd)-nInd 
   if (discrep>0){ #if too many Dams were sampled
     counted = which(IndPerHerd!=0) #idenitify Herds that have at least one Dam assigned to it
     if (length(counted)<abs(discrep)){ #if number of Herds with a Dam assigned is less than the number of Dams to be unassigned... 
       dHerds = sample(counted,abs(discrep),replace = T); t = table(dHerds) #sample the Herds, n = number of Dams to be unassigned, with replacement. Store the Herd IDs in table t
     }
     if (length(counted)>=abs(discrep)){ #if number of Herds with a Dam assigned is greater than the number of Dams to be unassigned...
       dHerds = sample(counted,abs(discrep),replace = F); t = table(dHerds) #sample the Herds, n = number of Dams to be unassigned, without replacement. Store the Herd IDs in table t
     }
     IndPerHerd[as.numeric(names(t))] = IndPerHerd[as.numeric(names(t))] - as.vector(t) #Calculate new numbers of individuals per Herd by subtracting the sampled Herd IDs
     while (all(IndPerHerd>=0)==FALSE){ # If this creates negative numbers of individuals for some Herd IDs
       IndPerHerd[IndPerHerd<0] = 0 ; negs = which(IndPerHerd<0) #Set Herds with negative numbers of individuals to 0 individuals & Check for Herds still with Negative counts
       counted = which(IndPerHerd!=0) #Re-check for Herds with numbers of individuals assigned greater than 0
       dHerds = sample(counted,length(negs),replace = F); t = table(dHerds) #sample the Herds, n = number of Dams to be unassigned, without replacement. Store the Herd IDs in table t
       IndPerHerd[as.numeric(names(t))] = IndPerHerd[as.numeric(names(t))] - as.vector(t) #Calculate new numbers of individuals per Herd by subtracting the sampled Herd IDs
       if(all(IndPerHerd>=0)==TRUE) break # repeat until no negative nmber of individuals per herd exist
     }
   }
   if (discrep<0){ #if too few Dams were sampled
     dHerds = sample(1:nHerds,abs(discrep),replace = T); t = table(dHerds) #sample the Herds, n = number of Dams to be assigned, with replacement. Store the Herd IDs in table t
     IndPerHerd[as.numeric(names(t))] = IndPerHerd[as.numeric(names(t))] + as.vector(t) #Calculate new numbers of individuals per Herd by adding the sampled Herd IDs
   }
    print(paste("Adjusted Individuals sampled is",sum(IndPerHerd),sep=" ")) #Print out the new total number of individuals
    Herds = sample(1:nHerds) #Generate Herd IDs
    Dam_Herd = sample(rep(Herds,IndPerHerd)); print(paste("Individuals Assigned to Herds is",length(Dam_Herd),sep=" ")) #Assign nummber of individuals per herd to each Herd
    if(length(Dam_Herd)!=nInd){stop("Incorrect Herd Assignment")} #If number of individuals sampled is greater than the number of Dams. Stop the Code!
    }
    return(Dam_Herd)
   }

  if (nHerds<=nInd){
   Dam_Herd = HerdAssign(nInd,nHerds,HS)
  }

  #if (nHerds==nInd){ #If the number of Herds is equal to the number of Dams
  #print("About to Sample 1 Dam Per Herd")
  #Herds = sample(1:nHerds) #generate & randomise Herd IDs
  #Dam_Herd = sample(Herds,nInd,replace=F) # sample Herd IDs, without replacement
  #print("Sampled 1 Dam per Herd")
  #}
  print("Herds Assigned")
  ########################################### SCALE HERD EFFECT VARIANCE ########################################################################################

  ### Rescale Effects - Sampling of nInd Per Herd and Herd IDs, will have created undesired effects
  mean=mean(Herd_Effect)
  sd_hVar=sqrt(Desired_Final_Hvar)
  new_mean = mean
  Scaled_Herd_Effect<-(sqrt(Desired_Final_Hvar)*((Herd_Effect-mean(Herd_Effect))/sd(Herd_Effect)))+new_mean
  New_Dam_Hvar = var(Scaled_Herd_Effect)

  ########################################### ASSIGN SCALED HERD EFFECTS TO IDS ########################################################################################

  Scaled_Dam_Herd_Effect = Scaled_Herd_Effect[Dam_Herd] #Assign Herd Effects for the Dams based on HerdID index
  
  print(paste("Number of Herds is",length(unique(Dam_Herd))))

  print(paste("Herd variance is",New_Dam_Hvar,sep=" "))
  print(paste("Mean Individuals Per Herd is",mean(unname(table(Dam_Herd))),sep=" "))
  print(paste("Variance of Individuals Per Herd is",var(unname(table(Dam_Herd))),sep=" "))

  hist(unname(table(Dam_Herd)),main = "Distribution of No. Ind Per Herd",xlab = "No. Ind Per Herd",ylab = "No. Herds"); mtext(paste("Mean =",round(mean(unname(table(Dam_Herd))),digits=3),",","Variance =",round(var(unname(table(Dam_Herd))),digits=3)), side=3); 
  dev.copy(pdf,'Distribution_of_No._Ind_Per_Herd.pdf'); dev.off()

  ########################################### SCALE RESIDUAL VARIANCE  ########################################################################################

  #Sample residual variance
  Desired_Final_Evar= 1 
  sd_vare=sqrt(Desired_Final_Evar)
  mean = 0
  Final_E = rnorm(nInd,mean,sd_vare)
  mean=mean(Final_E)
  new_mean=mean
  Final_Scaled_E = (sqrt(Desired_Final_Evar)*((Final_E-mean(Final_E))/sd(Final_E)))+new_mean
  New_Dam_Evar = var(Final_Scaled_E)
  Dam_Residual = Final_Scaled_E

  print(paste("Final Population Residual Variance is",New_Dam_Evar,sep=" "))

  New_Dam_Pvar = New_Dam_Avar + New_Dam_Hvar + New_Dam_Evar

  #Calculate phenotypes of Cows
  Dam_Phenotype=Scaled_Dam_BV+Scaled_Dam_Herd_Effect+Dam_Residual
  #Calculate True Breeding Values for Calves
  Dam_GV=Scaled_Dam_BV
  Dam_Heritability = New_Dam_Avar/New_Dam_Pvar
  Dam_Herd_Herit = New_Dam_Hvar/New_Dam_Pvar

  print(paste("Dam Population Narrow Sense Heritability is",Dam_Heritability,sep=" "))
  print(paste("Dam Population Herd Variance Ratio is",Dam_Herd_Herit,sep=" "))

  ################################### BUILD FINAL DATA FILES ###################################################

  Phenotype = Dam_Phenotype ; Dam = cbind(Dam,Phenotype)
  colnames(Dam) = c("ID","mother","father","HerdID","HerdEffect","GeneticValue","OriginalPheno","Phenotype")
  Dam$HerdID = Dam_Herd
  Dam$HerdEffect = Scaled_Dam_Herd_Effect
  Dam$GeneticValue = Scaled_Dam_BV
  Data = Dam[,-c(7)]

  Sons = Sons[,-c(4,5)]; colnames(Sons) = c("ID","mother","father","GeneticValue")
  Sons$GeneticValue = Scaled_Sons_BV
  SonData = Sons

  write.table(Data,paste0(HS,"/PedigreeAndScaledValues.txt"),col.names=T,row.names=F,quote = F)
  dir.create(paste0(HS,"/Sons"))
  write.table(SonData,paste0(HS,"/Sons","/PedigreeAndScaledValues.txt"),col.names=T,row.names=F,quote = F)

  SonsGeno = Genotype[off_index,] 
  GenoID = cbind(SonData[,1],SonsGeno)
  GenoID$ID <-format(GenoID$ID, scientific = FALSE)
  write.table(GenoID,paste0(HS,"/SonsGenotype.txt"),col.names=F,row.names=F,quote = F)

  SonData = Sons[,c(1,4)]
  SonData$GeneticValue=format(SonData$GeneticValue, scientific = FALSE)
  write.table(SonData,paste0(HS,"/SonsGV.txt"),col.names=F,row.names=F,quote = F)
}
