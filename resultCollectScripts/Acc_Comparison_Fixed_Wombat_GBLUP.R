library("data.table"); library('default')
default(hist.default) = list(col = "limegreen",border="white"); 
#default(plot.default) = list(col = "forestgreen")

args = commandArgs(trailingOnly=TRUE)

### ------------------------------------------------------------
### ------------------- DATA ENTRY -----------------------------
### ------------------------------------------------------------

### Load Genotypes and AlphaBayes Results ----------------------

#X = fread("Genotypes.txt"); X[1:5,1:5] ; X = X[,c(-1)]; X[1:5,1:5] 

# ### Fixed Effect Model Marker Estimates -----------------------
# 
# ab_f = read.table("FixedModel/MarkerEstimate.txt"); ab_f= as.vector(ab_f$V1)
# ab_fh = read.table("FixedModel/FixedEffectEstimate.txt"); ab_fh = ab_fh$V1
# 
# ab_febv = rowSums(sweep(X,MARGIN=2,ab_f,FUN='*')); ab_febv = ab_febv-mean(ab_febv)

### Fixed Effect Model Marker Estimates -----------------------
ab_f = read.table("GBlupFixedModel/RnSoln_animal.dat",header=F,skip=1); ab_f = ab_f$V4
ab_febv = tail(ab_f,n=as.numeric(args[1]))

### Load Truth Data Results ------------------------------------

true_data = read.table("PedigreeASV.txt",header=T); true_data = tail(true_data,n=as.numeric(args[1]))
tgv = true_data$GeneticValue; mean(tgv); tgv = tgv - mean(tgv)
herd_data = true_data[order(true_data$HerdID),]
thv = unique(herd_data$HerdEffect)


### ------------------------------------------------------------
### ------- CALCULATE EBV ACCURCIES & BIAS ---------------------
### ------------------------------------------------------------

 ### Fixed Model ---------------------------------------------
 
print(paste("Fixed Model : Accuracy of WOMBAT GEBVs is",round(cor(ab_febv,tgv),digits=3),sep=" ")) ; plot(ab_febv,tgv,main="WOMBAT GEBV vs TBV : Fixed Model",pch=16)
abline(lm(tgv~ab_febv),lty=2,col="limegreen")
 
f_acc = round(cor(ab_febv,tgv),digits=3)
fbv_reg= lm(tgv~ab_febv)
fbv_bias = fbv_reg$coefficients[2]-1
fbv_int = fbv_reg$coefficients[1]
# 
#summary(tgv)
max_tgv_index = which(tgv>0.229)
max_tgv = tgv[max_tgv_index]
max_febv = ab_rebv[max_tgv_index]
UQ_fcor = cor(max_tgv,max_febv)
# 
min_tgv_index = which(tgv>-0.242)
min_tgv = tgv[min_tgv_index]
min_febv = ab_rebv[min_tgv_index]
LQ_fcor = cor(min_tgv,min_febv)

### Random Model --------------------------------------------

#print(paste("Random Model : Accuracy of AB EBVs is",round(cor(ab_rebv,tgv),digits=3),sep=" ")) ; plot(ab_rebv,tgv,main="WOMBAT EBV vs TBV : Random Model",pch=16)
#abline(lm(tgv~ab_rebv),lty=2,col="limegreen")

#r_acc = round(cor(ab_rebv,tgv),digits=3)
#rbv_reg= lm(tgv~ab_rebv)
#rbv_bias = rbv_reg$coefficients[2]-1
#rbv_int = rbv_reg$coefficients[1]

#summary(tgv)
#max_tgv_index = which(tgv>0.229)
#max_tgv = tgv[max_tgv_index]
#max_rebv = ab_rebv[max_tgv_index]
#UQ_rcor = cor(max_tgv,max_rebv)

#min_tgv_index = which(tgv>-0.242)
#min_tgv = tgv[min_tgv_index]
#min_rebv = ab_rebv[min_tgv_index]
#LQ_rcor = cor(min_tgv,min_rebv)
#



