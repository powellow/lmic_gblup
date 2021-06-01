
data = data.frame()
nReps = c("Rep1")
nSires =c("1000Sires")
PopSize = c("8000")
HeritHerdRatio = c("0.2_0.8")
HerdSize = c("1","2","4","8","16","32")

for (rep in nReps){
  for (connectedness in nSires){
    for (pop in PopSize){
      for (ratio in HeritHerdRatio){
        for (size in HerdSize){
          filepath = paste("./simulationData",rep,connectedness,pop,ratio,size,sep="/")
	  print(filepath); setwd(filepath)
	  commandArgs <- function(...) pop
	  source("./resultCollectScripts/Acc_Comparison_Random_Wombat_GBLUP.R")
	  new_df = data.frame(matrix(NA, nrow = 1, ncol = 1)) ; new_df$Rep = rep ; new_df$nSires = connectedness ; new_df$PopSize = pop ; new_df$HHRatio = ratio ; new_df$HerdSize = size
	  new_df$EBV_Acc = r_acc; new_df$EBV_Slope = rbv_bias; new_df$Intercept = rbv_int; new_df$UQ_Acc = UQ_rcor ; new_df$LQ_Acc = LQ_rcor
	  new_df$Model = "Herd Random" ; data = rbind(data,new_df) ; new_df = data.frame()
          commandArgs <- function(...) pop
	  source("./resultCollectScripts/Acc_Comparison_Fixed_Wombat_GBLUP.R")
	  new_df = data.frame(matrix(NA, nrow = 1, ncol = 1)) ; new_df$Rep = rep ; new_df$nSires = connectedness ; new_df$PopSize = pop ; new_df$HHRatio = ratio ; new_df$HerdSize = size
          new_df$EBV_Acc = f_acc; new_df$EBV_Slope = fbv_bias; new_df$Intercept = fbv_int; new_df$UQ_Acc = UQ_fcor ; new_df$LQ_Acc = LQ_fcor
          new_df$Model = "Herd Fixed" ; data = rbind(data,new_df) ; new_df = data.frame()
          commandArgs <- function(...) pop
	  source("./resultCollectScripts/Acc_Comparison_NoHerd_Wombat_GBLUP.R")
          new_df = data.frame(matrix(NA, nrow = 1, ncol = 1)) ; new_df$Rep = rep ; new_df$nSires = connectedness ; new_df$PopSize = pop ; new_df$HHRatio = ratio ; new_df$HerdSize = size
          new_df$EBV_Acc = acc; new_df$EBV_Slope = bv_bias; new_df$Intercept = bv_int; new_df$UQ_Acc = UQ_cor ; new_df$LQ_Acc = LQ_cor
          new_df$Model = "No Herd" ; data = rbind(data,new_df) ; new_df = data.frame()
        }
      }
    }
  }
}
saveRDS(data,file="WombatGBLUP.rds")

