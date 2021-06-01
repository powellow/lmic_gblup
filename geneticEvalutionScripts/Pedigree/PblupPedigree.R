#GP = read.table("GlobalPedigree.txt",header = T)
damped = read.table("../PedigreeAndScaledValues.txt",header = T)
damped = damped[order(damped$ID),]
trimped = damped[,c(1:3)] ; names(trimped)=c("id","mother","father")
#ped = rbind(GP,trimped); rownames(ped) <- NULL

#write.table(ped,"Pedigree.txt",row.names=F,quote = F)
write.table(damped,"PedigreeASV.txt",row.names=F,quote = F)

