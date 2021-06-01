library(tidyr)

load('BreedingPopulation.RData')

ped <- SP$pedigree
id <- c(1:nrow(ped))
formatted_ped <- cbind(id,ped[,1:2])
formatted_ped <- formatted_ped %>% tail(.,n=480000) %>% head(.,n=400000) %>% data.frame()
names(formatted_ped) <- c("id","mother","father")
write.table(formatted_ped,"GlobalPedigree.txt",row.names=F,quote=F)
