#!/bin/sh

declare -a nReps=(Rep1 Rep2 Rep3 Rep4 Rep5 Rep6 Rep7 Rep8 Rep9 Rep10)
declare -a nSires=(1000)
declare -a nInd=(8000)
declare -a herits=(0.2_0.8 1_1.33 5_4)
declare -a nIndPerHerd=(1 2 4 8 16 32)
declare -a model=(Random) 
declare -a method=(PBLUP)

for reps in "${nReps[@]}"
do
	cd $reps	
	for sires in "${nSires[@]}"
	do
		cd $sires"Sires"
		for popsize in "${nInd[@]}"
		do
			for herit in "${herits[@]}"
			do
				for herdsize in "${nIndPerHerd[@]}"
				do
					cd ./simulationData/$reps/$sires"Sires"/$popsize/$herit/$herdsize
				    rm -r $method$model"Model"
					mkdir $method$model"Model"
					cd $method$model"Model"
					pwd
					cp ./../../../../../../../geneticEvalutionScripts/Pedigree/"wombat"$model".par" ./wombat.par
					var=(${herit//_/ })
					sed -i "s/GENVAR/${var[0]}/g" ./wombat.par
                    sed -i "s/HERDVAR/${var[1]}/g" ./wombat.par
					#Rscript ./../../../../../../../geneticEvalutionScripts/Pedigree/PblupPedigree.R
					awk '{if (NR!=1) {print}}' ./../PedigreeAndScaledValues.txt > ./Data.txt
					tail -n +2 ../../../../GlobalPedigree.txt > ./Ped.txt
					./../../../../../../../wombat --setup --noprune
					cp ./../../../../../../../geneticEvalutionScripts/Pedigree/"wombat"$model$method".par" ./wombat.par
					var=(${herit//_/ })
                    sed -i "s/GENVAR/${var[0]}/g" ./wombat.par
                    sed -i "s/HERDVAR/${var[1]}/g" ./wombat.par
					./../../../../../../../wombat --blup 
					cd ./../../../../../../../simulationData/$reps
				done
			done
		done
	done
done
