#!/bin/sh

declare -a nReps=(Rep1)
declare -a nSires=(1000)
declare -a nInd=(8000)
declare -a herits=(0.2_0.8)
declare -a nIndPerHerd=(1 2 4 8 16 32)
declare -a model=(Random) 
declare -a method=(HBLUP)

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
					cp ../../../DamGenotypes.txt ./Genotypes.txt
					
					tail -n +2 ./simulationData/$reps/$sires"Sires"/GlobalPedigree.txt > ./Ped.txt
					#Rscript ./geneticEvalutionScripts/Pedigree/PblupPedigree.R
					awk '{if (NR!=1) {print}}' ./../PedigreeAndScaledValues.txt > Data.txt
					awk '{print $1, 0, 0, 2}' Data.txt > TrimPed.txt
					cp ./../../../../../../../geneticEvalutionScripts/Genomics/"HinvMeta" ./
					cp ./../../../../../../../geneticEvalutionScripts/Genomics/"ginv.par" ./wombat.par
					./../../../../../../../wombat
					cp ./../../../../../../../geneticEvalutionScripts/Genomics/"wombat"$model$method".par" ./wombat.par
					var=(${herit//_/ })
                    sed -i "s/GENVAR/${var[0]}/g" ./wombat.par
                    sed -i "s/HERDVAR/${var[1]}/g" ./wombat.par
					./../../../../../../../wombat
					cd ./../../../../../../../simulationData/$reps
				done
			done
		done
	done
done
