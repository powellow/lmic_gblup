#!/bin/sh

declare -a nReps=($1) 
declare -a nSires=(1000)
declare -a nInd=(8000)
declare -a herits=(0.2_0.8) 

mkdir simulationData
cd ./simulationData
for reps in "${nReps[@]}"
do
	mkdir $reps
	cd $reps
	Rscript ./../../simulationScripts/CreateParents.R
	pwd
	for sires in "${nSires[@]}"
	do	
		mkdir $sires"Sires"
		Rscript ./../../simulationScripts/SimPops.R $sires
		echo $sires
		cd $sires"Sires"
		pwd
		Rscript ./../../../simulationScripts/GlobalPedigree.R
		for popsize in "${nInd[@]}"
		do
			mkdir $popsize
			awk '{if (NR!=1) {print}}' ./Dams/info.txt > ./Dams/info_trim.txt
			seq 40000 | shuf -n $popsize | awk 'NR==FNR{ z[$0]++;next}
			{if (FNR in z){ print >FILENAME"_subsetted"}}' - ./Dams/pheno.txt ./Dams/genotype.txt ./Dams/gv.txt ./Dams/info_trim.txt
			mv ./Dams/info_trim.txt_subsetted ./Dams/info.txt_subsetted ; rm ./Dams/info_trim.txt
			awk '{print $1}' ./Dams/info.txt_subsetted > ./Dams/id.txt_subsetted
			paste -d ' ' ./Dams/id.txt_subsetted ./Dams/genotype.txt_subsetted > ./$popsize/DamGenotypes.txt
 			mv ./Dams/*_subsetted ./$popsize/
			echo $popsize; cd ./$popsize
			pwd
			for herit in "${herits[@]}"
			do
				echo $herit
                mkdir $herit;  
				cd ./$herit ;  mkdir 1 2 4 8 16 32
                cp ../*_subsetted  ./1; cp ../*_subsetted  ./2; cp ../*_subsetted  ./4; cp ../*_subsetted  ./8; cp ../*_subsetted  ./16; cp ../*_subsetted  ./32
				Rscript ./../../../../../simulationScripts/Scale_Phenotypes.R $herit $sires $popsize
 				pwd
				cd ../
			done
			cd ../
		done
		cd ../
	done
	cd ../	
done
