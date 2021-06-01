#!/bin/sh
#$ -cwd
#$ -N Sim.Rep1
#$ -l h_vmem=32G
#$ -l h_rt=06:00:0

. /etc/profile.d/modules.sh

module load intel/2017u4

source ~/.bashrc

sh SimsUpdate.sh Rep1
