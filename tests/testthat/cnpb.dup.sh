#!/bin/bash
#$ -cwd
#$ -j y
#$ -l mem_free=2G
#$ -l h_vmem=3G
#$ -pe local 1
#$ -t 1-120
rlib="/users/mchou/R/x86_64-pc-linux-gnu-library/3.5"
WD=$PWD
cd $WD
R_LIBS_USER=$rlib Rscript ./duplication.script.R $SGE_TASK_ID 

