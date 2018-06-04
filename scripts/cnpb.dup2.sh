#!/bin/bash
#$ -cwd
#$ -j y
#$ -l mem_free=2G
#$ -l h_vmem=3G
#$ -pe local 1
#$ -t 1-24
rlib="/users/rscharpf/bin/R"
WD=$PWD
cd $WD
R_LIBS_USER=$rlib Rscript ./cnpbayes.trio.comparison.batch.dup.R $SGE_TASK_ID 

