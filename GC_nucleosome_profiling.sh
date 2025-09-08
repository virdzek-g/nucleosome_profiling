#! /bin/bash
#$ -e /Griffin/err
#$ -o /Griffin/out
#$ -pe smp 15 -R y -binding linear:15
#$ -l h_vmem=20G
#$ -l h_rt=72:00:00
#$ -N GV_GC_nucleosome_profiling

source /broad/software/scripts/useuse
reuse Anaconda3
#reuse .anaconda3-2022.05

#source activate /broad/hptmp/Virdzekova/Griffin_conda_env
source activate griffin_demo

cd /griffin_nucleosome_profiling
snakemake -s griffin_nucleosome_profiling.snakefile --cores 1
