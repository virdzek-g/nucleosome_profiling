#! /bin/bash
#$ -e ./Griffin/err
#$ -o ./Griffin/out
#$ -pe smp 15 -R y -binding linear:15
#$ -l h_vmem=20G
#$ -l h_rt=95:00:00
#$ -N GV_GC_mappability_correction

source /broad/software/scripts/useuse
reuse Anaconda3
#reuse .anaconda3-2022.05


source activate griffin_demo

cd /griffin_GC_and_mappability_correction
snakemake -s griffin_GC_and_mappability_correction.snakefile --cores 1
