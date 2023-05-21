#!/bin/bash -l
#SBATCH --time=90:00:00
#SBATCH --ntasks=8
#SBATCH --mem=245g
#SBATCH --tmp=200g
#SBATCH --job-name=rn21a

cd ~/git/nf/rn21a/rnaseq

conda activate nf

nextflow run $NXF_HOME/rnaseq -params-file genomes.yml -profile mangi