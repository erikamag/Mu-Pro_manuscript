#!/bin/bash -l
#SBATCH --time=80:00:00
#SBATCH --ntasks=6
#SBATCH --mem=150g
#SBATCH --tmp=24g
#SBATCH --job-name=rn21a.trinity.e2f13m1
#SBATCH --error=rn21a.trinity.e2f13m1.err
#SBATCH --output=rn21a.trinity.e2f13m1.out
    
module load trinity
module load trimmomatic/0.33

cd /scratch.global/rn21a/rn21a_fq/

Trinity --seqType fq --max_memory 150G --samples_file rn21a.samples.e2f13m1.txt --CPU 6 --trimmomatic --min_contig_length 200 --output C14.e2f13m1.trinity.out
