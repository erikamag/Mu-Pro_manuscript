#!/bin/bash -l
#SBATCH --time=80:00:00
#SBATCH --ntasks=6
#SBATCH --mem=150g
#SBATCH --tmp=24g
#SBATCH --job-name=rn21a.trinity.job.array
#SBATCH --error=rn21a.trinity_%A_%a.err
#SBATCH --output=rn21a.trinity_%A_%a.out
#SBATCH --array=1-32
echo "$SLURM_ARRAY_TASK_ID"

module load trinity
module load trimmomatic/0.33

cd /scratch.global/rn21a_fq/

file=$(ls *.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

Trinity --seqType fq --max_memory 150G --samples_file $file --CPU 6 --trimmomatic --min_contig_length 200 --output $file.trinity.out
