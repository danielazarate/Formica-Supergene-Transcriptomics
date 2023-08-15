#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30 # max cores per user on highmem is 32 
#SBATCH --mem-per-cpu=3G # 1 GB per 1M reads (for a total of 90 GB over 30 CPUs)
#SBATCH --time=0-05:00:00 # 3 days, 1 HR per 1M reads
#SBATCH --output=DESeq2.stdout
#SBATCH --mail-user=danielaz@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="DESeq2.job"
#SBATCH -p intel # Available paritions: intel, batch, highmem, gpu, short (each has walltime and memory limits)

# Print current date
date

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd /rhome/danielaz/bigdata/transcriptomics/raw_fastq

module load trinity-rnaseq/2.14.0

$TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl \
      --matrix Trinity.isoform.counts.matrix \
      --samples_file neoclara.samples.v2.txt  \
      --method DESeq2 \
      --output DESeq2_trans

# print name of node
hostname

# runtime: <45 min for 8 PW comparisons 
# runtime <2 hours for many more PW comparisons 
