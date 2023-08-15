#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30 # max cores per user on highmem is 32 
#SBATCH --mem-per-cpu=3G # 1 GB per 1M reads (for a total of 90 GB over 30 CPUs)
#SBATCH --time=0-20:00:00 # 3 days, 1 HR per 1M reads
#SBATCH --output=trinity.matrix.stdout
#SBATCH --mail-user=danielaz@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="transcript.matrix.job"
#SBATCH -p intel # Available paritions: intel, batch, highmem, gpu, short (each has walltime and memory limits)

# Print current date
date

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd /rhome/danielaz/bigdata/transcriptomics/raw_fastq

module load trinity-rnaseq/2.14.0

abundance_estimates_to_matrix.pl --est_method salmon \
--out_prefix Trinity --name_sample_by_basedir \
--quant_files quant_files.v1.list \
--gene_trans_map trinity_out_dir/Trinity.fasta.gene_trans_map

# print name of node
hostname

# Will output 11 files in the directory: 
# Trinity.v1.gene.counts.matrix
# Trinity.v1.gene.TPM.not_cross_norm           
# Trinity.v1.gene.TPM.not_cross_norm.TMM_info.txt  
# Trinity.v1.isoform.TMM.EXPR.matrix     
# Trinity.v1.isoform.TPM.not_cross_norm.runTMM.R
# Trinity.v1.gene.TMM.EXPR.matrix 
# Trinity.v1.gene.TPM.not_cross_norm.runTMM.R  
# Trinity.v1.isoform.counts.matrix                 
# Trinity.v1.isoform.TPM.not_cross_norm  
# Trinity.v1.isoform.TPM.not_cross_norm.TMM_info.txt

# runtime: <5 min
