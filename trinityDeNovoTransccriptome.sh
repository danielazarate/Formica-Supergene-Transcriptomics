#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32 # max cores per user on highmem
#SBATCH --mem-per-cpu=78G # 1 GB per 1M reads
#SBATCH --time=3-00:00:00 # 3 days, 1 HR per 1M reads
#SBATCH --output=trinity.stdout
#SBATCH --mail-user=danielaz@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="trinity1.0"
#SBATCH -p highmem # Available paritions: intel, batch, highmem, gpu, short (each has walltime and memory limits)

# Print current date
date

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd $SLURM_SUBMIT_DIR

module load trinity-rnaseq/2.14.0

RAW_FASTQ=/rhome/danielaz/bigdata/transcriptomics/raw_fastq/
SAMPLE_LIST=/rhome/danielaz/bigdata/transcriptomics/trinity/neoclara.deNovo.samples 

time Trinity  --seqType fq  --samples_file ${SAMPLE_LIST} \
    --min_contig_length 150 --CPU 32 --output deNovo_Neoclara_trinity


# --seqType : data input type, either fasta or fastq
# --min_contig_length : minimum assembled contig length to report (default=200)
# --CPU : number of CPUs to use (default=2)
# --output : name of directory to output files in 


# print name of node
hostname
