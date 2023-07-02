#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-08:00:00 
#SBATCH --output=PLACEHOLDER.stdout
#SBATCH --mail-user=danielaz@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="PLACEHOLDER-trimmomatic-log"
#SBATCH -p intel # Available paritions: intel, batch, highmem, gpu, short (each has walltime and memory limits)

# Print current date
date

# Trimmomatic is a fast, multithreaded command line tool that can be used to trim and crop \
# Illumina (FASTQ) data as well as to remove adapters. These adapters can pose a real problem \
# depending on the library preparation and downstream application. 
# Load software (version 0.39)
module load trimmomatic

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd $SLURM_SUBMIT_DIR

READ1=/rhome/danielaz/bigdata/transcriptomics/raw_fastq/PLACEHOLDER_R1_001.fastq.gz
READ2=/rhome/danielaz/bigdata/transcriptomics/raw_fastq/PLACEHOLDER_R2_001.fastq.gz
OUTPUT1=/rhome/danielaz/bigdata/transcriptomics/trim_fastq


trimmomatic PE ${READ1} ${READ2} \
 ${OUTPUT1}/PLACEHOLDER.forward.paired \
 ${OUTPUT2}/PLACEHOLDER.foward.unpaired \
 ${OUTPUT2}/PLACEHOLDER.reverse.paired \
 ${OUTPUT2}/PLACEHOLDER.reverse.unpaired \
 ILLUMINACLIP:TrueSeq3-PE.fa:2:30:10 \
 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Print name of node
hostname

# NexteraPE-PE is the fasta file’s name that contain the adapters sequence \
#(given with the program; you could also add your custom ones). You may have t \
# specify the path to it in certain conditions. Beware, Nextera adapters (works \
# for Nextera XT too) are always PE adapters (can be used for PE and SE).
# :2:30:10 are mismatch/accuracy treshold for adapter/reads pairing.
# LEADING:3 is the quality under which leading (hence the first, at the beginning of the read) nucleotide is trimmed.
# TRAILING:3 is the quality under which trailing (hence the last, at the end of the read) nucleotide is trimmed.
# SLIDINGWINDOW:4:15 Trimmomatic scans the reads in a 4 base window… If mean quality drops under 15, the read is trimmed.
# MINLEN:32 is the minimum length of trimmed/controled reads (here 32). If the read is smaller, it is discarded.
