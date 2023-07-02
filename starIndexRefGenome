#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=0-08:00:00 
#SBATCH --output=starIndex.stdout
#SBATCH --mail-user=danielaz@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="starIndex-log"
#SBATCH -p intel # Available paritions: intel, batch, highmem, gpu, short (each has walltime and memory limits)

# Print current date
date

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd $SLURM_SUBMIT_DIR

# Spliced Transcripts Alignment to a Reference © Alexander Dobin, 2009-2022
# STAR is an aligner designed to specifically address many of the challenges of RNA-seq data mapping \
# using a strategy to account for spliced alignments. STAR is shown to have high accuracy and outperforms \
#  other aligners by more than a factor of 50 in mapping speed, but it is memory intensive. 
# Load software, version: 2.7.10b
module load star

# STEP ONE - CREATE GENOME INDEX 

REFERENCE=/rhome/danielaz/bigdata/transcriptomics/glacialisREF/glac.v0.1.fa
STAR_DIR=/rhome/danielaz/bigdata/transcriptomics/starIndex

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir ${STAR_DIR} \
--genomeFastaFiles ${REFERENCE} \
--sjdbOverhang 99

# --genomeDir : Directory to store output 
# NOTE: In case of reads of varying length, the ideal value for --sjdbOverhang is max(ReadLength)-1. \
# In most cases, the default value of 100 will work similarly to the ideal value
# When generating genome without annotations, do not specify : --sjdbOverhang 

# STAR uses the genome file (FASTA) and gene annotation file (GTF or GFF3) to create the genome indices. \
# The gene annotation file is needed for creating the known splice junctions to improve the accuracy of \
# the genome mapping. Gene annotation file is optional, but it is highly recommended if it is available 
# (Note: you can also provide an annotation file during the mapping step).

# More on STAR:
# The algorithm achieves this highly efficient mapping by performing a two-step process: 1. Seed searching and 2. Clustering, \
# stitching, and sorting. For every read that STAR aligns, STAR will search for the longest sequence that exactly \
# matches one or more locations on the reference genome. These longest matching sequences are called the Maximal \
# Mappable Prefixes (MMPs): The different parts of the read that are mapped separately are called ‘seeds’. So the \
# first MMP that is mapped to the genome is called seed1. STAR will then search again for only the unmapped portion \
# of the read to find the next longest sequence that exactly matches the reference genome, or the next MMP, which will be seed2. \
# This sequential searching of only the unmapped portions of reads underlies the efficiency of the STAR algorithm. STAR uses \
# an uncompressed suffix array (SA) to efficiently search for the MMPs, this allows for quick searching against even the largest \
# reference genomes. Other slower aligners use algorithms that often search for the entire read sequence before splitting \
# reads and performing iterative rounds of mapping. If STAR does not find an exact matching sequence for each part of the read \
# due to mismatches or indels, the previous MMPs will be extended. If extension does not give a good alignment, then the poor \
# quality or adapter sequence (or other contaminating sequence) will be soft clipped. The separate seeds are stitched together \
# to create a complete read by first clustering the seeds together based on proximity to a set of ‘anchor’ seeds, or seeds that are not multi-mapping.
# Then the seeds are stitched together based on the best alignment for the read (scoring based on mismatches, indels, gaps, etc.).
