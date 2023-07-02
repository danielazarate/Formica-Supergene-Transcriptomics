#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=0-08:00:00 
#SBATCH --output=starAlign.stdout
#SBATCH --mail-user=danielaz@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="starAlign-log"
#SBATCH -p intel # Available paritions: intel, batch, highmem, gpu, short (each has walltime and memory limits)

# Print current date
date

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd $SLURM_SUBMIT_DIR

# Spliced Transcripts Alignment to a Reference © Alexander Dobin, 2009-2022
# STAR is an aligner designed to specifically address many of the challenges of RNA-seq data mapping \
# using a strategy to account for spliced alignments. STAR is shown to have high accuracy and outperforms \
# other aligners by more than a factor of 50 in mapping speed, but it is memory intensive. 
# Load software, version: 2.7.10b
module load star

STAR_INDEX=/rhome/danielaz/bigdata/transcriptomics/starIndex
DIR=/rhome/danielaz/bigdata/transcriptomics/raw_fastq

STAR --runThreadN 12 \
--readFilesIn PLACEHOLDER.forward.paired,PLACEHOLDER.foward.unpaired PLACEHOLDER.reverse.paired,PLACEHOLDER.reverse.paired \
--genomeDir ${STAR_INDEX} \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix PLACEHOLDER.map \
--outSAMunmapped Within

# --runThreadN : Number of threads (processors) for mapping reads to genome
# --readFilesIn : Read files for mapping to the genome.
# --genomeDir : PATH to the directory containing built genome indices
# --outSAMtype : Output coordinate sorted BAM file which is useful for many downstream analyses. This is optional.
# --outSAMunmapped : Output unmapped reads from the main SAM file in SAM format. This is optional
# --outFileNamePrefix : Provide output file prefix name
