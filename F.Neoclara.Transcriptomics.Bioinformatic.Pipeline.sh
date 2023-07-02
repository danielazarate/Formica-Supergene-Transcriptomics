# TRANSCRIPTOMICS of FORMICA NEOCLARA 
# Written by Daniela Zarate
________________________________________________________________________________________________________
# PROJECT GOAL: To determine if differential genotypes of Formica Neoclara show differential \
# RNA gene expression. 
# Dataset composed of 40 workers. 
________________________________________________________________________________________________________
# Important short cuts I should keep handy:
# branch to a new partition and work on an interactive command line 
srun -p short --pty bash -l 
squeue -u danielaz # check on status of jobs 
________________________________________________________________________________________________________
# full path to reference

REFERENCE=/bigdata/brelsfordlab/abrelsford/form_wgs/dovetail/glacialis/glac.v0.1.fa.gz
REFERENCE=/rhome/danielaz/bigdata/transcriptomics/glacialisREF/glac.v0.1.fa.gz

# Find the original raw data here: 
/rhome/danielaz/bigdata/transcriptomics/rawfastq

# Current personal working directory:
/rhome/danielaz/bigdata/transcriptomics
________________________________________________________________________________________________________

# Create a list of all the individuals, using only R1 reads:
ls neoclara_*_*_*_*_R1_001.fastq.gz > neoclara.R1.samples.txt
# Split the names of the files in order to have basenames: 
awk '{split($0,a,"_L"); print a[1]}' neoclara.R1.samples.txt > neoclara.samples.txt

__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# Acquire TruSeq3-PE.fa fasta file for adaptor removal:
https://github.com/timflutre/trimmomatic/blob/master/adapters/TruSeq3-PE.fa

# NOTE: Although it is uncertain whether this is optimal for RNA-seq
# NOTE: Although trimming itself seems like it might not be super important e.g. Liao & Shi (2020) 

vi TrueSeq3-PE.fa

>PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT

vi trimmomatic.sh 
chmod +x trimmomatic.sh 

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
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
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

# use loopsub to submit on command line 
while read i ; do sed "s/PLACEHOLDER/$i/g" trimmomatic.sh > trimmomatic.$i.sh; sbatch trimmomatic.$i.sh ; done<neoclara.samples.txt

# On average 98% retained reads 
__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# Create STAR index of reference genome:

vi starIndex.sh
chmod +x starIndex.sh

squeue -u danielaz 

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
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


~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

# submit script. 
sbatch starIndex.sh

__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# Create STAR alignment script

vi starAlign.sh
chmod +x starAlign.sh

neoclara_60_6_S40.forward.paired
neoclara_60_6_S40.foward.unpaired
neoclara_60_6_S40.reverse.paired
neoclara_60_6_S40.reverse.unpaired

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
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


~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

# use loopsub to submit on command line 
while read i ; do sed "s/PLACEHOLDER/$i/g" starAlign.sh > starAlign.$i.sh; sbatch starAlign.$i.sh ; done<neoclara.samples.txt
__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.







squeue -u danielaz 