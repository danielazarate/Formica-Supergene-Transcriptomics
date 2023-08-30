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
 ${OUTPUT2}/PLACEHOLDER.forward.unpaired \
 ${OUTPUT2}/PLACEHOLDER.reverse.paired \
 ${OUTPUT2}/PLACEHOLDER.reverse.unpaired \
 ILLUMINACLIP:TrueSeq3-PE.fa:2:30:10 \
 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Print name of node
hostname

# NexteraPE-PE is the fasta file’s name that contain the adapters sequence \
#(given with the program; you could also add your custom ones). You may have to \
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

# example file outputs: 
neoclara_60_6_S40.forward.paired
neoclara_60_6_S40.foward.unpaired
neoclara_60_6_S40.reverse.paired
neoclara_60_6_S40.reverse.unpaired
__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# Now we will generate a de novo transcriptome with the neoclara data using trinity

# Create a sample ID text file:

vi neoclara.deNovo.samples 

colony_29       sm_sp   neoclara_29_5_S4.forward.paired         neoclara_29_5_S4.reverse.paired
colony_58       sp_sp   neoclara_58_6_S32.forward.paired        neoclara_58_6_S32.reverse.paired
colony_60       sm_sm   neoclara_60_2_S38.forward.paired        neoclara_60_2_S38.reverse.paired
__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

vi trinity.sh
chmod +x trinity.sh

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30 # max cores per user on highmem is 32 
#SBATCH --mem-per-cpu=3G # 1 GB per 1M reads (for a total of 90 GB over 30 CPUs)
#SBATCH --time=5-00:00:00 # 3 days, 1 HR per 1M reads
#SBATCH --output=trinity.stdout
#SBATCH --mail-user=danielaz@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="trinity1.0"
#SBATCH -p intel # Available paritions: intel, batch, highmem, gpu, short (each has walltime and memory limits)

# Print current date
date

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd $SLURM_SUBMIT_DIR

module load trinity-rnaseq/2.14.0

SAMPLE_LIST=/rhome/danielaz/bigdata/transcriptomics/raw_fastq/neoclara.deNovo.samples 

time Trinity  --seqType fq  --samples_file ${SAMPLE_LIST} \
    --min_contig_length 150 --CPU 30 --max_memory 90G \
    --output deNovo_Neoclara_trinity


# --seqType : data input type, either fasta or fastq
# --min_contig_length : minimum assembled contig length to report (default=200)
# --CPU : number of CPUs to use (default=2)
# --output : name of directory to output files in 
# --max_memory : suggested max memory to use by Trinity where limiting can be enabled 


# print name of node
hostname

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

vi transcript.quant.sh
chmod +x transcript.quant.sh

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30 # max cores per user on highmem is 32 
#SBATCH --mem-per-cpu=3G # 1 GB per 1M reads (for a total of 90 GB over 30 CPUs)
#SBATCH --time=0-20:00:00 # 3 days, 1 HR per 1M reads
#SBATCH --output=transcript.quant.stdout
#SBATCH --mail-user=danielaz@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="transcript.quant.job"
#SBATCH -p intel # Available paritions: intel, batch, highmem, gpu, short (each has walltime and memory limits)

# Print current date
date

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd /rhome/danielaz/bigdata/transcriptomics/raw_fastq

module load trinity-rnaseq/2.14.0

FASTQ=/rhome/danielaz/bigdata/transcriptomics/raw_fastq
SAMPLES=/rhome/danielaz/bigdata/transcriptomics/trinity/neoclara.trinity.IDs.txt
DE_NOVO=/rhome/danielaz/bigdata/transcriptomics/trinity/neoclara.deNovo.fa


align_and_estimate_abundance.pl --seqType fq  \
    --samples_file ${SAMPLES}  --transcripts ${DE_NOVO} \
    --est_method salmon  --trinity_mode   --prep_reference  

# --seqType : data input type, either fasta or fastq
# --samples file : list of samples
# --transcripts : de novo transcriptome
# --output-dir : name of directory to output files in (default to name of sample)
# --est_method : abundance estimation method (alignment free - Salmon)
# --trinity_mode : will automatically generated the gene_trans_map
# --prep_reference : prep reference (builds target index)

# NOTE: the sample file's second column (replicate) will be the name of the output directory.
# OUTPUTS: A directory per individual provided in the sample file with these files: 
    # aux_info dir
    # cmd_info.json
    # lib_format_counts.json
    # libParams dir
    # logs dir 
    # quant.sf
    # quant.sf.genes 

    # The quant.sf file is the file of interest and has the following information:
        # NumReads’ corresponding to the number of RNA-Seq fragments predicted to be derived from that transcript
        # TPM: transcripts per million 

# print name of node
hostname

# run time < 3 hours 

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

sbatch transcript.quant.sh 

__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# Now, given the expression estimates for each of the transcripts in each of the samples, \
# pull together all values into matrices containing transcript IDs in the rows, and sample names in the columns.

# create a list of all the quant.sf files from each directory
find c*.rep.* -name "quant.sf" | tee quant_files.v2.list

__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# Use Trinity to generate the count and expression matrices for both the transcript isoforms and sepearate files for ‘gene’s.
# produce two matrices, one containing the estimated counts, and another containing the TPM expression values \
# that are cross-sample normalized using the TMM method. 

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

vi trinity.matrix.sh
chmod +x trinity.matrix.sh

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

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

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

sbatch trinity.matrix.sh 
__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# Files of interest:
# Trinity.isoform.counts.matrix’, which contains the counts of RNA-Seq fragments mapped to each transcript.
# Trinity.isoform.TMM.EXPR.matrix, which contains normalized expression matrix

Trinity.isoform.counts.matrix 
Trinity.isoform.TMM.EXPR.matrix
__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# Differential Expression Using DESeq2 of Isoforms 

# NOTE ON ISOFORMS: This data allows us to study expression both at the gene and isoform level, it's often \
# useful to study both, particularly in cases where differential transcript usage exists (isoform switching) \
# where differences in expression may not be apparent at the gene level. 

vi DESeq2.sh
chmod +x DESeq2.sh

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

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
      --matrix Trinity.v2.isoform.counts.matrix \
      --samples_file neoclara.samples.v2.txt  \
      --method DESeq2 \
      --output DESeq2_trans

# print name of node
hostname

# runtime: <45 min for 8 PW comparisons 
# runtime <2 hours for many more PW comparisons 

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

sbatch DESeq2.sh 

# Examine the contents of the DESeq2_trans/ directory.
# The _"results" file is the file of interest which contains the log fold change (log2FoldChange), \
#  mean expression (baseMean), P- value from an exact test, and false discovery rate (padj)
# It also produces MA and Volcano plots 
# MA plot: MA plot compares the log fold change against the mean of the normalized counts. 
# Volcano plot: The colored points are differentially expressed genes with

cd DESeq2_trans_v2  

# download volcano plots 
scp 'danielaz@cluster.hpcc.ucr.edu:/rhome/danielaz/bigdata/transcriptomics/raw_fastq/DESeq2_trans_v2/*.pdf' . 

# download matrix 

scp danielaz@cluster.hpcc.ucr.edu:/rhome/danielaz/bigdata/transcriptomics/raw_fastq/Trinity.v2.isoform.counts.matrix ./isoform.matrix

scp danielaz@cluster.hpcc.ucr.edu:/rhome/danielaz/bigdata/transcriptomics/raw_fastq/DESeq2_trans_v2/diffExpr.P1e-3_C2.matrix .

__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.


# Extract DE transcripts that are at least 4-fold differentially expressed at a significane of <=0.001 in any \
# of the pairwise comparisons. 

vi DESeq2_DE.sh
chmod +x DESeq2_DE.sh

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30 # max cores per user on highmem is 32 
#SBATCH --mem-per-cpu=3G # 1 GB per 1M reads (for a total of 90 GB over 30 CPUs)
#SBATCH --time=0-05:00:00 # 3 days, 1 HR per 1M reads
#SBATCH --output=DESeq_DE.stdout
#SBATCH --mail-user=danielaz@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="DESeq2_DE.job"
#SBATCH -p intel # Available paritions: intel, batch, highmem, gpu, short (each has walltime and memory limits)

# Print current date
date

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd /rhome/danielaz/bigdata/transcriptomics/raw_fastq/DESeq2_trans_v2  

module load trinity-rnaseq/2.14.0

$TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl \
      --matrix ../Trinity.v2.isoform.TMM.EXPR.matrix \
      --samples ../neoclara.samples.v2.txt  \
      -P 1e-3 -C 2 

# -P : False Discovery Rate threshold 
# -C : Fold change, where C is set to 2^(#) so when C is set to 2, it is 2^2 = 4 
# This will generate several files with diffExpr.P1e-3_C2’ prefix 


# print name of node
hostname

# runtime = < 5 minutes 

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

sbatch DESeq2_DE.sh 

# Observe the diffExpr.P1e-3_C2.matrix file, this is the subset of the FPKM matrix corresponding to the \
# DE transcripts identified at this threshold. 

wc -l diffExpr.P1e-3_C2.matrix

# DE transcripts observed: 793

# The file 'diffExpr.P1e-3_C2.matrix.log2.centered.genes_vs_samples_heatmap.pdf' is also produced which \
# shows transcripts clustered along the vertical axis and samples clustered along the horizontal axis. 
# Expression values are plotted in log2 space and mean centered (mean expression value for each feature is subtracted \
# from each of its expression values in that row) and upregulation is yellow and downregulation is purple. 

# download heatmaps
scp 'danielaz@cluster.hpcc.ucr.edu:/rhome/danielaz/bigdata/transcriptomics/raw_fastq/DESeq2_trans_v2/*_heatmap.pdf' . 
__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# Extract transcript clusters by expression profile by cutting dendrogram at a given percent of its height. 

vi DESeq2.dendrogram.sh
chmod +x DESeq2.dendrogram.sh

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30 # max cores per user on highmem is 32 
#SBATCH --mem-per-cpu=3G # 1 GB per 1M reads (for a total of 90 GB over 30 CPUs)
#SBATCH --time=0-05:00:00 # 3 days, 1 HR per 1M reads
#SBATCH --output=DESeq_dendrogram.stdout
#SBATCH --mail-user=danielaz@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="DESeq2_dendrogram.job"
#SBATCH -p intel # Available paritions: intel, batch, highmem, gpu, short (each has walltime and memory limits)

# Print current date
date

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd /rhome/danielaz/bigdata/transcriptomics/raw_fastq/DESeq2_trans_v2  

module load trinity-rnaseq/2.14.0

$TRINITY_HOME/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl \
       --Ptree 60 -R diffExpr.P1e-3_C2.matrix.RData


# Creates a directory with the individual transcript clusters, including a pdf file that summarizes \
# expression values for each cluster according to individual charts. 

# -Ptree: threshold for cutting the dendrogram

# print name of node
hostname

# runtime = 
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

sbatch DESeq2.dendrogram.sh

cd  diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_60   


ERRROR !!! 
Error in is.null(Rowv) || is.na(Rowv) : 
  'length = 2' in coercion to 'logical(1)'
Calls: heatmap.3
Execution halted
Error, cmd Rscript __tmp_define_clusters.R died with ret 256 at /opt/linux/rocky/8.x/x86_64/pkgs/trinity-rnaseq/2.14.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl line 213.



# download cluster expression profiles 
scp 'danielaz@cluster.hpcc.ucr.edu:/rhome/danielaz/bigdata/transcriptomics/raw_fastq/DESeq2_trans_v2/*_plots.pdf' . 

__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.
__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# Now repeat the DE analysis with genes instead of isoforms 

vi DESeq2.genes.sh
chmod +x DESeq2.genes.sh

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30 # max cores per user on highmem is 32 
#SBATCH --mem-per-cpu=3G # 1 GB per 1M reads (for a total of 90 GB over 30 CPUs)
#SBATCH --time=0-05:00:00 # 3 days, 1 HR per 1M reads
#SBATCH --output=DESeq2.genes.stdout
#SBATCH --mail-user=danielaz@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="DESeq2.genes.job"
#SBATCH -p intel # Available paritions: intel, batch, highmem, gpu, short (each has walltime and memory limits)

# Print current date
date

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd /rhome/danielaz/bigdata/transcriptomics/raw_fastq

module load trinity-rnaseq/2.14.0

$TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl \
      --matrix Trinity.v2.gene.counts.matrix \
      --samples_file neoclara.samples.v2.txt  \
      --method DESeq2 \
      --output DESeq2_gene

# print name of node
hostname

# runtime: <45 min for 8 PW comparisons 
# runtime <2 hours for many more PW comparisons 
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

sbatch DESeq2.genes.sh

cd DESeq2_gene

# download volcano plots 
scp 'danielaz@cluster.hpcc.ucr.edu:/rhome/danielaz/bigdata/transcriptomics/raw_fastq/DESeq2_gene/*.pdf' . 


__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# Extract DE transcripts that are at least 4-fold differentially expressed at a significane of <=0.001 in any \
# of the pairwise comparisons. 

vi DESeq2_genes_DE.sh
chmod +x DESeq2_genes_DE.sh

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30 # max cores per user on highmem is 32 
#SBATCH --mem-per-cpu=3G # 1 GB per 1M reads (for a total of 90 GB over 30 CPUs)
#SBATCH --time=0-05:00:00 # 3 days, 1 HR per 1M reads
#SBATCH --output=DESeq_genes_DE.stdout
#SBATCH --mail-user=danielaz@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="DESeq2_genes_DE.job"
#SBATCH -p intel # Available paritions: intel, batch, highmem, gpu, short (each has walltime and memory limits)

# Print current date
date

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd /rhome/danielaz/bigdata/transcriptomics/raw_fastq/DESeq2_gene

module load trinity-rnaseq/2.14.0

$TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl \
      --matrix ../Trinity.v2.gene.TMM.EXPR.matrix \
      --samples ../neoclara.samples.v2.txt  \
      -P 1e-3 -C 2 

# -P : False Discovery Rate threshold 
# -C : Fold change, where C is set to 2^(#) so when C is set to 2, it is 2^2 = 4 
# This will generate several files with diffExpr.P1e-3_C2’ prefix 


# print name of node
hostname

# runtime = < 5 minutes 

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

sbatch DESeq2_genes_DE.sh 

# Observe the diffExpr.P1e-3_C2.matrix file, this is the subset of the FPKM matrix corresponding to the \
# DE transcripts identified at this threshold. 

wc -l diffExpr.P1e-3_C2.matrix 

# DE transcripts observed: 393

# The file 'diffExpr.P1e-3_C2.matrix.log2.centered.genes_vs_samples_heatmap.pdf' is also produced which \
# shows transcripts clustered along the vertical axis and samples clustered along the horizontal axis. 
# Expression values are plotted in log2 space and mean centered (mean expression value for each feature is subtracted \
# from each of its expression values in that row) and upregulation is yellow and downregulation is purple. 

ERROR !!

[1] "Reading matrix file."
Error in is.null(Rowv) || is.na(Rowv) : 
  'length = 2' in coercion to 'logical(1)'
Calls: heatmap.3
Execution halted
Error, cmd: Rscript diffExpr.P1e-3_C2.matrix.R died with ret 256 at /bigdata/operations/pkgadmin/opt/linux/centos/8.x/x86_64/pkgs/trinity-rnaseq/2.14.0/Analysis/DifferentialExpression/PtR line 1703.
Error, Error, cmd: /bigdata/operations/pkgadmin/opt/linux/centos/8.x/x86_64/pkgs/trinity-rnaseq/2.14.0/Analysis/DifferentialExpression/PtR -m diffExpr.P1e-3_C2.matrix --log2 --heatmap --min_colSums 0 --min_ro
wSums 0 --gene_dist euclidean --sample_dist euclidean --sample_cor_matrix --center_rows --save  -s ../neoclara.samples.v2.txt died with ret 6400 at /opt/linux/rocky/8.x/x86_64/pkgs/trinity-rnaseq/2.14.0/Analy
sis/DifferentialExpression/analyze_diff_expr.pl line 272.


# download heatmaps
scp 'danielaz@cluster.hpcc.ucr.edu:/rhome/danielaz/bigdata/transcriptomics/raw_fastq/DESeq2_genes/*_heatmap.pdf' . 
__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# Extract transcript clusters by expression profile by cutting dendrogram at a given percent of its height. 

vi DESeq2_genes_dendrogram.sh
chmod +x DESeq2_genes_dendrogram.sh

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30 # max cores per user on highmem is 32 
#SBATCH --mem-per-cpu=3G # 1 GB per 1M reads (for a total of 90 GB over 30 CPUs)
#SBATCH --time=0-05:00:00 # 3 days, 1 HR per 1M reads
#SBATCH --output=DESeq_genes_dendrogram.stdout
#SBATCH --mail-user=danielaz@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="DESeq2_genes_dendrogram.job"
#SBATCH -p intel # Available paritions: intel, batch, highmem, gpu, short (each has walltime and memory limits)

# Print current date
date

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd /rhome/danielaz/bigdata/transcriptomics/raw_fastq/DESeq2_genes

module load trinity-rnaseq/2.14.0

$TRINITY_HOME/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl \
       --Ptree 60 -R diffExpr.P1e-3_C2.matrix.RData


# Creates a directory with the individual transcript clusters, including a pdf file that summarizes \
# expression values for each cluster according to individual charts. 

# -Ptree: threshold for cutting the dendrogram

# print name of node
hostname

# runtime = 
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

sbatch DESeq2_genes_dendrogram.sh

cd diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_60 
 
# download cluster expression profiles 
scp 'danielaz@cluster.hpcc.ucr.edu:/rhome/danielaz/bigdata/transcriptomics/raw_fastq/DESeq2_genes/*_plots.pdf' . 

__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# Functional annotation of assembled transcripts using Trinotate 

# Trinotate is a comprehensive annotation suite designed for automatic functional annotation of transcriptomes, \
# particularly de novo assembled transcriptomes, from model or non-model organisms. Trinotate makes use of a number \
# of different well referenced methods for functional annotation including homology search to known sequence data \
# (BLAST+/SwissProt), protein domain identification (HMMER/PFAM), protein signal peptide and transmembrane domain \
# prediction (signalP/tmHMM), and leveraging various annotation databases (eggNOG/GO/Kegg databases). All functional \
# annotation data derived from the analysis of transcripts is integrated into a SQLite database which allows fast \
# efficient searching for terms with specific qualities related to a desired scientific hypothesis or a means to 
# create a whole annotation report for a transcriptome.

mkdir Trinotate 
cd Trinotate 

# Identification of likely protein-coding regions in transcripts with TransDecoder - it identifies \
# long open reading frames (ORFs) within transcripts and scores them according to their sequence composition \
# ORFs that encode sequences with compositional properties (codon frequencies) consistent with coding transcripts \
# are reported. (ORFs do not include a STOP codon).


vi trinotate.ORFs.sh
chmod +x trinotate.ORFs.sh

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30 # max cores per user on highmem is 32 
#SBATCH --mem-per-cpu=3G # 1 GB per 1M reads (for a total of 90 GB over 30 CPUs)
#SBATCH --time=0-05:00:00 # 3 days, 1 HR per 1M reads
#SBATCH --output=trinotate.ORFs.stdout
#SBATCH --mail-user=danielaz@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="trinotate.ORFs.job"
#SBATCH -p intel # Available paritions: intel, batch, highmem, gpu, short (each has walltime and memory limits)

# Print current date
date

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd /rhome/danielaz/bigdata/transcriptomics/raw_fastq/Trinotate

module load transdecoder/5.5.0 

DE_NOVO_FA=/rhome/danielaz/bigdata/transcriptomics/trinity/neoclara.deNovo.fa   

# Run TransCoder step that identifies all long ORFs 
$TRANSDECODER_HOME/TransDecoder.LongOrfs -t ${DE_NOVO_FA} --output_dir . 

# Now run the step that predicts which ORFs are likely to be coding 
$TRANSDECODER_HOME/TransDecoder.Predict -t ${DE_NOVO_FA}  --output_dir . 

# print name of node
hostname

# runtime = < 30 minutes 
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

sbatch trinotate.ORFs.sh

# This will produce several files with the 'transdecoder' prefix 
# The file of importance is:

neoclara.deNovo.fa.transdecoder.pep

# which contains the protein sequences corresponding to the predicted coding regions within the transcripts. \
# This will contain information on what "type" the protein is - whether its complete (has both a step and start codon) \
# Or whether it's missing a start (5prime_partial) or stop (3prime_partial) codon. Or both (internal). Also, the + or - \
# indicator designates what strand the coding region is on. 

# This PEP file will be used for various sequence homology and other bioinformatic analyses below. 
__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# Sequence Homology Search
# Run BLAST against SWISSPROT database to identify full-length transcripts. 
# Run BLAST with both the full length de novo transcriptome and the predicrted protein sequences \
# generated by TransDecoder (pep file).

# download the SWISSPROT database and set it up:
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
makeblastdb -in uniprot_sprot.fasta -out uniprot_sprot.fasta -input_type fasta -dbtype prot 

# note: the -in and -out have to match exactly so BLAST identifies the index file exactly in later analysis. 

vi blast.swissprot.sh
chmod +x blast.swissprot.sh

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30 # max cores per user on highmem is 32 
#SBATCH --mem-per-cpu=3G # 1 GB per 1M reads (for a total of 90 GB over 30 CPUs)
#SBATCH --time=0-05:00:00 # 3 days, 1 HR per 1M reads
#SBATCH --output=blast.swissprot.stdout
#SBATCH --mail-user=danielaz@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="blast.swissprot.job"
#SBATCH -p intel # Available paritions: intel, batch, highmem, gpu, short (each has walltime and memory limits)

# Print current date
date

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd /rhome/danielaz/bigdata/transcriptomics/blast

module load ncbi-blast/2.14.0+ 

SWISSPROT_DB=/rhome/danielaz/bigdata/transcriptomics/blast/uniprot_sprot.fasta 
DE_NOVO_FA=/rhome/danielaz/bigdata/transcriptomics/trinity/neoclara.deNovo.fa  

blastx -db ${SWISSPROT_DB} \
         -query ${DE_NOVO_FA} -num_threads 2 \
         -max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
          > swissprot.blastx.outfmt6

# print name of node
hostname

# runtime = 

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

sbatch blast.swissprot.sh


# Now, let’s look for sequence homologies by just searching our predicted protein sequences rather \
# than using the entire transcript as a target:

vi blast.pep.swissprot.sh
chmod +x blast.pep.swissprot.sh

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30 # max cores per user on highmem is 32 
#SBATCH --mem-per-cpu=3G # 1 GB per 1M reads (for a total of 90 GB over 30 CPUs)
#SBATCH --time=0-05:00:00 # 3 days, 1 HR per 1M reads
#SBATCH --output=blast.pep.swisprot.stdout
#SBATCH --mail-user=danielaz@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="blast.pep.swissprot.job"
#SBATCH -p intel # Available paritions: intel, batch, highmem, gpu, short (each has walltime and memory limits)

# Print current date
date

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd /rhome/danielaz/bigdata/transcriptomics/blast

module load ncbi-blast/2.14.0+ 

SWISSPROT_DB=/rhome/danielaz/bigdata/transcriptomics/blast/uniprot_sprot.fasta 
PEP=/rhome/danielaz/bigdata/transcriptomics/raw_fastq/Trinotate/neoclara.deNovo.fa.transdecoder.pep

blastp -db ${SWISSPROT_DB} \
         -query ${PEP} -num_threads 2 \
         -max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
          > swissprot.blastp.pep.outfmt6

# print name of node
hostname

# runtime = 1 day, 3 hours 

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

sbatch blast.pep.swissprot.sh


# Using our predicted protein sequences, let’s also run a HMMER search against the Pfam database, and \
# identify conserved domains that might be indicative or suggestive of function:

# HMMER is mainly used for finding the domain structure represented in protein sequence. It uses a \
# statistical model to represent the sequence, is much more sensative, and can pick up more divergent hits \
# than BLAST. 

# HMMER has a position-specific scoring system which capitalises on the fact thtat certain positions in a \
# sequence alignment tend to have biases in which residues are most likely to occur, and are likely to differ \
# in their probability of having an insertion/deletion. BLAST approaches equally penalize substitutions, \
# insertions, or deletions, regardless of where in the alignment they occur. 

vi hmmmscan.sh
chmod +x hmmscan.sh

# download the PFAM-A database and set it up:
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz 
hmmpress Pfam-A.hmm 

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30 # max cores per user on highmem is 32 
#SBATCH --mem-per-cpu=3G # 1 GB per 1M reads (for a total of 90 GB over 30 CPUs)
#SBATCH --time=0-05:00:00 # 3 days, 1 HR per 1M reads
#SBATCH --output=hmmscan.stdout
#SBATCH --mail-user=danielaz@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="hmmscan.job"
#SBATCH -p intel # Available paritions: intel, batch, highmem, gpu, short (each has walltime and memory limits)

# Print current date
date

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd /rhome/danielaz/bigdata/transcriptomics/blast

module load hmmer/3.3.2 

SWISSPROT_DB=/rhome/danielaz/bigdata/transcriptomics/blast/uniprot_sprot.fasta 
PEP=/rhome/danielaz/bigdata/transcriptomics/raw_fastq/Trinotate/neoclara.deNovo.fa.transdecoder.pep
PFAM=/rhome/danielaz/bigdata/transcriptomics/blast/Pfam-A.hmm  


hmmscan --cpu 2 --domtblout TrinotatePFAM.out \
          ${PFAM} ${PEP}

#  --domtblout <f>  : save parseable table of per-domain hits to file <f>

# print name of node
hostname

# runtime = 

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

sbatch hmmmscan.sh
__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# To determine the genes on the social supergene, align the de novo transcriptome to F. selysi

vi alignment.sh
chmod +x alignment.sh

~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30 # max cores per user on highmem is 32 
#SBATCH --mem-per-cpu=3G # 1 GB per 1M reads (for a total of 90 GB over 30 CPUs)
#SBATCH --time=0-05:00:00 # 3 days, 1 HR per 1M reads
#SBATCH --output=alignment.stdout
#SBATCH --mail-user=danielaz@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="alignment.job"
#SBATCH -p intel # Available paritions: intel, batch, highmem, gpu, short (each has walltime and memory limits)

# Print current date
date

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd /rhome/danielaz/bigdata/transcriptomics/chr3

module load bwa-mem2

REFERENCE=/rhome/danielaz/bigdata/polyergus/fselysi/f_selysi_v02.fasta
DE_NOVO_FA=/rhome/danielaz/bigdata/transcriptomics/trinity/neoclara.deNovo.fa 

bwa-mem2 mem -t 8 ${REFERENCE} ${DE_NOVO_FA} > neoclara.de.novo.sam

# Print name of node
hostname
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

sbatch alignment.sh

# ouputs: 
neoclara.de.novo.sam
__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

# Extract scaffold 3 from the sam file: 

python3 chrExtractor.py --input neoclara.de.novo.sam

wc -l neoclara.de.novo.sam.scaffExtract 

# 16,306 lines out of total ~325,381 transcripts 

# Compare the list of transcripts mapping to Scaffold 3 with the list of DE isoform transcripts that are 4-fold DE'd. (isoforms)
cp /rhome/danielaz/bigdata/transcriptomics/raw_fastq/DESeq2_trans_v2/diffExpr.P1e-3_C2.matrix . 
python3  UniqScaffMatchr.py  --input1 diffExpr.P1e-3_C2.matrix --input2 neoclara.de.novo.sam.scaffExtract 

neoclara.de.novo.sam.scaffExtract.scaffsMatched
diffExpr.P1e-3_C2.matrix.scaffsMatched

# Compare the list of transcripts mapping to Scaffold 3 with the list of DE gene transcripts that are 4-fold DE'd. (genes)
cd /rhome/danielaz/bigdata/transcriptomics/raw_fastq/DESeq2_gene
python  --file1 diffExpr.P1e-3_C2.matrix --file2 neoclara.de.novo.sam.scaffExtract  

# 


__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.
## BELOW THIS LINE DEBUGGING NOTES 
squeue -u danielaz

