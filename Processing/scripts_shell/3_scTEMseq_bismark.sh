#!/usr/bin/env bash
#
#PBS -l select=2:ncpus=4:mem=64GB 
#PBS -l walltime=72:00:00            
#PBS -j oe
#PBS -q xeon4q
#

# STEP 3: Bismark processing to .cov files
## This is where you need to change the directory to where your indexed genome assembly file is for bismark read alignment.
genomeDir="/home/bioinf/genomes/hs/bowtie2/old"

cd $PBS_O_WORKDIR


outDir="${PBS_O_WORKDIR}/outdir"

sample="${outDir}/1.scTEM_demux/*.fastq.gz"
ls -ll

mkdir -p ${outDir}

resFolder3="${outDir}/3.trim"
resFolder4="${outDir}/4.bismark"
resFolder5="${outDir}/5.methylation_extraction"

mkdir -p ${resFolder3}
mkdir -p ${resFolder4}
mkdir -p ${resFolder5}

# STEP 3A: Trimming with trimgalore
echo "==========================================================================="
echo "STEP 3A: Trimming raw reads"
date +"%Y%m%d.%H%M%S %A"

nThreads=2

module load trimgalore/0.6.5
trim_galore --output_dir ${resFolder3} --fastqc --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10 --paired ${sample}

date +"%Y%m%d.%H%M%S %A"
echo "====================================="

ls ${resFolder3}/

# STEP 3B: Read alignment with bismark
echo "==========================================================================="
echo "STEP 3B: Bismark read alignment"
date +"%Y%m%d.%H%M%S %A"

module load bismark

nThreads=4

module load bismark
bismark --version

for k in ${resFolder3}/*1_val_1.fq.gz
do
	sample=${k/1_val_1.fq.gz/}
	echo ${sample}
    bismark --genome_folder ${genomeDir} \
	--output_dir ${resFolder4} \
	--unmapped --ambiguous \
	-p ${nThreads}\
	--non_directional \
    -1 ${sample}1_val_1.fq.gz	\
    -2 ${sample}2_val_2.fq.gz 
    done

date +"%Y%m%d.%H%M%S %A"
echo "====================================="

ls ${resFolder3}/

# Step 3C: Methylation extrtaction with bismark
echo "==========================================================================="
echo "Step 3C: Methylation extraction"
date +"%Y%m%d.%H%M%S %A"

bismark_methylation_extractor -o ${resFolder5} --paired-end --parallel ${nThreads} --bedgraph ${resFolder4}/*.bam

date +"%Y%m%d.%H%M%S %A"
echo "====================================="

exit 0