#!/usr/bin/env bash
#
#PBS -l select=1:ncpus=2:mem=4GB 
#PBS -l walltime=1:00:00            
#PBS -j oe
#PBS -q xeon4q
#


###Script for the quantifying of scTEM-seq demultiplexed fastq files. 
#
##Run this script from the parent directory above 'outDir'. 
##Input should all fastq files in the output directory of the scTEM-seq demultiplex script (cutadapt_demux.sh). 
##Output will be a txt file with the fastq file name and read count for each sample. 
#
##Good for comparing read counts between samples and to the original pool fastq file. 
#

cd $PBS_O_WORKDIR

echo "Counting reads in files demultiplexed using custom scTEM-seq indexes"
date +"%Y%m%d.%H%M%S %A"

outDir="${PBS_O_WORKDIR}/outdir"

mkdir -p ${outDir}
resFolder2="${outDir}/2.QC"

mkdir -p ${resFolder2}

sample="${outDir}/1.scTEM_demux/*.fastq.gz"

echo "number of reads"
for k in ${sample}
do

	echo -n $(basename ${k}) " " >> readcounts.txt; echo $(zcat ${k}|wc -l)/4|bc >> readcounts.txt
	
	done	

mv readcounts.txt ${resFolder2}

echo "Done"
date +"%Y%m%d.%H%M%S %A"

exit 0 