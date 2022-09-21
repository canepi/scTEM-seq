#!/usr/bin/env bash
#
#PBS -l select=1:ncpus=1:mem=16GB
#PBS -l walltime=2:00:00
#PBS -j oe
#PBS -q xeon4q
#


#####Script for the demultplexing of pooled scTEM-seq samples using secondary custom indexes in the scTEM-seq Alu or Line1 oligo sets.
#
###Input set in stagein should be fastq.gz files for each scTEM-seq pool (ie. plate) that was indexed with a single NebNext index pair.
#
###Output will be fastq files including  the i5 index, well number, i7 index, pool name and read 1 or 2 in the name 
#(eg. "i501_A1_i701_211011_KG1a_Alu_pool_S1_L001.1.fastq.gz"). With standard barcode lists, The script will output files for each position
#on the index plate if reads are found with the correct index pairs. Check size of demultiplexed fastq files as a sanity check.
#
###Script assumes a full plate of samples has been indexed with the scTEM-seq index plate. For less samples, the script will scan all index pairs
#from the plate and output those present. Keep in mind the output will be named based on the well the index came from- not necessarily the same 
#well as the sample. Check size of demultiplexed fastq files as a sanity check (majority of reads should be assigned to the correct index pairs).
#
###Check notes in script before running
#


cd $PBS_O_WORKDIR

echo "Demultiplexing scTEM-seq fastq files with internal barcodes"
date +"%Y%m%d.%H%M%S %A"

outDir="${PBS_O_WORKDIR}/outdir"
tmpDir="${PBS_O_WORKDIR}/tmp"

mkdir -p ${outDir}

resFolder1="${outDir}/1.scTEM_demux"

mkdir -p ${resFolder1}

module load cutadapt

thisSample="*_R1_001.fastq.gz"
ls -ll

##### Copy barcode/index list files to working directory or set directory to barcode/index lists here. 
cp ./barcodes/barcodes_fwd.fasta ${PBS_O_WORKDIR}
cp ./barcodes/barcodes_rev.fasta ${PBS_O_WORKDIR}

for k in ${thisSample}
do
	##get base name for each input fastq file
	sample=${k/_R1_001.fastq.gz/}
	
	##Where the magic happens. See Cutadapt user guide demultiplexing section for info about this command.
	cutadapt \
	-e 0.15 --discard-untrimmed \
	-g file:barcodes_fwd.fasta \
	-G file:barcodes_rev.fasta \
	-o {name1}{name2}_${sample}.1.fastq.gz -p {name1}{name2}_${sample}.2.fastq.gz \
	${sample}_R1_001.fastq.gz ${sample}_R2_001.fastq.gz
	
	
	mv -v i5[0-1][0-9]_*_i7[0-1][0-9]_${sample}.[1-2].fastq.gz ${resFolder1}
	
	done


echo "Done"
date +"%Y%m%d.%H%M%S %A"

exit 0 
