#!/usr/bin/env bash
#
#PBS -l select=1:ncpus=2:mem=16GB 
#PBS -l walltime=5:00:00            
#PBS -j oe
#PBS -k oe
#PBS -q xeon4q
#


cd $PBS_O_WORKDIR

outDir="${PBS_O_WORKDIR}/outdir"
mkdir -p ${outDir}

resFolder2="${outDir}/2.QC"
mkdir -p ${resFolder2}

echo "Output folder: ${resFolder2}"
echo ">> Runing on " `pwd`

samples="${outDir}/1.scTEM_demux/*.fastq.gz"

module load fastqc
fastqc --outdir ${resFolder2} ${samples}


exit 0 
