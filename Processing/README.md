scTEM-seq (methylation processing)
=========
The steps here are also explained in the [protocol]() and [method](https://www.nature.com/articles/s41598-022-09765-x) publications.

Content:
--------
* `/barcodes/`: Contains the example secondary indexes list that are used in 'STEP 1B: Demultiplexing secondary indexes'. Modify as required for your setup.
* `/scripts_shell/`: Example shell scripts to automate processing all scTEM-seq samples on a HPC and/or unix environment. In the future, if there is enough interest, we may provide a snakemake pipeline to integrate and automate processing and coverage summary analysis.

Processing:
--------
The commands used to generate the coverage files are outline below. This assumes you have unix environment, and have already created the genome index files using Bowtie 2. Hisat2 indexed genomes are also usable by bismark, but needs the additional '--hisat2' flag. If in doubt, check the very well written [Bismark documentation for setup and general usage](https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html).

Tools/versions used:
- bcl2fastq/2.20
- cutadapt/2.10-python3.6
- bcftools/1.10.2
- trimgalore/0.6.5
- fastqc/0.11.8
- bismark/0.22.3
N.B. In theory, newer versions should work fine with the scripts, but the main source of any potential clashes of versions is likely from **cutadapt**. In which case, you'd need to install the above version.


### STEP 1: Demultiplexing
#### STEP 1A: Initial demultiplexing
If you already have fastq files, skip to 'STEP 1B' (Add link). The sample sheet provided for Illumina sequencing should only include the sequencing indexes used. These are used to perform the initial demultiplexing of the pool(s) sequenced, which may or may not already be done for you by your sequencer or sequencing company. If you only have bcl files, you can use 'bcl2fastq' as shown below:

```
module load bcl2fastq

bcl2fastq --runfolder-dir /path/to/top_level/of/sequencing_run/ \
--sample-sheet /path/to/top_level/of/sequencing_run/SampleSheet.csv \
--output-dir /path/to/your/new_favourite/folder_name/fastq/ \
--create-fastq-for-index-reads
```
Output = An R1 and R2 fastq file per library pool sequenced.
Troubleshoot = If all pools end up in 'unaligned' double check the barcode sequences in the the 'SampleSheet.csv'. Often it's because you need to modify the second barcode to the 'reverse compliment' sequence.

#### STEP 1B: Demultiplex secondary indexes with cutadapt
This is where you will utilise the barcodes/indexes that were added to your samples before being pooled. Utilise or modify the files in '/barcodes/' that are in fasta format for this step.

```
cutadapt -e 0.15 --discard-untrimmed --no-indels \
-g file:barcodes_fwd.fasta -G file:barcodes_rev.fasta \
-o {name1}{name2}_input.1.fastq.gz -p {name1}{name2}_input.2.fastq.gz \
input_R1_001.fastq.gz input_R2_001.fastq.gz
```
Output = 'Secondary demultiplexed' fastq files. R1 and R2 for every sample that was sequenced with a secondary index.
Troubleshoot = If there isn't any (or way too few) R1 and R2 fastq files comapared to your samples, then either the initial demultiplexing was wrong, the barcodes in the fasta files are wrong, or sample prep error. If only a few sample samples are missing, and the rest of the steps work, then it's likely it's a sample prep issue and the sample was lost prior to pooling or didn't receive an index.

### STEP 2: QC check
After demultiplexing, General QC metrics including read length and basecall quality should be checked with software such as Fastqc. Read counts for all samples should be collated and compared to check indexing was executed as expected.

#### STEP 2A: Count read length/depth
This step essentially counts the number of reads (lines) in the demultiplexed fastq files

```
cd /path/to/secondary_demultiplexed_only/
echo "Counting reads in files demultiplexed using custom scTEM-seq indexes"
date +"%Y%m%d.%H%M%S %A"
sample="*.fastq.gz"

for k in ${sample}
do

	echo -n ${k}" " >> readcounts.txt; echo $(zcat ${k}|wc -l)/4|bc >> readcounts.txt

done

echo "Done - Didn't take too long!?"
date +"%Y%m%d.%H%M%S %A"

```
Output = A single text file with a read count (and sample name) per line.

#### STEP 2B: FastQC

```
module load fastqc
cd /path/to/secondary_demultiplexed_only/

fastqc /path/to/secondary_demultiplexed_only/*.fastq.gz
```
Output = fastqc file per sample. Use multi-qc to collate them into a single file.

### STEP 3: Bismark processing to .cov files
#### STEP 3A: Trimming with trimgalore
Trimming is performed with Trim Galore. For 150bp read lengths (used in Miseq kits), 10bp is removed from the 5’ and 3’ ends of both forward and reverse reads to ensure index sequences are removed.

```
trim_galore --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10 --paired input.file

```
Output = Trimmed 'fastq.gz' R1 and R2 files.

#### STEP 3B: Read alignment with Bismark
Despite the conserved TE sequences targeted, scTEM-seq reads effectively map to the genome. Genome mapping is performed using Bismark. Reads are aligned in paired end and non-directional mode. *Unique mapping rate for scTEM-seq reads should be >60%*
```
bismark --genome_folder /genome/folder/ --non_directional -1 input_1.fq.gz -2 input_2.fq.gz
```
Output = One '.bam' file per sample.

#### Step 3C: Methylation extrtaction with bismark
Bam files are then used for methylation extraction with the Bismark methylation extraction tool. Methylation extraction produces site specific methylation data that is easily readable for quantification of average methylation levels.

```
bismark_methylation_extractor --paired-end --bedgraph input.bam
```
Output = One coverage (.cov.gz) file per sample.

### Step 4: Review your QC files and size of your .cov files  



Output:
--------
This should produce individual coverage files per sample that can then be used to obtain 'global methylation levels' [following the steps in the parent directory](../). The size of the .cov file depends on the sequencing depth you've used.

You can use either only the samples clearly passing QC (from STEP 2) or all samples to obtain global methylation levels, which can also provide additional QC metrics if you use the SINE-Alu annotation file. If no SINE-Alu annotation is used, you can only use fastqc (STEP 2B) and read count (STEP 2A) files to assess QC.

### Automating for all samples
The subdirectory './scripts/shell' contains example scripts that we developed for use on the UoN HPC which uses 'slurm' submission. Adjust as required for your HPC or unix environment/setup. The scripts are setup to be run (qsub) directly from the folder containing the 'pooled' fastq file. For '3_scTEMseq_bismark.sh' you need to modify the location of 'genomeDir' (line 10) to the directory containing the bowtie2 indexed genome you've created. You can also use Hisat2 indexed genome if you add '--hisat2' flag to the bismark command. For more details check the very well written [Bismark documentation for setup and general usage](https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html).

Acknowledgements/Contact:
--------
Computational analysis: Dr Sean Burnard (sean.burnard@newcastle.edu.au)
Experimental protocol: Kooper Hunt, Dr Danielle Bond and Dr Heather Lee (heather.lee@newcastle.edu.au)  




TO DO:
- Finish generalising some of the code (step 1B) and directory paths for the example shell scripts.
	- make it clearer which sections need to be adjusted to run on other users HPC system.
- Add conda environment (yaml) to make it easier for users to install the same program versions.
- If there's enough interest, add automated snakemake pipeline to process sequence files and obtain global methylation estimates from cov files. This could be made to work locally and on HPC platforms.
