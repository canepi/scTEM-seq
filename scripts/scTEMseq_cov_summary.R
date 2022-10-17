##########################################################################
# Purpose: Function to evaluate scTEMseq coverage files in a given directory
# Output will be a table containing mean methylation levels and QC metrics per sample.
# Date: 02.Sep.22
# Version: 0.1.0
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# To do:
##########################################################################

#  Load packages
library(data.table)
library(GenomicRanges)
library(dplyr)
library(foreach)
#library(doParallel) # If analysing a lot of samples and running parallelised.

scTEMseq_cov_summary <- function(cov.dir, output.dir = cov.dir, cov.suffix = "cov.gz", SINE.Alu.anno = "NULL", run.parallel = FALSE){

if(output.dir == "") output.dir = cov.dir


run.parallel <- run.parallel
    
if(SINE.Alu.anno != "NULL"){
### Read in annotations of interest ###
  
  cat("Reading in annotation file: ", SINE.Alu.anno)

  SINEAluData=fread(SINE.Alu.anno)
  SINEAluGR=makeGRangesFromDataFrame(SINEAluData, start.field=c("start","begin"), end.field=c("end","stop"), seqnames.field = c("chr","Chromosome"), ignore.strand = T)
  seqlevelsStyle(SINEAluGR)<-"NCBI"
}
  
# List and read in all cov files

covfiles=list.files(cov.dir, pattern = cov.suffix, full.names = TRUE) %>% .[order(.)]
cat("\nFound", length(covfiles), "cov files to import.")
if(length(covfiles)>300) message("If analysing a lot of samples (>300), it's worth running in parallel. Either way, consider testing the time taken to run a subset of files.")

#covfiles=covfiles[1:10] ### To test the script on a small number of samples, unhash this line ###
Cov_summary <- NULL

if(SINE.Alu.anno == "NULL"){
  message("\nIt is highly recommend that you obtain the SINE-Alu annotation prior to running this to obtain all the relevant QC metrics.")
  message("\nFor now, this will still obtain the 'global methylation level' for each sample, without the additional QC info.")
  message("This will be run on a single core.")
  
  Cov_summary <- foreach(cov = covfiles, .inorder=FALSE, .combine = 'rbind', .export ='data.table') %do% { # To make this parallelised, change this to '%dopar%'
    library(data.table)
    library(dplyr)
    # Read in cov file
    tmp <- fread(cov, col.names = c("chromosome","start","end","methylation_percentage","count_methylated","count_unmethylated"))
    # setwd(cov.dir)
    
    tmp <- tmp %>% 
      filter(chromosome %in% c(1:22,"X", "Y", "y", "x")) %>%
      summarise(sample_id = basename(cov),
                Number_of_cytosines = length(methylation_percentage),
                Average_methylation = mean(methylation_percentage), # This is preferable, because sum of cytosines (methylated and unmethylated) gives too much weight to positions covered more than once.
                Percent_digital = (Number_of_cytosines - sum(methylation_percentage != 100 & methylation_percentage !=0))/ Number_of_cytosines *100,
                Number_of_duplicated_cytosines_in_coverage_file = nrow(.[duplicated(.[,1:3])]))
                # Mean_methylation_across_SINE_Alu_only = mean(methylation_percentage[overlap@from]))
  }
  
} else if(run.parallel == FALSE){
  message("\nLooping the import and global methylation estimation on a single core. This may take some time... Enjoy a coffee or tea!")

  Cov_summary <- foreach(cov.name = covfiles, .inorder=FALSE, .combine = 'rbind', .export ='data.table') %do% { # To make this parallelised, change this to '%dopar%'
    library(data.table)
    library(dplyr)
    
    #setwd(cov.dir)
    
  # Read in cov file
    cov <- fread(cov.name, data.table = T, col.names = c("chromosome","start","end","methylation_percentage","count_methylated","count_unmethylated"))
  # Find overlap
    GR1=makeGRangesFromDataFrame(cov, start.field="start", end.field="end", seqnames.field = "chromosome", ignore.strand = T)
    GR1=keepSeqlevels(GR1, value=intersect(seqlevels(GR1),seqlevels(SINEAluGR)), pruning.mode = "tidy")
    overlap=findOverlaps(GR1,SINEAluGR, minoverlap=1, maxgap= -1)

    tmp <- cov %>% 
      filter(chromosome %in% c(1:22,"X", "Y", "y", "x")) %>%
      summarise(sample_id = basename(cov.name),
                Number_of_cytosines = length(methylation_percentage),
                Average_methylation = mean(methylation_percentage), # This is preferable, because sum of cytosines (methylated and unmethylated) gives too much weight to positions covered more than once.
                Percent_digital = (Number_of_cytosines - sum(methylation_percentage != 100 & methylation_percentage !=0))/ Number_of_cytosines *100,
                Number_of_SINE_Alu_annotations = length(overlap@to[!duplicated(overlap@to)]),
                Number_of_cytosines_in_SINE_Alu_annotations = length(overlap@from[!duplicated(overlap@from)]),
                Percent_cytosines_in_SINE_Alu_annotations =Number_of_cytosines_in_SINE_Alu_annotations / Number_of_cytosines*100,
                Number_of_duplicated_cytosines_in_coverage_file = nrow(.[duplicated(.[,1:3])]),
                Mean_methylation_across_SINE_Alu_only = mean(methylation_percentage[overlap@from])) 
    
  }
  
} else if(run.parallel == TRUE){
  library(doParallel)
  N_CORES <- parallel::detectCores() # Detect how many cores are available
  my.cl <- parallel::makeCluster((N_CORES-1), type = "PSOCK") # Make cluster using 'socket' which works across unix and windows! Using 1 core less than totally available to allow for other work to be done while waiting.
  doParallel::registerDoParallel(cl = my.cl) # Register parallel options
  
  Cov_summary <- foreach(cov.name = covfiles, .inorder=FALSE, .combine = 'rbind', .export ='data.table') %dopar% { # To make this parallelised, change this to '%dopar%'
    library(data.table)
    library(dplyr)
    library(GenomicRanges)
    
    # Read in cov file
    cov <- fread(cov.name, col.names = c("chromosome","start","end","methylation_percentage","count_methylated","count_unmethylated"))
    # Find overlap
    GR1=makeGRangesFromDataFrame(cov, start.field="start", end.field="end", seqnames.field = "chromosome", ignore.strand = T)
    GR1=keepSeqlevels(GR1, value=intersect(seqlevels(GR1),seqlevels(SINEAluGR)), pruning.mode = "tidy")
    overlap=findOverlaps(GR1,SINEAluGR, minoverlap=1, maxgap= -1)
    
    tmp <- cov %>% 
      filter(chromosome %in% c(1:22,"X", "Y", "y", "x")) %>%
      summarise(sample_id = basename(cov.name),
                Number_of_cytosines = length(methylation_percentage),
                Average_methylation = mean(methylation_percentage), # This is preferable, because sum of cytosines (methylated and unmethylated) gives too much weight to positions covered more than once.
                Percent_digital = (Number_of_cytosines - sum(methylation_percentage != 100 & methylation_percentage !=0))/ Number_of_cytosines *100,
                Number_of_SINE_Alu_annotations = length(overlap@to[!duplicated(overlap@to)]),
                Number_of_cytosines_in_SINE_Alu_annotations = length(overlap@from[!duplicated(overlap@from)]),
                Percent_cytosines_in_SINE_Alu_annotations =Number_of_cytosines_in_SINE_Alu_annotations / Number_of_cytosines*100,
                Number_of_duplicated_cytosines_in_coverage_file = nrow(.[duplicated(.[,1:3])]),
                Mean_methylation_across_SINE_Alu_only = mean(methylation_percentage[overlap@from]))
    
  }
}


# Save 
message("\nFinished processing ", length(covfiles), " cov files.")
cat("\nSaving results to: ", output.dir,"Coverage_Summary.txt")
write.table(Cov_summary, file = paste0(output.dir,"Coverage_Summary.txt"), sep ="\t", quote = F, row.names = F)

}
########################################################