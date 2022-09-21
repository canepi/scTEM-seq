##########################################################################
# Purpose: Basic information, and script to obtain SINE-Alu annotations from the latest RepBase annotation.
# Output: text file containing SINE.Alu annotations (for use in scTEMseq_cov_summary.R function).
# Date: 02.Sep.22
# Version: 0.1.0
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# To do:
##########################################################################
# These two steps have been separated to allow users to more flexibly change where and how they want to download the annotation file. i.e. doesn't need to be done in R.
# However, if the annotation file is downloaded from a different source, the filtering function will not work for a different format/source and will need to modified accordingly.
# Future releases may enable a more flexible approach/options for filtering annotation files and sources.

library(data.table)
library(dplyr)


### Automate Download

Download.Annotation <- function(URL = "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.out.gz", 
                                annotation.dir = "./Annotations/"){
  
###### Create folder ########
  
  mkdirs <- function(fp) {
    if(!file.exists(fp)) {
        mkdirs(dirname(fp))
        dir.create(fp)
    }
  } 
  
  mkdirs(annotation.dir)
  
####################################

  options(timeout=10000) # Default in R is too short for most internet speeds to download the large file within the permitted time. Change this if it is still too short for your connection.

cat("Download annotation file from: ", URL)
message("\nThis may take a while depending on your connection. If R times out/takes too long on your computer, change the script 'timeout' option in the script (line 36). Or download it outside of R.")

  download.file(url=URL, destfile=paste0(annotation.dir,'hg38.fa.out.gz'), method='libcurl')

# Recommend noting the date/version you downloaded for future reference as the above points 

  Line1 = paste0("downloaded.file.name:\t",annotation.dir, basename(URL))
  Line2 = paste0("download.site:\t", URL)
  Line3 = paste0("download.date:\t", Sys.Date())
  
  save.download.info.name <- paste0(annotation.dir,basename(URL),"_download_info.txt")
  fileConn<-file(save.download.info.name)
  writeLines(c(Line1, Line2,Line3), fileConn)
  close(fileConn)

cat("Finished and saved download information in: ", save.download.info.name)
  
}





### Automate Filtering

Create.Filtered.Annotation <- function(anno.file = "./Annotations/hg38.fa.out.gz", 
                                       annotation.dir = "./Annotations/",
                                       filtered.anno.name = "SINE.Alu.anno.txt"){
  

###### Create folder ########
  
  mkdirs <- function(fp) {
    if(!file.exists(fp)) {
      mkdirs(dirname(fp))
      dir.create(fp)
    }
  } 
  
  mkdirs(annotation.dir)
  
####################################
  
# Read in anno file
cat("Reading in: ", anno.file)
  full_anno <- as.data.frame(fread(anno.file, skip = 3, col.names = c("SW_score", "perc_div", "perc_del", "perc_ins", "chr", "start", "end", "left", "strand?", "Name", "Class.Family", "Repeat.Pos.Begin", "Repeat.Pos.End","Repeat.Pos.Left","Sequential.ID.Num")))

# Filter for SINE-Alu elements only
cat("\nFiltering annotation file for 'SINE/Alu' positions")
  SINE.Alu <- full_anno %>% filter(Class.Family == "SINE/Alu") %>%
    select(Name,Class.Family, chr, start, end) 
  
# Save
cat("\nSaving filtered annotation file into: ", annotation.dir, filtered.anno.name)
  write.csv(SINE.Alu, file = paste0(annotation.dir, filtered.anno.name ), quote = F, row.names = F)
cat("\nFinished. The saved file can now be used with the scTEMseq cov summary script")  
  

}

