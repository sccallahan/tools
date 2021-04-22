###################################################################################
# Author: Carson Callahan
#
# Date: 2021-04-21
#
# Purpose: Functions for collecting and scaling spike in PROseq samples
#
# Notes:
# 1) Uses the peppro pipeline for processing: http://peppro.databio.org/en/latest/
# 2) Requires using your spike in genome as a "prealignment"
#################################################################################

#' @title generate_pull_peppro_stats_command
#' 
#' @description Generates the rsync command needed to pull only the stats.tsv files from peppro outputs
#' 
#' @param username Your username on seadragon
#' @param path Top-level directory containing your peppro outputs.
#' @param savedir Where you'd like to save the files on your local disk
#' 
#' @note The resulting command can simply be copy/pasted into the terminal. You will then be prompted for your password.
#'
#' @return command in terminal output
#' @export
#'
generate_pull_peppro_stats_command <- function(username, path, savedir){
  
  # what is the parent folder?
  from <- paste0(username, "@seadragon.mdanderson.edu:", path)
  # where are we saving these files?
  to <- savedir
  # generalized rsync command for peppro output structure
  command <- paste0("rsync -avhPn --include=*/ --include='stats.tsv' --exclude=* --prune-empty-dirs ")
  # message
  message("Please copy/paste the following command into the terminal to collect your stats.tsv files!")
  # final command
  noquote(paste0(command, from, " ", to))
}

#' @title generate_pull_peppro_bams_command
#' 
#' @description Generates the rsync command needed to pull only the stats.tsv files from peppro outputs
#' 
#' @param genome Main genome used for aligments, entered as a string
#' @param username Your username on seadragon
#' @param path Top-level directory containing your peppro outputs.
#' @param savedir Where you'd like to save the files on your local disk
#' 
#' @note The resulting command can simply be copy/pasted into the terminal. You will then be prompted for your password.
#'
#' @return command in terminal output
#' @export
#'
generate_pull_peppro_bams_command <- function(genome, username, path, savedir){
  
  # what is the parent folder?
  from <- paste0(username, "@seadragon.mdanderson.edu:", path)
  # where are we saving these files?
  to <- savedir
  # generalized rsync command for peppro output structure
  command <- paste0("rsync -avhPn --include=*/ --include=*/aligned_", genome, "/* --exclude=* --prune-empty-dirs ")
  # message
  message("Please copy/paste the following command into the terminal to collect your bams!")
  # final command
  noquote(paste0(command, from, " ", to))
}

#' @title make_proseq_scale_factors
#' 
#' @description Generates the scaling factors based on number of spike-in reads aligned
#' 
#' @param path Top-level directory where the stats.tsv files are contained. Should ONLY contain the sample folders with the stats.tsv files.
#'
#' @note Spike-in reads are hard-coded to row 17 of the stats.tsv files. This can be changed if needed.
#' 
#' @return A data frame containing the sample name, number of spike in reads aligned, and the resulting scaling factors
#' @export
#'
make_proseq_scale_factors <- function(path){
  
  # grab all the folder names
  folders <- list.files(path = path, full.names = TRUE)
  
  # loop through and pull fly reads
  reads <- c()
  for (i in 1:length(folders)){
    file <- list.files(path = folders[i], full.names = TRUE)
    print(file)
    dat <- read.table(file, sep = '\t', header = FALSE)
    # this pulls row 17 (dm6 aligned reads) and column 2 (value) and may need to be changed
    dat <- dat[17, 2, drop = FALSE]
    rownames(dat) <- folders[i]
    dat <- rownames_to_column(dat, var = "sample")
    reads <- rbind(reads, dat)
  }
  
  # convert to numeric
  reads$V2 <- as.numeric(as.character(reads$V2))
  
  # calculate scale factor
  # based on https://www.activemotif.com/documents/1977.pdf
  reads$scale_factor <- min(reads$V2)/reads$V2
  
  # shorten sample names
  reads$sample_short <- gsub(pattern = paste0(path, "/(.*)"),
                             replacement = "\\1",
                             reads$sample)
  
  return(reads)
}




















