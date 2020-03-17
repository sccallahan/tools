#!/usr/bin/env Rscript
###############################################
# Author: Carson Callahan
# Purpose: combine list of csv files by column
# Date: 2020-03-16
# Notes:
###############################################

#### libraries ####
# install.packages('optparse')
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(plyr))

#### Make options ####
option_list = list(
  make_option(c("-p", "--path"), default = NA, type = "character",
              help = "full path to files"),
  make_option(c("-e", "--extension"), default = NA, type = "character",
              help = "extension for files of interest (e.g. '.txt')" ),
  make_option(c("-t", "--type"), default = NA, type = "character",
              help = "type of binding, one of {column, row}"),
  make_option(c("-f", "--filename"), default = NA, type = "character",
              help = "name of merged file to save WITHOUT extension (e.g. 'my_data'). Saved in {path}.")
)

opt_parser = OptionParser(usage = "Usage: %prog [-p path] [-e extension] [-t type]. Files should not have header.",
                          option_list = option_list)
opt = parse_args(opt_parser)

if (is.na(opt$p) | is.na(opt$e) | is.na(opt$t) | is.na(opt$f)){
  stop("Required argument missing. Use -h or --help for help.", call. = FALSE)
  print_help(opt_parser)
}


#### functions ####
# Taken from: https://stackoverflow.com/questions/7962267/cbind-a-dataframe-with-an-empty-dataframe-cbind-fill
cbind.fill <- function(...) {                                                                                                                                                       
  transpoted <- lapply(list(...),t)                                                                                                                                                 
  transpoted_dataframe <- lapply(transpoted, as.data.frame)                                                                                                                         
  return (data.frame(t(plyr::rbind.fill(transpoted_dataframe))))                                                                                                                          
} 

#### Read in the file list ####
# filenames <- list.files(opt$p)
# filenames_spec
filenames <- Sys.glob(file.path(opt$p, paste0("*",opt$e)))


#### combine files into one df ####
if (!opt$t %in% c("column", "row")){
  stop("'type' must be one of {column, row}")
}
all_csv <- lapply(filenames, function(i){read.csv(i, header = FALSE)})
if (opt$t == "column"){
  combined_df <- do.call(cbind.fill, all_csv)
  colnames(combined_df) <- filenames
}
if (opt$t == "row"){
  combined_df <- do.call(rbind.data.frame, all_csv)
}

#### save file ####
# print(paste0(opt$p, "/", opt$f, ".csv"))
write.csv(combined_df, file = paste0(opt$p, "/", opt$f, ".csv"), row.names = FALSE)
