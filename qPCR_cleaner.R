#!/usr/bin/env Rscript
#######################################
# Author: Carson Callahan
# Purpose: Clean qPCR data
# Date: 2019-01-15
# Notes: Updated 2019-08-30
#######################################


#### libraries ####
# install.packages('optparse')
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tools))

#### Make options ####
option_list = list(
  make_option(c("-f", "--file"), default = NA, type = "character",
              help = "File to be cleaned up"),
  make_option(c("-c", "--columns"), default = NA, type = "numeric",
              help = "Number of columns used"),
  make_option(c("-r", "--rows"), default = NA, type = "numeric",
              help = "Number of rows used")
)

opt_parser = OptionParser(usage = "Usage: %prog [-f filename] [-c column numbers] [-r row numbers]", option_list = option_list)
opt = parse_args(opt_parser)

if (is.na(opt$f) | is.na(opt$c) | is.na(opt$r)){
  stop("Required argument missing. Use -h or --help for help.", call. = FALSE)
  print_help(opt_parser)
}


#### Read in the data ####
## Make sure it's a csv
if (file_ext(opt$f) != "csv"){
  stop("File is not a csv!", call. = FALSE)
}

## Make sure the file exists
if (file.exists(opt$f)){
  raw <- read.csv(file = opt$f, header=T)
} else {
  stop("File does not exist!", call. = FALSE)
}

## Remove rows we don't need
raw <- raw%>%
  filter(!grepl("Reference", Threshold..dRn.))
## Check we loaded the file

## Filter for well number and letter
number <- c(1:opt$c) # Which well numbers (columns) were used
letter <- LETTERS[1:opt$r] # Which well rows were used
## Collect only the parts of the dataframe occupied with CTs
# paste argument basically collapses things to a grepl-friendly format of "x|y|z" 
clean <- raw %>% 
  filter(grepl(paste(letter, collapse = "|"), raw$Well)) %>% 
  filter(grepl(paste(number, collapse = "|"), Well))
if (opt$c <= 9){
  clean <- clean[str_length(clean$Well) == 2, ] # gets rid of numbers, e.g. "12" when you only want "2"
} 
clean <- clean[, c(1,4)] # collects just the Well ID and CT value
write.csv(x = clean, file = paste0(gsub(".csv", '', opt$f), "_processed.csv"), row.names = F)         



